#!/usr/bin/env venv/bin/python
"""Simple web interface to browse current and archived simulation runs."""

from flask import Flask, render_template_string, send_file, request, jsonify, session, Response, stream_with_context
from pathlib import Path
from datetime import datetime
import hashlib
import zipfile
import tempfile
import shutil
import secrets

app = Flask(__name__)
app.secret_key = secrets.token_hex(16)  # For session management

# Add custom Jinja2 filter for formatting numbers with commas
@app.template_filter('commafy')
def commafy_filter(value):
    """Format a number with comma separators."""
    if value is None:
        return 'N/A'
    try:
        return f"{int(value):,}"
    except (ValueError, TypeError):
        return value

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / 'data'
ARCHIVE_DIR = DATA_DIR  # Runs stored directly in data directory
EXPORT_EXTENSION = '.nanosim'

# Cache for completed plot HTML, keyed by (dirname, theme).
# Populated by /plot-stream; consumed by /plot.
_plot_cache: dict = {}

# Cache for completed exports, keyed by token.
# Populated by /export-stream; consumed by /export-download.
_export_cache: dict = {}

FAVICON_SVG = "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'><rect width='32' height='32' rx='6' fill='%23ab63fa'/><polyline points='8,26 8,6 24,26 24,6' fill='none' stroke='white' stroke-width='4' stroke-linecap='round' stroke-linejoin='round'/></svg>"
FAVICON_VERSION = hashlib.md5(FAVICON_SVG.encode()).hexdigest()[:8]


@app.route('/favicon.ico')
@app.route('/favicon.svg')
def favicon():
    """Serve an inline SVG favicon."""
    svg = FAVICON_SVG.replace('%23', '#')
    resp = Response(svg, mimetype='image/svg+xml')
    resp.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    resp.headers['Pragma'] = 'no-cache'
    resp.headers['Expires'] = '0'
    return resp


HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Simulation Browser</title>
    <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={{ favicon_version }}">
    <script>
        // Set theme immediately to prevent flash
        (function() {
            const savedTheme = localStorage.getItem('theme') || 'dark';
            document.documentElement.setAttribute('data-theme', savedTheme);
        })();
    </script>
    <style>
        :root {
            --bg-primary: #1e1e1e;
            --bg-secondary: #2d2d2d;
            --bg-hover: #3a3a3a;
            --text-primary: #e0e0e0;
            --text-secondary: #b0b0b0;
            --border-color: #404040;
            --shadow: rgba(0,0,0,0.3);
            --param-bg: #383838;
            --msg-info: #5dade2;
            --msg-success: #28a745;
            --msg-error: #dc3545;
            --msg-muted: #9aa0a6;
        }
        
        [data-theme="light"] {
            --bg-primary: #f5f5f5;
            --bg-secondary: #ffffff;
            --bg-hover: #fafafa;
            --text-primary: #333333;
            --text-secondary: #666666;
            --border-color: #eeeeee;
            --shadow: rgba(0,0,0,0.1);
            --param-bg: #f8f9fa;
            --msg-info: #0056b3;
            --msg-success: #1e7e34;
            --msg-error: #b02a37;
            --msg-muted: #6c757d;
        }
        
        body { 
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            max-width: 1200px; 
            margin: 40px auto; 
            padding: 0 20px;
            background: var(--bg-primary);
            color: var(--text-primary);
            transition: background-color 0.3s ease, color 0.3s ease;
        }
        h1 { 
            color: var(--text-primary);
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .header-controls {
            display: flex;
            gap: 10px;
            align-items: center;
        }
        .theme-toggle {
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            padding: 8px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 18px;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            height: 38px;
        }
        .theme-toggle:hover {
            background: var(--bg-hover);
        }
        .refresh-btn {
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            padding: 8px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 18px;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            height: 38px;
        }
        .refresh-btn:hover {
            background: var(--bg-hover);
        }
        .combine-btn {
            background: linear-gradient(90deg, #ab63fa 0%, #00b8d4 100%);
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 16px;
            font-weight: 600;
            display: flex;
            align-items: center;
            gap: 6px;
            min-height: 38px;
        }
        .combine-btn:hover {
            background: linear-gradient(90deg, #9b4fe8 0%, #0099b8 100%);
        }
        .archive-list {
            background: transparent;
            box-shadow: none;
        }
        .archive-item { 
            padding: 15px;
            margin-bottom: 12px;
            border-radius: 8px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 4px var(--shadow);
        }
        .archive-item:last-child {
            border-bottom: 1px solid var(--border-color);
        }
        .archive-item:hover { background: var(--bg-hover); }
        .timestamp { 
            font-size: 1.2em; 
            font-weight: bold; 
            color: var(--text-primary);
            margin-bottom: 8px;
        }
        .status { 
            display: inline-block;
            padding: 4px 12px;
            border-radius: 4px;
            font-weight: 600;
            font-size: 0.85em;
            margin-bottom: 8px;
            margin-right: 8px;
        }
        .status-completed { background: #28a745; color: white; }
        .status-running { background: #ff9800; color: white; }
        .status-interrupted { background: #ffc107; color: black; }
        .status-error { background: #dc3545; color: white; }
        .status-unknown { background: #e2e3e5; color: #383d41; }
        .combined-chip {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 4px;
            font-weight: 600;
            font-size: 0.85em;
            margin-bottom: 8px;
            background: linear-gradient(90deg, #ab63fa 0%, #00b8d4 100%);
            color: white;
        }
        .dynamics-chip {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 4px;
            font-weight: 600;
            font-size: 0.85em;
            margin-bottom: 8px;
        }
        .dynamics-chip.reversible {
            background: #ab63fa;
            color: white;
        }
        .dynamics-chip.irreversible {
            background: #00b8d4;
            color: black;
        }
        .details { 
            color: var(--text-secondary); 
            font-size: 0.9em;
            line-height: 1.6;
        }
        .details-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 8px;
            margin: 8px 0;
        }
        .param { 
            font-family: 'Courier New', monospace;
            background: var(--param-bg);
            padding: 4px 8px;
            border-radius: 4px;
        }
        .no-archives {
            text-align: center;
            padding: 60px 20px;
            color: var(--text-secondary);
        }
        a { color: #007bff; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .edit-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #007bff;
            border-radius: 3px;
            font-size: 1em;
            transition: all 0.2s;
            text-align: center;
        }
        .edit-link:hover {
            background: #007bff;
            color: white;
            text-decoration: none;
        }
        .view-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #17a2b8;
            border-radius: 3px;
            font-size: 1em;
            color: #17a2b8;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .view-link:hover {
            background: #17a2b8;
            color: white;
            text-decoration: none;
        }
        .notes-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #007bff;
            border-radius: 3px;
            font-size: 1em;
            color: #007bff;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .notes-link:hover {
            background: #007bff;
            color: white;
            text-decoration: none;
        }
        .plot-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #ab63fa;
            border-radius: 3px;
            font-size: 1em;
            color: #ab63fa;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .plot-link:hover {
            background: #ab63fa;
            color: white;
            text-decoration: none;
        }
        .delete-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #dc3545;
            border-radius: 3px;
            font-size: 1em;
            color: #dc3545;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .delete-link:hover {
            background: #dc3545;
            color: white;
            text-decoration: none;
        }
        
        .export-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #28a745;
            border-radius: 3px;
            font-size: 1em;
            color: #28a745;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .export-link:hover {
            background: #28a745;
            color: white;
            text-decoration: none;
        }
        
        .import-btn {
            background: #28a745;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 16px;
            font-weight: 600;
            transition: background 0.2s;
            gap: 6px;
            min-height: 38px;
        }
        
        .import-btn:hover {
            background: #218838;
        }
        
        .import-status {
            margin-top: 15px;
            font-size: 14px;
            text-align: center;
        }
        .import-status .msg {
            display: inline-block;
        }
        .import-status .msg-info {
            color: var(--msg-info);
        }
        .import-status .msg-success {
            color: var(--msg-success);
        }
        .import-status .msg-error {
            color: var(--msg-error);
        }
        .import-status .msg-muted {
            color: var(--msg-muted);
        }
        
        .drop-zone {
            border: 2px dashed var(--border-color);
            border-radius: 8px;
            padding: 40px 20px;
            text-align: center;
            background: var(--bg-primary);
            transition: all 0.3s ease;
            cursor: pointer;
            margin: 20px 0;
        }
        
        .drop-zone:hover {
            border-color: #28a745;
            background: var(--bg-hover);
        }
        
        .drop-zone.drag-over {
            border-color: #28a745;
            background: var(--bg-hover);
            border-style: solid;
        }
        
        .drop-zone input[type="file"] {
            display: none;
        }
        
        .drop-zone-text {
            color: var(--text-primary);
            font-size: 16px;
            margin-bottom: 8px;
        }
        
        .drop-zone-hint {
            color: var(--text-secondary);
            font-size: 14px;
        }
        
        .file-selected {
            color: #28a745;
            font-weight: 600;
            margin-top: 10px;
        }
        .hidden {
            display: none;
        }
        
        /* Toast notification */
        .toast {
            position: fixed;
            bottom: 20px;
            right: 20px;
            background: #155724;
            color: white;
            padding: 12px 24px;
            border-radius: 4px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
            z-index: 2000;
            animation: slideIn 0.3s ease-out;
        }
        @keyframes slideIn {
            from {
                transform: translateX(400px);
                opacity: 0;
            }
            to {
                transform: translateX(0);
                opacity: 1;
            }
        }
        
        /* Modal styles */
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.6);
        }
        .modal-content {
            background-color: var(--bg-secondary);
            color: var(--text-primary);
            margin: 10% auto;
            padding: 20px;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            width: 80%;
            max-width: 600px;
            max-height: 80vh;
            overflow-y: auto;
        }
        .modal-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        .modal-header h2 {
            margin: 0;
            font-size: 1.3em;
        }
        .close-btn {
            font-size: 28px;
            font-weight: bold;
            color: #aaa;
            cursor: pointer;
            border: none;
            background: none;
        }
        .close-btn:hover {
            color: #000;
        }
        .notes-textarea {
            width: 100%;
            min-height: 150px;
            padding: 10px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            font-family: inherit;
            font-size: 14px;
            resize: vertical;
            background: var(--bg-primary);
            color: var(--text-primary);
        }
        .modal-footer {
            margin-top: 15px;
            display: flex;
            gap: 10px;
            justify-content: flex-end;
        }
        .btn {
            padding: 8px 16px;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
        }
        .btn-primary {
            background: #007bff;
            color: white;
        }
        .btn-primary:hover {
            background: #0056b3;
        }
        .btn-secondary {
            background: #6c757d;
            color: white;
        }
        .btn-secondary:hover {
            background: #545b62;
        }
        .btn-success {
            background: #28a745;
            color: white;
        }
        .btn-success:hover {
            background: #218838;
        }
        .btn-success:disabled {
            background: #6c757d;
            cursor: not-allowed;
        }
        .combine-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin: 20px 0;
        }
        .sim-column {
            border: 2px solid var(--border-color);
            border-radius: 8px;
            padding: 15px;
            background: var(--bg-primary);
        }
        .sim-column:first-child {
            border-color: #ab63fa;
        }
        .sim-column:last-child {
            border-color: #00d9ff;
        }
        .sim-column h3 {
            margin-top: 0;
            margin-bottom: 10px;
            font-size: 1.1em;
            padding-bottom: 8px;
            border-bottom: 2px solid var(--border-color);
        }
        .sim-column:first-child h3 {
            color: #ab63fa;
            border-bottom-color: #ab63fa;
        }
        .sim-column:last-child h3 {
            color: #00d9ff;
            border-bottom-color: #00d9ff;
        }
        .sim-list {
            max-height: 400px;
            overflow-y: auto;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            background: var(--bg-secondary);
        }
        .sim-option {
            padding: 12px;
            border-bottom: 1px solid var(--border-color);
            cursor: pointer;
            transition: all 0.2s;
            background: var(--bg-secondary);
        }
        .sim-option:last-child {
            border-bottom: none;
        }
        .sim-option:hover:not(.disabled) {
            background: var(--bg-hover);
        }
        .rev-option.selected {
            background: #ab63fa;
            color: white;
        }
        .rev-option.selected:hover {
            background: #9a4eeb;
        }
        .irr-option.selected {
            background: #00d9ff;
            color: #1e1e1e;
        }
        .irr-option.selected:hover {
            background: #00c4e6;
        }
        .sim-option.disabled {
            opacity: 0.4;
            cursor: not-allowed;
        }
        .sim-option-title {
            font-weight: bold;
            margin-bottom: 4px;
        }
        .sim-option-details {
            font-size: 0.85em;
            opacity: 0.8;
        }
        .sim-option.selected .sim-option-details {
            opacity: 0.9;
        }
        .modal-content-wide {
            max-width: 800px;
        }
        .section-pad {
            padding: 20px;
        }
        .helper-text {
            color: var(--text-secondary);
            margin-bottom: 15px;
        }
        .center-row {
            text-align: center;
        }
        .empty-sim-msg {
            padding: 20px;
            text-align: center;
            color: var(--text-secondary);
        }
        .archive-chip-row {
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
            margin-bottom: 8px;
        }
        .mismatch-panel {
            margin-bottom: 12px;
            padding: 10px;
            background: var(--param-bg);
            border-radius: 4px;
        }
        .warning-banner {
            background: #ff9800;
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            margin-bottom: 8px;
            font-size: 0.9em;
        }
        .meta-block {
            margin-top: 8px;
        }
        .notes-preview {
            font-style: italic;
            white-space: pre-wrap;
            margin-top: 4px;
            overflow: hidden;
            line-height: 1.5;
            max-height: 4.5em;
            transition: max-height 0.3s ease;
        }
        .notes-toggle {
            font-size: 0.9em;
            display: none;
        }
        .action-buttons {
            margin-top: 12px;
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 10px;
        }
        .export-status {
            font-size: 1.05em;
            color: var(--text-secondary);
            margin-top: 6px;
            text-align: center;
        }
        .export-substatus {
            font-size: 0.9em;
            color: var(--text-secondary);
            margin-top: 6px;
            text-align: center;
            min-height: 1.2em;
        }
        .export-help {
            color: var(--text-secondary);
            margin-top: 12px;
            font-size: 0.85em;
            text-align: center;
        }
        .export-close-btn {
            display: none;
            margin-top: 12px;
            margin-left: auto;
            margin-right: auto;
        }
        .export-log-box {
            margin-top: 14px;
            max-height: 220px;
            overflow-y: auto;
            background: var(--bg-primary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 10px 12px;
            font-family: 'Menlo', 'Consolas', monospace;
            font-size: 12px;
            line-height: 1.6;
            text-align: left;
            color: var(--text-primary);
        }
        .export-log-box .log-line { margin: 0; white-space: pre-wrap; word-break: break-word; }
        .export-log-box .log-error { color: var(--msg-error); }
        .atom-spinner.export-spinner {
            width: 64px;
            height: 64px;
            margin: 12px auto;
            position: relative;
        }
        .atom-spinner.export-spinner .nucleus {
            position: absolute;
            top: 50%;
            left: 50%;
            width: 10px;
            height: 10px;
            background: #28a745;
            border-radius: 50%;
            transform: translate(-50%, -50%);
            box-shadow: 0 0 10px #28a745;
        }
        .atom-spinner.export-spinner .orbit {
            position: absolute;
            top: 50%;
            left: 50%;
            border-radius: 50%;
            opacity: 0.4;
        }
        .atom-spinner.export-spinner .orbit-1 {
            width: 32px;
            height: 32px;
            margin: -16px 0 0 -16px;
            border: 2px solid #28a745;
            animation: rotate 1.5s linear infinite;
        }
        .atom-spinner.export-spinner .orbit-2 {
            width: 48px;
            height: 48px;
            margin: -24px 0 0 -24px;
            border: 2px solid #3ecf5e;
            animation: rotate 2s linear infinite;
            animation-delay: -0.66s;
        }
        .atom-spinner.export-spinner .orbit-3 {
            width: 64px;
            height: 64px;
            margin: -32px 0 0 -32px;
            border: 2px solid #66e08a;
            animation: rotate 2.5s linear infinite;
        }
        .atom-spinner.export-spinner .electron {
            position: absolute;
            top: -4px;
            left: 50%;
            margin-left: -3px;
            width: 6px;
            height: 6px;
            border-radius: 50%;
            background: #28a745;
            box-shadow: 0 0 8px #28a745;
        }
        
        /* Responsive styles for mobile */
        @media (max-width: 768px) {
            body {
                margin: 20px auto;
                padding: 0 15px;
            }
            h1 {
                font-size: 1.5em;
                flex-direction: column;
                gap: 15px;
                align-items: flex-start;
            }
            .header-controls {
                width: 100%;
                justify-content: space-between;
                flex-wrap: wrap;
            }
            .import-btn {
                order: 1;
            }
            .combine-btn {
                order: 2;
            }
            .refresh-btn {
                order: 3;
                margin-left: auto;
            }
            #themeToggle {
                order: 4;
            }
            .timestamp {
                font-size: 1em;
            }
            .details-grid {
                grid-template-columns: 1fr;
            }
            .modal-content {
                width: 95%;
                margin: 5% auto;
                padding: 15px;
            }
            .combine-grid {
                grid-template-columns: 1fr;
            }
            .sim-list {
                max-height: 250px;
            }
            .action-buttons {
                flex-wrap: wrap;
                row-gap: 8px;
            }
            .action-buttons a,
            .action-buttons button {
                font-size: 0.8em;
                padding: 4px 10px;
            }
            .edit-link,
            .view-link,
            .notes-link,
            .plot-link,
            .delete-link,
            .export-link {
                font-size: 0.8em;
                padding: 4px 10px;
            }
        }
        
        @media (max-width: 480px) {
            body {
                padding: 0 10px;
            }
            h1 {
                font-size: 1.3em;
            }
            .header-controls {
                gap: 8px;
            }
            .header-controls button,
            .header-controls .combine-btn,
            .header-controls .import-btn {
                padding: 6px 12px;
                font-size: 14px;
            }
            .status, .combined-chip, .dynamics-chip {
                font-size: 0.75em;
                padding: 3px 8px;
            }
            .edit-link,
            .view-link,
            .notes-link,
            .plot-link,
            .delete-link,
            .export-link {
                font-size: 0.8em;
                padding: 4px 10px;
            }
        }
    </style>
    <script>
        let currentDirname = null;
        let currentExportDirname = null;
        let exportEventSource = null;
        
        // Theme management
        function initTheme() {
            const savedTheme = localStorage.getItem('theme') || 'dark';
            document.documentElement.setAttribute('data-theme', savedTheme);
            updateThemeIcon();
        }
        
        function toggleTheme() {
            const currentTheme = document.documentElement.getAttribute('data-theme');
            const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
            document.documentElement.setAttribute('data-theme', newTheme);
            localStorage.setItem('theme', newTheme);
            updateThemeIcon();
        }
        
        function updateThemeIcon() {
            const theme = document.documentElement.getAttribute('data-theme');
            const btn = document.getElementById('themeToggle');
            btn.textContent = theme === 'dark' ? '☀️' : '🌙';
            btn.title = theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';
        }
        
        // Initialize theme on page load
        document.addEventListener('DOMContentLoaded', initTheme);
        
        function addThemeToUrl(event, link) {
            event.preventDefault();
            const theme = document.documentElement.getAttribute('data-theme');
            window.location.href = link.href + '?theme=' + theme;
        }
        
        function openNotesModal(dirname, currentNotes) {
            currentDirname = dirname;
            document.getElementById('notesTextarea').value = currentNotes || '';
            document.getElementById('notesModal').style.display = 'block';
        }
        
        function openNotesModalFromLink(linkElement) {
            const dirname = linkElement.getAttribute('data-dirname');
            const notes = linkElement.getAttribute('data-notes');
            openNotesModal(dirname, notes);
        }
        
        function closeNotesModal() {
            document.getElementById('notesModal').style.display = 'none';
        }

        function openExportModal(dirname) {
            currentExportDirname = dirname;
            const modal = document.getElementById('exportModal');
            document.getElementById('exportStatus').textContent = 'Preparing export...';
            document.getElementById('exportProgress').textContent = '';
            document.getElementById('exportLog').innerHTML = '';
            document.getElementById('exportHelp').style.display = 'block';
            document.getElementById('exportCloseBtn').style.display = 'none';
            document.getElementById('exportSpinner').style.display = 'block';
            modal.style.display = 'block';
            startExportStream(dirname);
        }

        function closeExportModal() {
            if (exportEventSource) {
                exportEventSource.close();
                exportEventSource = null;
            }
            document.getElementById('exportModal').style.display = 'none';
        }

        function appendExportLog(text, isError) {
            const logBox = document.getElementById('exportLog');
            const p = document.createElement('p');
            p.className = 'log-line' + (isError ? ' log-error' : '');
            p.textContent = text;
            logBox.appendChild(p);
            logBox.scrollTop = logBox.scrollHeight;
            while (logBox.children.length > 300) logBox.removeChild(logBox.firstChild);
        }

        function startExportStream(dirname) {
            if (exportEventSource) {
                exportEventSource.close();
            }

            const statusEl = document.getElementById('exportStatus');
            const progressEl = document.getElementById('exportProgress');
            const helpTextEl = document.getElementById('exportHelp');
            const closeBtnEl = document.getElementById('exportCloseBtn');
            const spinnerEl = document.getElementById('exportSpinner');

            const theme = document.documentElement.getAttribute('data-theme') || localStorage.getItem('theme') || 'dark';
            const streamUrl = '/export-stream/' + dirname + '?theme=' + theme;
            exportEventSource = new EventSource(streamUrl);

            exportEventSource.onopen = function() {
                statusEl.textContent = 'Preparing export...';
                progressEl.textContent = '';
            };

            exportEventSource.onmessage = function(e) {
                const line = e.data;
                appendExportLog(line, false);
                if (line.includes('Found')) statusEl.textContent = 'Scanning files...';
                else if (line.includes('Packed')) statusEl.textContent = 'Hashing and packaging files...';
                else if (line.includes('Writing MANIFEST')) statusEl.textContent = 'Writing manifest...';
                else if (line.includes('Writing SIGNATURE')) statusEl.textContent = 'Signing export...';

                const packedMatch = line.match(/Packed\s+([\d,]+)\/([\d,]+)/);
                if (packedMatch) {
                    const done = parseInt(packedMatch[1].replace(/,/g, ''), 10);
                    const total = parseInt(packedMatch[2].replace(/,/g, ''), 10);
                    if (total > 0) {
                        const pct = Math.min(100, Math.round((done / total) * 100));
                        progressEl.textContent = `${pct}% complete (${done.toLocaleString()}/${total.toLocaleString()} files)`;
                    }
                }
            };

            exportEventSource.addEventListener('done', function(e) {
                exportEventSource.close();
                exportEventSource = null;

                const token = (e.data || '').trim();
                const downloadId = 'dl_' + Math.random().toString(36).slice(2);
                statusEl.textContent = 'Done! Starting download...';
                progressEl.textContent = '100% complete';
                appendExportLog('Export complete. Starting download now...', false);

                const iframe = document.createElement('iframe');
                iframe.style.display = 'none';
                iframe.src = '/export-download/' + encodeURIComponent(token) + '?download_id=' + encodeURIComponent(downloadId);
                document.body.appendChild(iframe);

                const cookieName = 'download_complete_' + downloadId + '=1';
                const startedAt = Date.now();
                const poll = setInterval(() => {
                    if (document.cookie.includes(cookieName)) {
                        clearInterval(poll);
                        statusEl.textContent = 'Download complete.';
                        appendExportLog('Download complete. You can close this dialog.', false);
                        helpTextEl.style.display = 'none';
                        closeBtnEl.style.display = 'inline-block';
                        spinnerEl.style.display = 'none';
                        document.cookie = 'download_complete_' + downloadId + '=; Max-Age=0; path=/';
                        return;
                    }
                    if (Date.now() - startedAt > 120000) {
                        clearInterval(poll);
                        statusEl.textContent = 'Download started (completion confirmation timed out).';
                        appendExportLog('Download likely started, but completion confirmation timed out.', false);
                        closeBtnEl.style.display = 'inline-block';
                    }
                }, 250);
            });

            exportEventSource.addEventListener('error', function(e) {
                exportEventSource.close();
                exportEventSource = null;
                statusEl.textContent = 'Export failed.';
                progressEl.textContent = '';
                appendExportLog('ERROR: ' + (e.data || 'export error'), true);
                closeBtnEl.style.display = 'inline-block';
            });

            exportEventSource.onerror = function() {
                statusEl.textContent = 'Connection lost. Check server logs.';
            };
        }
        
        // Import dialog functions
        function openImportModal() {
            document.getElementById('importModal').style.display = 'block';
            document.getElementById('importStatus').innerHTML = '';
            document.getElementById('fileSelected').style.display = 'none';
            setupDropZone();
        }
        
        function closeImportModal() {
            document.getElementById('importModal').style.display = 'none';
            document.getElementById('importFile').value = '';
            document.getElementById('importStatus').innerHTML = '';
            document.getElementById('fileSelected').style.display = 'none';
        }
        
        function setupDropZone() {
            const dropZone = document.getElementById('dropZone');
            
            dropZone.addEventListener('dragover', (e) => {
                e.preventDefault();
                e.stopPropagation();
                dropZone.classList.add('drag-over');
            });
            
            dropZone.addEventListener('dragleave', (e) => {
                e.preventDefault();
                e.stopPropagation();
                dropZone.classList.remove('drag-over');
            });
            
            dropZone.addEventListener('drop', (e) => {
                e.preventDefault();
                e.stopPropagation();
                dropZone.classList.remove('drag-over');
                
                const files = e.dataTransfer.files;
                if (files.length > 0) {
                    const file = files[0];
                    const lower = file.name.toLowerCase();
                    if (lower.endsWith('.zip') || lower.endsWith('.nanosim')) {
                        document.getElementById('importFile').files = files;
                        handleFileSelect({ target: { files: files } });
                    } else {
                        document.getElementById('importStatus').innerHTML = '<span class="msg msg-error">Please drop a .nanosim/.zip file or use the browse button to select a directory</span>';
                    }
                }
            });
        }
        
        function handleFileSelect(event) {
            const files = event.target.files;
            const fileSelectedDiv = document.getElementById('fileSelected');
            
            if (files.length === 0) {
                return;
            }
            
            // Check if it's a directory (multiple files)
            if (files.length > 1) {
                // Directory selected
                fileSelectedDiv.textContent = '✓ Selected directory: ' + files[0].webkitRelativePath.split('/')[0];
                fileSelectedDiv.style.display = 'block';
                document.getElementById('importStatus').innerHTML = '';
            } else {
                // Single file selected
                fileSelectedDiv.textContent = '✓ Selected: ' + files[0].name;
                fileSelectedDiv.style.display = 'block';
                document.getElementById('importStatus').innerHTML = '';
            }
        }
        
        // Combine dialog functions
        let selectedRev = null;
        let selectedIrr = null;
        
        function openCombineModal() {
            document.getElementById('combineModal').style.display = 'block';
            updateCombineButton();
        }
        
        function closeCombineModal() {
            document.getElementById('combineModal').style.display = 'none';
            selectedRev = null;
            selectedIrr = null;
            document.querySelectorAll('.sim-option.selected').forEach(el => el.classList.remove('selected'));
            // Reset all hidden states
            document.querySelectorAll('.sim-option').forEach(el => {
                el.style.display = 'block';
            });
            // Hide all no-match messages
            document.querySelectorAll('.no-match-message').forEach(msg => {
                msg.style.display = 'none';
            });
        }
        
        function selectRev(dirname) {
            const selectedEl = event.target.closest('.sim-option');
            
            // If clicking already selected item, deselect it
            if (selectedRev === dirname) {
                selectedRev = null;
                selectedEl.classList.remove('selected');
                
                // Show all irr options again
                document.querySelectorAll('.irr-option').forEach(el => {
                    el.style.display = 'block';
                });
                
                // Hide no-match message
                const irrList = document.getElementById('irrList');
                const noMatchMsg = irrList.querySelector('.no-match-message');
                if (noMatchMsg) {
                    noMatchMsg.style.display = 'none';
                }
                
                updateCombineButton();
                return;
            }
            
            selectedRev = dirname;
            document.querySelectorAll('.rev-option').forEach(el => el.classList.remove('selected'));
            selectedEl.classList.add('selected');
            
            // Get parameters from selected rev
            const revParams = {
                n: selectedEl.dataset.n,
                sweeps: selectedEl.dataset.sweeps,
                radius: selectedEl.dataset.radius,
                runs: selectedEl.dataset.runs
            };
            
            // Update irr options - hide incompatible ones
            let hasCompatible = false;
            document.querySelectorAll('.irr-option').forEach(el => {
                const irrParams = {
                    n: el.dataset.n,
                    sweeps: el.dataset.sweeps,
                    radius: el.dataset.radius,
                    runs: el.dataset.runs
                };
                
                const isCompatible = revParams.n === irrParams.n &&
                                   revParams.sweeps === irrParams.sweeps &&
                                   revParams.radius === irrParams.radius &&
                                   revParams.runs === irrParams.runs;
                
                if (isCompatible) {
                    el.style.display = 'block';
                    hasCompatible = true;
                } else {
                    el.style.display = 'none';
                    if (el.classList.contains('selected')) {
                        el.classList.remove('selected');
                        selectedIrr = null;
                    }
                }
            });
            
            // Show/hide no match message
            const irrList = document.getElementById('irrList');
            let noMatchMsg = irrList.querySelector('.no-match-message');
            if (!hasCompatible) {
                if (!noMatchMsg) {
                    noMatchMsg = document.createElement('div');
                    noMatchMsg.className = 'no-match-message';
                    noMatchMsg.style.cssText = 'padding: 20px; text-align: center; color: var(--text-secondary); font-style: italic;';
                    noMatchMsg.textContent = 'No irreversible simulations with matching parameters';
                    irrList.appendChild(noMatchMsg);
                }
                noMatchMsg.style.display = 'block';
            } else if (noMatchMsg) {
                noMatchMsg.style.display = 'none';
            }
            
            updateCombineButton();
        }
        
        function selectIrr(dirname) {
            const selectedEl = event.target.closest('.sim-option');
            
            // If clicking already selected item, deselect it
            if (selectedIrr === dirname) {
                selectedIrr = null;
                selectedEl.classList.remove('selected');
                
                // Show all rev options again
                document.querySelectorAll('.rev-option').forEach(el => {
                    el.style.display = 'block';
                });
                
                // Hide no-match message
                const revList = document.getElementById('revList');
                const noMatchMsg = revList.querySelector('.no-match-message');
                if (noMatchMsg) {
                    noMatchMsg.style.display = 'none';
                }
                
                updateCombineButton();
                return;
            }
            
            selectedIrr = dirname;
            document.querySelectorAll('.irr-option').forEach(el => el.classList.remove('selected'));
            selectedEl.classList.add('selected');
            
            // Get parameters from selected irr
            const irrParams = {
                n: selectedEl.dataset.n,
                sweeps: selectedEl.dataset.sweeps,
                radius: selectedEl.dataset.radius,
                runs: selectedEl.dataset.runs
            };
            
            // Update rev options - hide incompatible ones
            let hasCompatible = false;
            document.querySelectorAll('.rev-option').forEach(el => {
                const revParams = {
                    n: el.dataset.n,
                    sweeps: el.dataset.sweeps,
                    radius: el.dataset.radius,
                    runs: el.dataset.runs
                };
                
                const isCompatible = revParams.n === irrParams.n &&
                                   revParams.sweeps === irrParams.sweeps &&
                                   revParams.radius === irrParams.radius &&
                                   revParams.runs === irrParams.runs;
                
                if (isCompatible) {
                    el.style.display = 'block';
                    hasCompatible = true;
                } else {
                    el.style.display = 'none';
                    if (el.classList.contains('selected')) {
                        el.classList.remove('selected');
                        selectedRev = null;
                    }
                }
            });
            
            // Show/hide no match message
            const revList = document.getElementById('revList');
            let noMatchMsg = revList.querySelector('.no-match-message');
            if (!hasCompatible) {
                if (!noMatchMsg) {
                    noMatchMsg = document.createElement('div');
                    noMatchMsg.className = 'no-match-message';
                    noMatchMsg.style.cssText = 'padding: 20px; text-align: center; color: var(--text-secondary); font-style: italic;';
                    noMatchMsg.textContent = 'No reversible simulations with matching parameters';
                    revList.appendChild(noMatchMsg);
                }
                noMatchMsg.style.display = 'block';
            } else if (noMatchMsg) {
                noMatchMsg.style.display = 'none';
            }
            
            updateCombineButton();
        }
        
        function updateCombineButton() {
            const btn = document.getElementById('doCombineBtn');
            btn.disabled = !(selectedRev && selectedIrr);
        }
        
        function doCombine() {
            if (!selectedRev || !selectedIrr) return;
            
            fetch('/combine-sims', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    rev: selectedRev,
                    irr: selectedIrr
                })
            })
            .then(response => response.json())
            .then(data => {
                if (data.success) {
                    sessionStorage.setItem('toastMessage', 'Simulations combined successfully: ' + data.archive_name);
                    window.location.reload();
                } else {
                    alert('Error combining simulations: ' + data.error);
                }
            })
            .catch(error => {
                alert('Error combining simulations: ' + error);
            });
        }
        
        function saveNotes() {
            const notes = document.getElementById('notesTextarea').value;
            
            fetch('/notes/' + currentDirname, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({notes: notes})
            })
            .then(response => response.json())
            .then(data => {
                if (data.success) {
                    closeNotesModal();
                    window.location.reload();
                } else {
                    alert('Error saving notes: ' + data.error);
                }
            })
            .catch(error => {
                alert('Error saving notes: ' + error);
            });
        }
        
        // Close modal when clicking outside
        window.onclick = function(event) {
            const notesModal = document.getElementById('notesModal');
            const importModal = document.getElementById('importModal');
            const combineModal = document.getElementById('combineModal');
            const exportModal = document.getElementById('exportModal');
            if (event.target == notesModal) closeNotesModal();
            if (event.target == importModal) closeImportModal();
            if (event.target == combineModal) closeCombineModal();
            if (event.target == exportModal) closeExportModal();
        }
        
        function showToast(message) {
            const toast = document.createElement('div');
            toast.className = 'toast';
            toast.textContent = message;
            document.body.appendChild(toast);
            setTimeout(() => {
                toast.remove();
            }, 3000);
        }
        
        // Show toast on page load if there's a pending message
        window.addEventListener('DOMContentLoaded', () => {
            const toastMessage = sessionStorage.getItem('toastMessage');
            if (toastMessage) {
                sessionStorage.removeItem('toastMessage');
                showToast(toastMessage);
            }
        });
        
        function toggleNotesExpand(dirname) {
            const notesDiv = document.getElementById('notesContent_' + dirname);
            const toggleBtn = document.getElementById('toggleNotes_' + dirname);
            
            if (notesDiv.style.maxHeight === 'none') {
                // Collapse
                notesDiv.style.maxHeight = '4.5em';
                toggleBtn.textContent = '▼ Show more';
            } else {
                // Expand
                notesDiv.style.maxHeight = 'none';
                toggleBtn.textContent = '▲ Show less';
            }
        }
        
        function deleteArchive(dirname) {
            if (confirm('Are you sure you want to delete this archived run? This cannot be undone.\\n\\nArchive: ' + dirname)) {
                fetch('/delete/' + dirname, {
                    method: 'POST'
                })
                .then(response => response.json())
                .then(data => {
                    if (data.success) {
                        sessionStorage.setItem('toastMessage', 'Archive deleted successfully');
                        window.location.reload();
                    } else {
                        alert('Error deleting archive: ' + data.error);
                    }
                })
                .catch(error => {
                    alert('Error deleting archive: ' + error);
                });
            }
        }
        
        function exportArchive(dirname) {
            openExportModal(dirname);
        }
        
        function importArchive(overwrite = false) {
            const fileInput = document.getElementById('importFile');
            const statusDiv = document.getElementById('importStatus');
            
            if (!fileInput.files || fileInput.files.length === 0) {
                statusDiv.innerHTML = '<span class="msg msg-error">Please select a file or directory</span>';
                return;
            }
            
            const formData = new FormData();
            
            // Check if it's a directory (multiple files) or single zip file
            if (fileInput.files.length > 1) {
                // Directory selected - send all files with their paths
                for (let i = 0; i < fileInput.files.length; i++) {
                    const file = fileInput.files[i];
                    // Use webkitRelativePath to preserve directory structure
                    const relativePath = file.webkitRelativePath || file.name;
                    // Create a new file with the relative path as the name
                    const blob = new Blob([file], { type: file.type });
                    formData.append('files', blob, relativePath);
                }
                formData.append('is_directory', 'true');
            } else {
                // Single zip file
                formData.append('file', fileInput.files[0]);
                formData.append('is_directory', 'false');
            }
            
            formData.append('overwrite', overwrite.toString());
            
            statusDiv.innerHTML = '<span class="msg msg-info">Importing and validating...</span>';
            
            fetch('/import', {
                method: 'POST',
                body: formData
            })
            .then(async response => {
                if (response.status === 409) {
                    // Conflict - duplicate exists
                    return response.json().then(data => {
                        if (confirm(`Archive "${data.archive_name}" already exists.\\n\\nDo you want to overwrite it? This cannot be undone.`)) {
                            // Retry with overwrite flag
                            importArchive(true);
                        } else {
                            statusDiv.innerHTML = '<span class="msg msg-muted">Import cancelled</span>';
                        }
                        return null; // Don't continue to .then
                    });
                }

                if (!response.ok) {
                    if (response.status === 413) {
                        throw new Error('Upload too large (HTTP 413). Increase nginx client_max_body_size and retry.');
                    }

                    const contentType = (response.headers.get('content-type') || '').toLowerCase();
                    if (contentType.includes('application/json')) {
                        const errData = await response.json();
                        throw new Error(errData.error || `HTTP ${response.status} ${response.statusText}`);
                    }

                    throw new Error(`HTTP ${response.status} ${response.statusText}`);
                }

                return response.json();
            })
            .then(data => {
                if (!data) return; // Cancelled or handled above
                
                if (data.success) {
                    statusDiv.innerHTML = '<span class="msg msg-success">✓ ' + data.message + '</span>';
                    setTimeout(() => {
                        closeImportModal();
                        window.location.reload();
                    }, 1500);
                } else {
                    let errorMsg = data.error;
                    if (data.details) {
                        errorMsg += '<br>' + data.details.join('<br>');
                    }
                    statusDiv.innerHTML = '<span class="msg msg-error">✗ ' + errorMsg + '</span>';
                }
            })
            .catch(error => {
                const msg = (error && error.message) ? error.message : String(error);
                if (msg.includes('Load failed') || msg.includes('Failed to fetch') || msg.includes('NetworkError')) {
                    statusDiv.innerHTML = '<span class="msg msg-error">✗ Import failed before reaching app (network/proxy). If using nginx, raise upload/time limits (client_max_body_size, proxy_read_timeout) and retry.</span>';
                } else {
                    statusDiv.innerHTML = '<span class="msg msg-error">✗ Import failed: ' + msg + '</span>';
                }
            });
        }
    </script>
</head>
<body>
    <!-- Notes Modal -->
    <div id="notesModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h2>Edit Notes</h2>
                <button class="close-btn" onclick="closeNotesModal()">&times;</button>
            </div>
            <textarea id="notesTextarea" class="notes-textarea" placeholder="Add notes about this simulation run..."></textarea>
            <div class="modal-footer">
                <button class="btn btn-secondary" onclick="closeNotesModal()">Cancel</button>
                <button class="btn btn-primary" onclick="saveNotes()">Save</button>
            </div>
        </div>
    </div>
    
    <div id="importModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h2>Import Simulation</h2>
                <button class="close-btn" onclick="closeImportModal()">&times;</button>
            </div>
            <div class="section-pad">
                <p class="helper-text">Import a previously exported simulation archive (.nanosim, .zip, or directory)</p>
                <div class="drop-zone" id="dropZone" onclick="document.getElementById('importFile').click()">
                    <input type="file" id="importFile" accept=".nanosim,.zip" onchange="handleFileSelect(event)" webkitdirectory directory multiple>
                    <div class="drop-zone-text">📦 Drop .nanosim/.zip here or click to browse</div>
                    <div class="drop-zone-hint">Drop: .nanosim or .zip • Browse: .nanosim, .zip, or directory</div>
                    <div id="fileSelected" class="file-selected hidden"></div>
                </div>
                <div class="center-row">
                    <button onclick="importArchive()" class="btn btn-primary">Import & Validate</button>
                </div>
                <div id="importStatus" class="import-status"></div>
            </div>
        </div>
    </div>
    
    <div id="combineModal" class="modal">
        <div class="modal-content modal-content-wide">
            <div class="modal-header">
                <h2>Combine Simulations</h2>
                <button class="close-btn" onclick="closeCombineModal()">&times;</button>
            </div>
            <p class="helper-text">Select one reversible and one irreversible simulation to combine into a new archive.</p>
            <div class="combine-grid">
                <div class="sim-column">
                    <h3>Reversible Simulations</h3>
                    <div class="sim-list" id="revList">
                        {% for archive in rev_archives %}
                        <div class="sim-option rev-option" onclick="selectRev('{{ archive.dirname }}')"
                             data-dirname="{{ archive.dirname }}"
                             data-n="{{ archive.params.n if archive.params else '' }}"
                             data-sweeps="{{ archive.params.sweeps if archive.params else '' }}"
                             data-radius="{{ archive.params.radius if archive.params else '' }}"
                             data-runs="{{ archive.params.runs if archive.params else '' }}">
                            <div class="sim-option-title">{{ archive.display_time }}</div>
                            <div class="sim-option-details">
                                {% if archive.params %}
                                n={{ archive.params.n|commafy }}, s={{ archive.params.sweeps|commafy }}, r={{ archive.params.radius|commafy }}
                                {% endif %}
                            </div>
                        </div>
                        {% endfor %}
                        {% if not rev_archives %}
                        <div class="empty-sim-msg">No completed reversible simulations found</div>
                        {% endif %}
                    </div>
                </div>
                <div class="sim-column">
                    <h3>Irreversible Simulations</h3>
                    <div class="sim-list" id="irrList">
                        {% for archive in irr_archives %}
                        <div class="sim-option irr-option" onclick="selectIrr('{{ archive.dirname }}')"
                             data-dirname="{{ archive.dirname }}"
                             data-n="{{ archive.params.n if archive.params else '' }}"
                             data-sweeps="{{ archive.params.sweeps if archive.params else '' }}"
                             data-radius="{{ archive.params.radius if archive.params else '' }}"
                             data-runs="{{ archive.params.runs if archive.params else '' }}">
                            <div class="sim-option-title">{{ archive.display_time }}</div>
                            <div class="sim-option-details">
                                {% if archive.params %}
                                n={{ archive.params.n|commafy }}, s={{ archive.params.sweeps|commafy }}, r={{ archive.params.radius|commafy }}
                                {% endif %}
                            </div>
                        </div>
                        {% endfor %}
                        {% if not irr_archives %}
                        <div class="empty-sim-msg">No completed irreversible simulations found</div>
                        {% endif %}
                    </div>
                </div>
            </div>
            <div class="modal-footer">
                <button class="btn btn-secondary" onclick="closeCombineModal()">Cancel</button>
                <button id="doCombineBtn" class="btn btn-success" onclick="doCombine()" disabled>Combine</button>
            </div>
        </div>
    </div>

    <div id="exportModal" class="modal">
        <div class="modal-content modal-content-wide">
            <div class="modal-header">
                <h2>Export Simulation</h2>
                <button class="close-btn" onclick="closeExportModal()">&times;</button>
            </div>
            <div class="section-pad">
                <div class="atom-spinner export-spinner" id="exportSpinner">
                    <div class="nucleus"></div>
                    <div class="orbit orbit-1"><div class="electron"></div></div>
                    <div class="orbit orbit-2"><div class="electron"></div></div>
                    <div class="orbit orbit-3"><div class="electron"></div></div>
                </div>
                <div id="exportStatus" class="export-status">Preparing export...</div>
                <div id="exportProgress" class="export-substatus"></div>
                <div id="exportLog" class="export-log-box"></div>
                <p id="exportHelp" class="export-help">Large datasets can take time to hash and package.</p>
                <button id="exportCloseBtn" class="btn btn-secondary export-close-btn" onclick="closeExportModal()">Close</button>
            </div>
        </div>
    </div>
    
    <h1>
        <span>📊 Simulation Browser</span>
        <div class="header-controls">
            <button class="import-btn" onclick="openImportModal()" title="Import a previously exported simulation archive">Import</button>
            <button class="combine-btn" onclick="openCombineModal()" title="Combine reversible and irreversible simulations into one archive">Combine</button>
            <button class="refresh-btn" onclick="window.location.reload()" title="Refresh the page to see latest simulations">🔄</button>
            <button id="themeToggle" class="theme-toggle" onclick="toggleTheme()" title="Switch between dark and light themes">🌙</button>
        </div>
    </h1>
    
    {% if archives %}
    <div class="archive-list">\n        {% for archive in archives %}
        <div class="archive-item">
            <div class="timestamp">{{ archive.display_time }}</div>
            <div class="archive-chip-row">
                {% if archive.is_combined %}
                <span class="combined-chip">COMBINED</span>
                {% elif archive.params %}
                <span class="dynamics-chip {{ 'irreversible' if archive.params.flag == 'i' else 'reversible' }}">{{ 'IRREVERSIBLE' if archive.params.flag == 'i' else 'REVERSIBLE' }}</span>
                {% endif %}
                <span class="status status-{{ archive.status_class }}">{{ archive.status }}</span>
            </div>
            <div class="details">
                {% if archive.is_combined %}
                {% if archive.params_mismatch %}
                <div class="mismatch-panel">
                    <div class="warning-banner">
                        ⚠️ Warning: Reversible and irreversible parameters don't match
                    </div>
                    {% if archive.rev_params and archive.irr_params %}
                    <div class="details-grid">
                        <div><strong>Reversible:</strong></div>
                        <div><span class="param">n={{ archive.rev_params.n|commafy }}, s={{ archive.rev_params.sweeps|commafy }}, r={{ archive.rev_params.radius|commafy }}, m={{ archive.rev_params.runs|commafy }}</span></div>
                        <div><strong>Irreversible:</strong></div>
                        <div><span class="param">n={{ archive.irr_params.n|commafy }}, s={{ archive.irr_params.sweeps|commafy }}, r={{ archive.irr_params.radius|commafy }}, m={{ archive.irr_params.runs|commafy }}</span></div>
                    </div>
                    {% endif %}
                </div>
                {% else %}
                {% if archive.rev_params %}
                <div class="details-grid">
                    <div><strong>Lattice:</strong> <span class="param">n={{ archive.rev_params.n|commafy }}</span></div>
                    <div><strong>Sweeps:</strong> <span class="param">s={{ archive.rev_params.sweeps|commafy }}</span></div>
                    <div><strong>Radius:</strong> <span class="param">r={{ archive.rev_params.radius|commafy }}</span></div>
                    <div><strong>Runs:</strong> <span class="param">m={{ archive.rev_params.runs|commafy }}</span></div>
                    <div><strong>Total sims:</strong> <span class="param">{{ archive.rev_params.total|commafy }}</span></div>
                </div>
                {% endif %}
                {% endif %}
                {% elif archive.params %}
                <div class="details-grid">
                    <div><strong>Lattice:</strong> <span class="param">n={{ archive.params.n|commafy }}</span></div>
                    <div><strong>Sweeps:</strong> <span class="param">s={{ archive.params.sweeps|commafy }}</span></div>
                    <div><strong>Radius:</strong> <span class="param">r={{ archive.params.radius|commafy }}</span></div>
                    <div><strong>Runs:</strong> <span class="param">m={{ archive.params.runs|commafy }}</span></div>
                    <div><strong>Total sims:</strong> <span class="param">{{ archive.params.total|commafy }}</span></div>
                </div>
                {% endif %}
                {% if archive.completion_info %}
                <div class="meta-block">
                    <strong>Runtime:</strong> {{ archive.completion_info.total_time }}<br>
                    <strong>Throughput:</strong> {{ archive.completion_info.throughput }}
                </div>
                {% endif %}
                {% if archive.progress %}
                <div class="meta-block">
                    <strong>Progress:</strong> {{ archive.progress }}
                </div>
                {% endif %}
                {% if archive.notes %}
                <div class="meta-block">
                    <strong>Notes:</strong> 
                    <a href="#" class="edit-link" data-dirname="{{ archive.dirname }}" data-notes="{{ archive.notes|escape }}" onclick="openNotesModalFromLink(this); return false;">Edit</a>
                    <div id="notesContent_{{ archive.dirname }}" class="notes-preview">{{ archive.notes }}</div>
                    <a href="#" id="toggleNotes_{{ archive.dirname }}" class="notes-toggle" onclick="toggleNotesExpand('{{ archive.dirname }}'); return false;">▼ Show more</a>
                </div>
                <script>
                    (function() {
                        const notesDiv = document.getElementById('notesContent_{{ archive.dirname }}');
                        const toggleBtn = document.getElementById('toggleNotes_{{ archive.dirname }}');
                        
                        // Check if content overflows (more than 3 lines)
                        if (notesDiv.scrollHeight > notesDiv.clientHeight) {
                            toggleBtn.style.display = 'inline';
                        }
                    })();
                </script>
                {% endif %}
                <div class="action-buttons">
                    <a href="/plot-loading/{{ archive.dirname }}" class="plot-link" onclick="addThemeToUrl(event, this)" title="View interactive plots and analysis">Plot data</a>
                    <a href="/view/{{ archive.dirname }}" class="view-link" title="Browse all files in this archive">View files</a>
                    {% if not archive.notes %}
                    <a href="#" class="notes-link" data-dirname="{{ archive.dirname }}" data-notes="" onclick="openNotesModalFromLink(this); return false;" title="Add notes about this simulation">Add notes</a>
                    {% endif %}
                    <a href="#" class="export-link" onclick="exportArchive('{{ archive.dirname }}'); return false;" title="Export this simulation as a validated zip file">Export</a>
                    <a href="#" class="delete-link" onclick="deleteArchive('{{ archive.dirname }}'); return false;" title="Permanently delete this archive">Delete</a>
                </div>
            </div>
        </div>
        {% endfor %}
    </div>
    {% else %}
    <div class="no-archives">
        <h2>No archived runs found</h2>
        <p>Archives will appear here after running simulations with <code>make sim</code></p>
    </div>
    {% endif %}
</body>
</html>
"""

FILE_LIST_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Archive Files - {{ dirname }}</title>
    <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={{ favicon_version }}">
    <script>
        // Set theme immediately to prevent flash
        (function() {
            const savedTheme = localStorage.getItem('theme') || 'dark';
            document.documentElement.setAttribute('data-theme', savedTheme);
        })();
    </script>
    <style>
        :root[data-theme="light"] {
            --bg-primary: #ffffff;
            --bg-secondary: #f8f9fa;
            --text-primary: #333333;
            --text-secondary: #666666;
            --border-color: #eeeeee;
            --link-color: #007bff;
        }
        
        :root[data-theme="dark"] {
            --bg-primary: #1a1a1a;
            --bg-secondary: #2a2a2a;
            --text-primary: #e0e0e0;
            --text-secondary: #999999;
            --border-color: #404040;
            --link-color: #4dabf7;
        }
        
        body { 
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            max-width: 1200px; 
            margin: 40px auto; 
            padding: 0 20px;
            background-color: var(--bg-primary);
            color: var(--text-primary);
            transition: background-color 0.3s, color 0.3s;
        }
        
        h1 { color: var(--text-primary); }
        
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
        }
        
        .theme-toggle {
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 18px;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            height: 38px;
        }
        
        .theme-toggle:hover {
            background: var(--bg-secondary);
        }
        
        .file-list { 
            list-style: none; 
            padding: 0; 
        }
        
        .file-item { 
            padding: 12px;
            border-bottom: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background-color 0.2s;
        }
        
        .file-item:hover { 
            background: var(--bg-secondary); 
        }
        
        .file-name { 
            font-family: 'Courier New', monospace;
            color: var(--text-primary);
        }
        
        .file-size { 
            color: var(--text-secondary); 
            font-size: 0.9em; 
        }
        
        a { 
            color: var(--link-color); 
            text-decoration: none; 
        }
        
        a:hover { 
            text-decoration: underline; 
        }
        
        .back-link { 
            margin-bottom: 20px; 
        }
        
        .back-link a {
            color: #007bff;
            font-size: 34px;
            font-weight: 800;
            line-height: 1;
            text-decoration: none;
            padding: 0px 14px 4px 14px;
            border-radius: 4px;
            align-items: center;
            justify-content: center;
            transition: all 0.2s;
        }
        
        .back-link a:hover {
            background: #007bff;
            color: white;
            text-decoration: none;
        }
        
        /* Responsive styles for mobile */
        @media (max-width: 768px) {
            body {
                padding: 0 15px;
                margin: 20px auto;
            }
            h1 {
                font-size: 1.5em;
            }
            .file-item {
                border: 1px solid var(--border-color);
                border-radius: 8px;
                margin-bottom: 12px;
                padding: 15px;
                background: var(--bg-secondary);
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                flex-direction: column;
                align-items: flex-start;
                gap: 8px;
            }
            .file-item:hover {
                box-shadow: 0 4px 8px rgba(0,0,0,0.15);
            }
            .file-name {
                word-break: break-all;
            }
            .file-size {
                font-size: 0.85em;
            }
        }
        
        @media (max-width: 480px) {
            body {
                padding: 0 10px;
                margin: 15px auto;
            }
            h1 {
                font-size: 1.3em;
            }
            .file-item {
                padding: 12px;
                margin-bottom: 10px;
            }
        }
        }
    </style>
</head>
<body>
    <div class="header">
        <div class="back-link"><a href="/" title="Back to simulation browser">‹</a></div>
        <button id="themeToggle" class="theme-toggle" onclick="toggleTheme()" title="Toggle theme">🌙</button>
    </div>
    <h1>📂 Archive: {{ dirname }}</h1>
    <ul class="file-list">
        {% for file in files %}
        <li class="file-item">
            <span class="file-name">{{ file.path }}</span>
            <span class="file-size">{{ file.size }}</span>
        </li>
        {% endfor %}
    </ul>
    
    <script>
        // Load theme from localStorage or default to dark
        const savedTheme = localStorage.getItem('theme') || 'dark';
        document.documentElement.setAttribute('data-theme', savedTheme);
        updateThemeIcon();
        
        function toggleTheme() {
            const currentTheme = document.documentElement.getAttribute('data-theme');
            const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
            document.documentElement.setAttribute('data-theme', newTheme);
            localStorage.setItem('theme', newTheme);
            updateThemeIcon();
        }
        
        function updateThemeIcon() {
            const theme = document.documentElement.getAttribute('data-theme');
            const btn = document.getElementById('themeToggle');
            btn.textContent = theme === 'dark' ? '☀️' : '🌙';
            btn.title = theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';
        }
    </script>
</body>
</html>
"""

def parse_status_file(archive_path):
    """Parse sim_status.txt to get status info."""
    status_file = find_status_file(archive_path, 'sim_status.txt')
    if not status_file:
        return None
    
    info = {}
    with open(status_file) as f:
        for line in f:
            if ':' in line:
                key, value = line.split(':', 1)
                info[key.strip()] = value.strip()
    return info

def parse_start_file(archive_path):
    """Parse sim_started.txt to get parameters."""
    start_file = find_status_file(archive_path, 'sim_started.txt')
    if not start_file:
        return None
    
    params = {}
    with open(start_file) as f:
        for line in f:
            if 'Parameters:' in line:
                # Parse: n=10, sweeps=100, flag=c, radius=11, runs=5
                param_str = line.split('Parameters:')[1].strip()
                for param in param_str.split(','):
                    key, value = param.strip().split('=')
                    params[key] = value
            elif 'Total simulations:' in line:
                params['total'] = line.split(':')[1].strip()
    return params

def parse_start_file_root_first(archive_path):
    """Parse sim_started.txt checking root first, then subdirectories."""
    start_file = find_status_file_root_first(archive_path, 'sim_started.txt')
    if not start_file:
        return None
    
    params = {}
    with open(start_file) as f:
        for line in f:
            if 'Parameters:' in line:
                # Parse: n=10, sweeps=100, flag=c, radius=11, runs=5
                param_str = line.split('Parameters:')[1].strip()
                for param in param_str.split(','):
                    key, value = param.strip().split('=')
                    params[key] = value
            elif 'Total simulations:' in line:
                params['total'] = line.split(':')[1].strip()
    return params

def parse_completion_file(archive_path):
    """Parse sim_completed.txt to get completion info."""
    completion_file = find_status_file(archive_path, 'sim_completed.txt')
    if not completion_file:
        return None
    
    info = {}
    with open(completion_file) as f:
        for line in f:
            if 'Total time:' in line:
                info['total_time'] = line.split(':', 1)[1].strip()
            elif 'Throughput:' in line:
                info['throughput'] = line.split(':', 1)[1].strip()
    return info

def find_status_file(archive_path, filename):
    """Find a status file in archive_path, checking rev/ and irr/ subdirectories."""
    # Check in order: irr/, rev/, then root
    for subdir in ['irr', 'rev', '']:
        check_path = archive_path / subdir / filename if subdir else archive_path / filename
        if check_path.exists():
            return check_path
    return None

def find_status_file_root_first(archive_path, filename):
    """Find a status file checking root first, then subdirectories."""
    # Check root first
    root_path = archive_path / filename
    if root_path.exists():
        return root_path
    
    # Then check subdirectories
    for subdir in ['rev', 'irr']:
        check_path = archive_path / subdir / filename
        if check_path.exists():
            return check_path
    
    return None

def read_notes(archive_path):
    """Read sim_notes.txt to get user notes."""
    # Check if this is a combined archive
    is_combined = (archive_path / 'rev').exists() and (archive_path / 'irr').exists()
    
    if is_combined:
        # For combined archives, read from root
        notes_file = archive_path / 'sim_notes.txt'
    else:
        # For rev/irr archives, check subdirectories first
        notes_file = find_status_file(archive_path, 'sim_notes.txt')
        if not notes_file:
            return None
    
    if not notes_file or not notes_file.exists():
        return None
    
    try:
        with open(notes_file, 'r') as f:
            return f.read().strip()
    except:
        return None

def format_size(size_bytes):
    """Format file size in human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"

@app.route('/')
def index():
    """Main page listing all archived runs."""
    archives = []
    
    # Add archived runs
    if ARCHIVE_DIR.exists():
        for archive_dir in sorted(ARCHIVE_DIR.iterdir(), reverse=True):
            if not archive_dir.is_dir():
                continue
            
            # Skip init_fin directory
            if archive_dir.name == 'init_fin':
                continue
            
            # Parse timestamp from directory name (YYYYMMDD_HHMMSS)
            dirname = archive_dir.name
            try:
                dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
                display_time = dt.strftime('%b %d, %Y %H:%M:%S')
            except ValueError:
                display_time = dirname
            
            # Get status
            status_info = parse_status_file(archive_dir)
            if status_info:
                status = status_info.get('Status', 'UNKNOWN')
                progress = status_info.get('Completed', None)
            else:
                # If no status file, check if simulation has started but not completed
                start_file = find_status_file(archive_dir, 'sim_started.txt')
                completion_file = find_status_file(archive_dir, 'sim_completed.txt')
                
                if start_file and not completion_file:
                    status = 'RUNNING'
                else:
                    status = 'UNKNOWN'
                progress = None
            
            # Map status to CSS class
            status_class = status.lower()
            
            # Check if this is a combined archive (has both rev and irr subdirectories)
            is_combined = (archive_dir / 'rev').exists() and (archive_dir / 'irr').exists()
            
            # Get parameters
            params = parse_start_file(archive_dir)
            rev_params = None
            irr_params = None
            
            if is_combined:
                # Check root first for combined archives
                root_params = parse_start_file_root_first(archive_dir)
                
                if root_params:
                    # Root has params - use them for both rev and irr
                    rev_params = root_params
                    irr_params = root_params
                else:
                    # Root doesn't have params - check subdirectories
                    rev_params = parse_start_file(archive_dir / 'rev')
                    irr_params = parse_start_file(archive_dir / 'irr')
                
                # Check if parameters match
                params_mismatch = False
                if rev_params and irr_params:
                    # Compare key parameters (excluding flag which is expected to differ)
                    for key in ['n', 'sweeps', 'radius', 'runs']:
                        if rev_params.get(key) != irr_params.get(key):
                            params_mismatch = True
                            break
            
            # Get completion info
            completion_info = parse_completion_file(archive_dir)
            
            # Get notes
            notes = read_notes(archive_dir)
            
            archives.append({
                'dirname': dirname,
                'display_time': display_time,
                'status': status,
                'status_class': status_class,
                'params': params,
                'progress': progress,
                'completion_info': completion_info,
                'notes': notes,
                'is_combined': is_combined,
                'rev_params': rev_params,
                'irr_params': irr_params,
                'params_mismatch': params_mismatch if is_combined else False
            })
    
    # Collect completed rev and irr archives for combine dialog
    rev_archives = []
    irr_archives = []
    
    # Add completed archives from archive directory
    if ARCHIVE_DIR.exists():
        for archive_dir in sorted(ARCHIVE_DIR.iterdir(), reverse=True):
            if not archive_dir.is_dir():
                continue
            
            # Skip init_fin directory
            if archive_dir.name == 'init_fin':
                continue
            
            # Skip already combined archives (those with both rev/ and irr/ subdirectories)
            is_combined = (archive_dir / 'rev').exists() and (archive_dir / 'irr').exists()
            if is_combined:
                continue
            
            # Get status - only include completed
            status_info = parse_status_file(archive_dir)
            if not status_info or status_info.get('Status') != 'COMPLETED':
                continue
            
            # Get parameters to determine if rev or irr
            params = parse_start_file(archive_dir)
            if not params:
                continue
            
            # Parse timestamp from directory name
            dirname = archive_dir.name
            try:
                dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
                display_time = dt.strftime('%b %d, %Y %H:%M:%S')
            except ValueError:
                display_time = dirname
            
            archive_info = {
                'dirname': dirname,
                'display_time': display_time,
                'params': params
            }
            
            # Determine if rev or irr - check flag parameter (supports both old 0/1 and new r/i values)
            flag = params.get('flag', 'r')
            # Map old numeric flags: 0 = reversible, 1 = irreversible
            if flag in ['r', 'rev', '0', 0]:
                rev_archives.append(archive_info)
            elif flag in ['i', 'irr', '1', 1]:
                irr_archives.append(archive_info)
    
    return render_template_string(HTML_TEMPLATE, archives=archives, rev_archives=rev_archives, irr_archives=irr_archives, favicon_version=FAVICON_VERSION)

@app.route('/notes/<dirname>', methods=['POST'])
def update_notes(dirname):
    """Update notes for a run."""
    from flask import jsonify, request
    
    archive_path = ARCHIVE_DIR / dirname
    if not archive_path.exists():
        return jsonify({'success': False, 'error': 'Archive not found'}), 404
    
    # Check if this is a combined archive
    is_combined = (archive_path / 'rev').exists() and (archive_path / 'irr').exists()
    
    if is_combined:
        # For combined archives, store in root
        notes_path = archive_path / 'sim_notes.txt'
    else:
        # For rev/irr archives, store in subdirectory
        if (archive_path / 'rev').exists():
            notes_path = archive_path / 'rev' / 'sim_notes.txt'
        elif (archive_path / 'irr').exists():
            notes_path = archive_path / 'irr' / 'sim_notes.txt'
        else:
            notes_path = archive_path / 'sim_notes.txt'
    
    try:
        data = request.get_json()
        notes = data.get('notes', '')
        
        # Write notes to file
        with open(notes_path, 'w') as f:
            f.write(notes)
        
        return jsonify({'success': True})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/combine-sims', methods=['POST'])
def combine_sims():
    """Combine a reversible and irreversible simulation into a new archive."""
    from flask import jsonify, request
    import shutil
    from datetime import datetime
    import time
    
    data = request.get_json()
    rev_dirname = data.get('rev')
    irr_dirname = data.get('irr')
    
    if not rev_dirname or not irr_dirname:
        return jsonify({'success': False, 'error': 'Both rev and irr must be specified'}), 400
    
    try:
        # Both should be in archives
        rev_path = ARCHIVE_DIR / rev_dirname
        irr_path = ARCHIVE_DIR / irr_dirname
        
        if not rev_path.exists() or not irr_path.exists():
            return jsonify({'success': False, 'error': 'One or both archives not found'}), 404
        
        # Create new combined archive with current timestamp (guaranteed unique after sleep)
        combined_timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        combined_path = ARCHIVE_DIR / combined_timestamp
        
        # Create the combined archive directory
        combined_path.mkdir(parents=True, exist_ok=True)
        
        # Copy rev data to combined/rev
        # Check if source has a rev subdirectory, otherwise copy the whole archive
        rev_source = rev_path / 'rev' if (rev_path / 'rev').exists() else rev_path
        rev_dest = combined_path / 'rev'
        shutil.copytree(rev_source, rev_dest)
        
        # Copy irr data to combined/irr
        # Check if source has an irr subdirectory, otherwise copy the whole archive
        irr_source = irr_path / 'irr' if (irr_path / 'irr').exists() else irr_path
        irr_dest = combined_path / 'irr'
        shutil.copytree(irr_source, irr_dest)
        
        # Create a note documenting the combination
        notes_file = combined_path / 'sim_notes.txt'
        with open(notes_file, 'w') as f:
            f.write(f"Combined simulation created {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Reversible: {rev_dirname}\n")
            f.write(f"Irreversible: {irr_dirname}\n")
            f.write(f"\n")
            
            # Read and append notes from rev archive
            rev_notes = read_notes(rev_path)
            if rev_notes:
                f.write(f"=== Reversible Notes ===\n")
                f.write(f"{rev_notes}\n")
                f.write(f"\n")
            
            # Read and append notes from irr archive
            irr_notes = read_notes(irr_path)
            if irr_notes:
                f.write(f"=== Irreversible Notes ===\n")
                f.write(f"{irr_notes}\n")
        
        return jsonify({'success': True, 'archive_name': combined_timestamp})
    except Exception as e:
        # Clean up on error
        if 'combined_path' in locals() and combined_path.exists():
            shutil.rmtree(combined_path)
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/delete/<dirname>', methods=['POST'])
def delete_archive(dirname):
    """Delete an archived run."""
    from flask import jsonify
    import shutil
    
    archive_path = ARCHIVE_DIR / dirname
    if not archive_path.exists():
        return jsonify({'success': False, 'error': 'Archive not found'}), 404
    
    if not archive_path.is_dir():
        return jsonify({'success': False, 'error': 'Not a directory'}), 400
    
    try:
        shutil.rmtree(archive_path)
        return jsonify({'success': True})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/view/<dirname>')
def view_files(dirname):
    """List files in an archive."""
    if dirname == 'current':
        archive_path = DATA_DIR
        if not archive_path.exists():
            return "Data directory not found", 404
    else:
        archive_path = ARCHIVE_DIR / dirname
        if not archive_path.exists():
            return "Archive not found", 404
    
    files = []
    for file_path in sorted(archive_path.rglob('*')):
        if file_path.is_file():
            rel_path = file_path.relative_to(archive_path)
            size = format_size(file_path.stat().st_size)
            files.append({'path': str(rel_path), 'size': size})
    
    return render_template_string(FILE_LIST_TEMPLATE, dirname=dirname, files=files, favicon_version=FAVICON_VERSION)

@app.route('/export-loading/<dirname>')
def export_loading(dirname):
    """Show loading page while an export zip is prepared."""
    from flask import request

    theme = request.args.get('theme', 'dark')

    if dirname == 'current':
        title = 'Current Run'
    else:
        try:
            dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
            title = dt.strftime('%b %d, %Y %H:%M:%S')
        except ValueError:
            title = dirname

    loading_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Preparing Export...</title>
        <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={FAVICON_VERSION}">
        <script>
            (function() {{
                const urlParams = new URLSearchParams(window.location.search);
                const urlTheme = urlParams.get('theme');
                const savedTheme = urlTheme || localStorage.getItem('theme') || 'dark';
                document.documentElement.setAttribute('data-theme', savedTheme);
            }})();
        </script>
        <style>
            :root {{
                --bg-primary: #1e1e1e;
                --bg-secondary: #2d2d2d;
                --text-primary: #e0e0e0;
                --text-secondary: #b0b0b0;
                --border-color: #404040;
                --log-bg: rgba(0,0,0,0.35);
                --log-border: rgba(255,255,255,0.12);
                --log-text: #c8e6c9;
                --log-error: #ff8a80;
            }}

            [data-theme="light"] {{
                --bg-primary: #ffffff;
                --bg-secondary: #f8f9fa;
                --text-primary: #333333;
                --text-secondary: #666666;
                --border-color: #dddddd;
                --log-bg: #f8f9fa;
                --log-border: #d8dee4;
                --log-text: #1f2937;
                --log-error: #b42318;
            }}

            body {{
                margin: 0;
                padding: 20px;
                font-family: Arial, sans-serif;
                background: var(--bg-primary);
                color: var(--text-primary);
                display: flex;
                flex-direction: column;
                align-items: center;
                justify-content: center;
                min-height: 100vh;
            }}

            .loading-container {{
                text-align: center;
                max-width: 700px;
            }}

            .atom-spinner {{
                position: relative;
                width: 80px;
                height: 80px;
                margin: 30px auto;
            }}

            .nucleus {{
                position: absolute;
                top: 50%;
                left: 50%;
                width: 12px;
                height: 12px;
                background: #28a745;
                border-radius: 50%;
                transform: translate(-50%, -50%);
                margin: 2px 0 0 2px;
                box-shadow: 0 0 10px #28a745, 0 0 20px #28a745;
            }}

            .orbit {{
                position: absolute;
                top: 50%;
                left: 50%;
                border-radius: 50%;
                opacity: 0.4;
            }}

            .orbit-1 {{
                width: 40px;
                height: 40px;
                margin: -20px 0 0 -20px;
                border: 2px solid #28a745;
                animation: rotate 1.5s linear infinite;
                box-shadow: 0 0 8px #28a745;
            }}

            .orbit-2 {{
                width: 60px;
                height: 60px;
                margin: -30px 0 0 -30px;
                border: 2px solid #3ecf5e;
                animation: rotate 2s linear infinite;
                animation-delay: -0.66s;
                box-shadow: 0 0 8px #3ecf5e;
            }}

            .orbit-3 {{
                width: 80px;
                height: 80px;
                margin: -40px 0 0 -40px;
                border: 2px solid #66e08a;
                animation: rotate 2.5s linear infinite;
                box-shadow: 0 0 8px #66e08a;
            }}

            .electron {{
                position: absolute;
                width: 6px;
                height: 6px;
                border-radius: 50%;
            }}

            .orbit-1 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
                background: #28a745;
                box-shadow: 0 0 5px #28a745, 0 0 10px #28a745;
            }}

            .orbit-2 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
                background: #3ecf5e;
                box-shadow: 0 0 5px #3ecf5e, 0 0 10px #3ecf5e;
            }}

            .orbit-3 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
                background: #66e08a;
                box-shadow: 0 0 5px #66e08a, 0 0 10px #66e08a;
            }}

            @keyframes rotate {{
                0% {{ transform: rotate(0deg); }}
                100% {{ transform: rotate(360deg); }}
            }}

            .status {{
                font-size: 1.2em;
                color: var(--text-secondary);
                margin-top: 20px;
            }}

            .sub-status {{
                font-size: 0.95em;
                color: var(--text-secondary);
                margin-top: 8px;
                min-height: 1.2em;
            }}

            .back-link {{
                position: absolute;
                top: 20px;
                left: 20px;
            }}

            .back-link a {{
                color: #007bff;
                font-size: 34px;
                font-weight: 800;
                line-height: 1;
                text-decoration: none;
                padding: 0px 14px 4px 14px;
                border-radius: 4px;
                transition: all 0.2s;
            }}

            .back-link a:hover {{
                background: #007bff;
                color: white;
            }}

            #log-box {{
                margin-top: 24px;
                width: 700px;
                max-width: 92vw;
                max-height: 260px;
                overflow-y: auto;
                background: var(--log-bg);
                border: 1px solid var(--log-border);
                border-radius: 6px;
                padding: 10px 14px;
                font-family: 'Menlo', 'Consolas', monospace;
                font-size: 12px;
                line-height: 1.6;
                color: var(--log-text);
                text-align: left;
            }}
            #log-box .log-line {{ margin: 0; white-space: pre-wrap; word-break: break-word; }}
            #log-box .log-error {{ color: var(--log-error); }}

            .help-text {{
                color: var(--text-secondary);
                margin-top: 16px;
                font-size: 0.85em;
            }}

            .done-close-btn {{
                display: none;
                margin-top: 16px;
                padding: 8px 16px;
                border-radius: 4px;
                border: 1px solid var(--border-color);
                background: var(--bg-secondary);
                color: var(--text-primary);
                cursor: pointer;
                font-size: 14px;
                font-weight: 600;
            }}

            .done-close-btn:hover {{
                opacity: 0.9;
            }}
        </style>
    </head>
    <body>
        <div class="back-link"><a href="/" title="Back to simulation browser">‹</a></div>
        <div class="loading-container">
            <h1>Exporting Archive - {title}</h1>
            <div class="atom-spinner" id="atom-spinner">
                <div class="nucleus"></div>
                <div class="orbit orbit-1"><div class="electron"></div></div>
                <div class="orbit orbit-2"><div class="electron"></div></div>
                <div class="orbit orbit-3"><div class="electron"></div></div>
            </div>
            <div class="status" id="status-text">Connecting...</div>
            <div class="sub-status" id="progress-text"></div>
            <div id="log-box"></div>
            <p class="help-text" id="help-text">Large datasets can take time to hash and package.</p>
            <button id="done-close-btn" class="done-close-btn" onclick="window.location.href='/'">Close</button>
        </div>
        <script>
            (function() {{
                const urlParams = new URLSearchParams(window.location.search);
                const theme = urlParams.get('theme') || localStorage.getItem('theme') || 'dark';
                const streamUrl = '/export-stream/{dirname}?theme=' + theme;

                const statusEl = document.getElementById('status-text');
                const progressEl = document.getElementById('progress-text');
                const logBox   = document.getElementById('log-box');
                const helpTextEl = document.getElementById('help-text');
                const doneCloseBtn = document.getElementById('done-close-btn');
                const atomSpinnerEl = document.getElementById('atom-spinner');

                function appendLog(text, isError) {{
                    const p = document.createElement('p');
                    p.className = 'log-line' + (isError ? ' log-error' : '');
                    p.textContent = text;
                    logBox.appendChild(p);
                    logBox.scrollTop = logBox.scrollHeight;
                    while (logBox.children.length > 300) logBox.removeChild(logBox.firstChild);
                }}

                const es = new EventSource(streamUrl);

                es.onopen = function() {{
                    statusEl.textContent = 'Preparing export...';
                    progressEl.textContent = '';
                }};

                es.onmessage = function(e) {{
                    const line = e.data;
                    appendLog(line, false);
                    if (line.includes('Found')) statusEl.textContent = 'Scanning files...';
                    else if (line.includes('Packed')) statusEl.textContent = 'Hashing and packaging files...';
                    else if (line.includes('Writing MANIFEST')) statusEl.textContent = 'Writing manifest...';
                    else if (line.includes('Writing SIGNATURE')) statusEl.textContent = 'Signing export...';

                    const packedMatch = line.match(/Packed\s+([\d,]+)\/([\d,]+)/);
                    if (packedMatch) {{
                        const done = parseInt(packedMatch[1].replace(/,/g, ''), 10);
                        const total = parseInt(packedMatch[2].replace(/,/g, ''), 10);
                        if (total > 0) {{
                            const pct = Math.min(100, Math.round((done / total) * 100));
                            progressEl.textContent = `${{pct}}% complete (${{done.toLocaleString()}}/${{total.toLocaleString()}} files)`;
                        }}
                    }}
                }};

                es.addEventListener('done', function(e) {{
                    es.close();
                    const token = (e.data || '').trim();
                    const downloadId = 'dl_' + Math.random().toString(36).slice(2);
                    statusEl.textContent = 'Done! Starting download...';
                    progressEl.textContent = '100% complete';
                    appendLog('Export complete. Starting download now...', false);

                    const iframe = document.createElement('iframe');
                    iframe.style.display = 'none';
                    iframe.src = '/export-download/' + encodeURIComponent(token) + '?download_id=' + encodeURIComponent(downloadId);
                    document.body.appendChild(iframe);

                    const cookieName = 'download_complete_' + downloadId + '=1';
                    const startedAt = Date.now();
                    const poll = setInterval(() => {{
                        if (document.cookie.includes(cookieName)) {{
                            clearInterval(poll);
                            statusEl.textContent = 'Download complete.';
                            appendLog('Download complete. You can return to the simulation browser.', false);
                            helpTextEl.style.display = 'none';
                            doneCloseBtn.style.display = 'inline-block';
                            atomSpinnerEl.style.display = 'none';
                            document.cookie = 'download_complete_' + downloadId + '=; Max-Age=0; path=/';
                            return;
                        }}
                        if (Date.now() - startedAt > 120000) {{
                            clearInterval(poll);
                            statusEl.textContent = 'Download started (completion confirmation timed out).';
                            appendLog('Download likely started, but completion confirmation timed out.', false);
                        }}
                    }}, 250);
                }});

                es.addEventListener('error', function(e) {{
                    es.close();
                    statusEl.textContent = 'Export failed.';
                    progressEl.textContent = '';
                    appendLog('ERROR: ' + (e.data || 'export error'), true);
                }});

                es.onerror = function() {{
                    statusEl.textContent = 'Connection lost. Check server logs.';
                    progressEl.textContent = '';
                }};
            }})();
        </script>
    </body>
    </html>
    """

    return loading_html


@app.route('/export-stream/<dirname>')
def export_stream(dirname):
    """SSE endpoint: builds export zip with progress updates and returns a download token."""

    if dirname == 'current':
        def current_not_allowed():
            yield "event: error\ndata: Cannot export current run\n\n"
        return Response(stream_with_context(current_not_allowed()),
                        mimetype='text/event-stream',
                        headers={'X-Accel-Buffering': 'no', 'Cache-Control': 'no-cache'})

    archive_path = ARCHIVE_DIR / dirname

    @stream_with_context
    def generate():
        if not archive_path.exists():
            yield "event: error\ndata: Archive not found\n\n"
            return

        token = secrets.token_urlsafe(24)
        with tempfile.NamedTemporaryFile(mode='w+b', suffix='.zip', delete=False) as tmp_zip:
            zip_path = Path(tmp_zip.name)

        try:
            yield "data: Starting export build...\n\n"

            file_list = sorted([p for p in archive_path.rglob('*') if p.is_file()])
            total_files = len(file_list)
            yield f"data: Found {total_files:,} files to export\n\n"

            manifest = {
                'export_version': '1.0',
                'archive_name': dirname,
                'export_date': datetime.now().isoformat(),
                'files': {}
            }

            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for index, file_path in enumerate(file_list, start=1):
                    rel_path = file_path.relative_to(archive_path)

                    sha256_hash = hashlib.sha256()
                    with open(file_path, 'rb') as f:
                        for chunk in iter(lambda: f.read(4096), b''):
                            sha256_hash.update(chunk)

                    manifest['files'][str(rel_path)] = {
                        'sha256': sha256_hash.hexdigest(),
                        'size': file_path.stat().st_size
                    }

                    zipf.write(file_path, arcname=str(rel_path))

                    if index == 1 or index == total_files or index % 25 == 0:
                        yield f"data: Packed {index:,}/{total_files:,}: {rel_path}\n\n"

                yield "data: Writing MANIFEST.json...\n\n"
                import json
                manifest_json = json.dumps(manifest, indent=2)
                zipf.writestr('MANIFEST.json', manifest_json)

                yield "data: Writing SIGNATURE...\n\n"
                manifest_hash = hashlib.sha256(manifest_json.encode()).hexdigest()
                zipf.writestr('SIGNATURE', manifest_hash)

            _export_cache[token] = {
                'path': zip_path,
                'filename': f'ember_nanosim_{dirname}{EXPORT_EXTENSION}'
            }

            yield "data: Export package complete\n\n"
            yield f"event: done\ndata: {token}\n\n"

        except Exception as e:
            if zip_path.exists():
                zip_path.unlink()
            yield f"event: error\ndata: Export failed: {str(e)}\n\n"

    return Response(generate(),
                    mimetype='text/event-stream',
                    headers={'X-Accel-Buffering': 'no', 'Cache-Control': 'no-cache'})


@app.route('/export-download/<token>')
def export_download(token):
    """Download a previously prepared export file and clean it up afterwards."""
    from flask import send_file, after_this_request, request

    entry = _export_cache.pop(token, None)
    if not entry:
        return "Export token expired or invalid", 404

    zip_path = entry['path']
    filename = entry['filename']

    if not zip_path.exists():
        return "Export file no longer available", 404

    @after_this_request
    def cleanup(response):
        try:
            if zip_path.exists():
                zip_path.unlink()
        except Exception:
            pass
        return response

    response = send_file(
        zip_path,
        as_attachment=True,
        download_name=filename,
        mimetype='application/octet-stream'
    )

    download_id = request.args.get('download_id')
    if download_id:
        response.set_cookie(
            f'download_complete_{download_id}',
            '1',
            max_age=120,
            path='/',
            samesite='Lax'
        )

    return response

@app.route('/export/<dirname>')
def export_archive(dirname):
    """Export an archive as a validated zip file."""
    from flask import send_file
    import json
    
    if dirname == 'current':
        return "Cannot export current run", 400
    
    archive_path = ARCHIVE_DIR / dirname
    if not archive_path.exists():
        return "Archive not found", 404
    
    # Create a temporary zip file
    with tempfile.NamedTemporaryFile(mode='w+b', suffix='.zip', delete=False) as tmp_zip:
        zip_path = Path(tmp_zip.name)
    
    try:
        # Create manifest with file hashes for validation
        manifest = {
            'export_version': '1.0',
            'archive_name': dirname,
            'export_date': datetime.now().isoformat(),
            'files': {}
        }
        
        # Create zip and calculate hashes
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in archive_path.rglob('*'):
                if file_path.is_file():
                    rel_path = file_path.relative_to(archive_path)
                    
                    # Calculate SHA256 hash of file
                    sha256_hash = hashlib.sha256()
                    with open(file_path, 'rb') as f:
                        for chunk in iter(lambda: f.read(4096), b''):
                            sha256_hash.update(chunk)
                    
                    # Store file and hash in manifest
                    manifest['files'][str(rel_path)] = {
                        'sha256': sha256_hash.hexdigest(),
                        'size': file_path.stat().st_size
                    }
                    
                    # Add file to zip
                    zipf.write(file_path, arcname=str(rel_path))
            
            # Add manifest to zip
            manifest_json = json.dumps(manifest, indent=2)
            zipf.writestr('MANIFEST.json', manifest_json)
            
            # Add manifest hash as signature
            manifest_hash = hashlib.sha256(manifest_json.encode()).hexdigest()
            zipf.writestr('SIGNATURE', manifest_hash)
        
        # Send the zip file
        return send_file(
            zip_path,
            as_attachment=True,
            download_name=f'ember_nanosim_{dirname}{EXPORT_EXTENSION}',
            mimetype='application/octet-stream'
        )
    except Exception as e:
        if zip_path.exists():
            zip_path.unlink()
        return f"Error creating export: {str(e)}", 500

@app.route('/import', methods=['POST'])
def import_archive():
    """Import and validate a simulation archive from zip file or directory."""
    import json
    
    is_directory = request.form.get('is_directory', 'false') == 'true'
    
    if is_directory:
        # Handle directory upload
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files provided'}), 400
        
        files = request.files.getlist('files')
        if not files:
            return jsonify({'success': False, 'error': 'No files in directory'}), 400
        
        # Create temporary directory to store files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Reconstruct directory structure
            for file in files:
                # Get relative path from webkitRelativePath
                relative_path = file.filename
                print(f"DEBUG: Uploading file with filename: {relative_path}")
                file_path = temp_path / relative_path
                file_path.parent.mkdir(parents=True, exist_ok=True)
                file.save(str(file_path))
            
            # Find the archive directory
            # When a directory is selected, files come with paths like "dirname/file.txt"
            # So we need to find the common root directory
            all_items = list(temp_path.iterdir())
            
            # If there's exactly one directory, that's our archive
            if len(all_items) == 1 and all_items[0].is_dir():
                archive_dir = all_items[0]
            # If there are multiple files/dirs at root level, use temp_path itself
            elif len(all_items) > 0:
                archive_dir = temp_path
            else:
                return jsonify({'success': False, 'error': 'No files found in upload'}), 400
            
            # Debug: Check what's actually in archive_dir
            print(f"DEBUG: archive_dir = {archive_dir}")
            print(f"DEBUG: archive_dir contents = {list(archive_dir.iterdir())[:10]}")  # First 10 items
            
            # Check for required validation files
            manifest_file = archive_dir / 'MANIFEST.json'
            signature_file = archive_dir / 'SIGNATURE'
            
            print(f"DEBUG: manifest exists = {manifest_file.exists()}, signature exists = {signature_file.exists()}")
            
            if not manifest_file.exists() or not signature_file.exists():
                return jsonify({'success': False, 'error': 'Missing validation files (MANIFEST.json or SIGNATURE)'}), 400
            
            try:
                # Read and verify manifest
                with open(manifest_file, 'r') as f:
                    manifest_json = f.read()
                    manifest = json.loads(manifest_json)
                
                with open(signature_file, 'r') as f:
                    expected_signature = f.read().strip()
                
                actual_signature = hashlib.sha256(manifest_json.encode()).hexdigest()
                
                if expected_signature != actual_signature:
                    return jsonify({'success': False, 'error': 'Manifest signature verification failed'}), 400
                
                # Verify all files match their hashes
                validation_errors = []
                for file_name, file_info in manifest['files'].items():
                    file_path = archive_dir / file_name
                    if not file_path.exists():
                        validation_errors.append(f"Missing file: {file_name}")
                        continue
                    
                    with open(file_path, 'rb') as f:
                        file_data = f.read()
                    
                    actual_hash = hashlib.sha256(file_data).hexdigest()
                    
                    if actual_hash != file_info['sha256']:
                        validation_errors.append(f"Hash mismatch: {file_name}")
                    
                    if len(file_data) != file_info['size']:
                        validation_errors.append(f"Size mismatch: {file_name}")
                
                if validation_errors:
                    return jsonify({
                        'success': False,
                        'error': 'File validation failed',
                        'details': validation_errors
                    }), 400
                
                # Determine archive name
                archive_name = manifest.get('archive_name', 'imported_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
                target_path = ARCHIVE_DIR / archive_name
                
                # Check for duplicate
                overwrite = request.form.get('overwrite', 'false') == 'true'
                
                if target_path.exists() and not overwrite:
                    return jsonify({
                        'success': False,
                        'conflict': True,
                        'archive_name': archive_name,
                        'message': f'Archive "{archive_name}" already exists.'
                    }), 409
                
                # If overwriting, delete existing
                if target_path.exists() and overwrite:
                    shutil.rmtree(target_path)
                
                # Copy files (excluding MANIFEST.json and SIGNATURE)
                target_path.mkdir(parents=True, exist_ok=True)
                for file_name in manifest['files'].keys():
                    src = archive_dir / file_name
                    dst = target_path / file_name
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dst)
                
                return jsonify({
                    'success': True,
                    'archive_name': archive_name,
                    'message': f'Successfully imported as {archive_name}'
                })
            
            except json.JSONDecodeError:
                return jsonify({'success': False, 'error': 'Invalid manifest format'}), 400
            except Exception as e:
                return jsonify({'success': False, 'error': f'Import failed: {str(e)}'}), 500
    
    else:
        # Handle zip file upload (original logic)
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'}), 400
        
        lower_name = file.filename.lower()
        if not (lower_name.endswith('.zip') or lower_name.endswith(EXPORT_EXTENSION)):
            return jsonify({'success': False, 'error': 'File must be a .nanosim or .zip'}), 400
        
        # Save uploaded file temporarily
        with tempfile.NamedTemporaryFile(mode='w+b', suffix='.zip', delete=False) as tmp_file:
            file.save(tmp_file.name)
            zip_path = Path(tmp_file.name)
        
        try:
            # Verify it's a valid zip
            if not zipfile.is_zipfile(zip_path):
                return jsonify({'success': False, 'error': 'Invalid zip file'}), 400
            
            with zipfile.ZipFile(zip_path, 'r') as zipf:
                # Check for required files - handle both root level and subdirectory
                zip_contents = zipf.namelist()
                
                # Find all MANIFEST.json files
                manifest_files = [f for f in zip_contents if f.endswith('MANIFEST.json') and not f.startswith('__MACOSX')]
                
                if len(manifest_files) == 0:
                    return jsonify({'success': False, 'error': 'Missing MANIFEST.json file'}), 400
                elif len(manifest_files) > 1:
                    return jsonify({'success': False, 'error': 'Multiple MANIFEST.json files found. Please ensure zip contains only one archive.'}), 400
                
                manifest_path = manifest_files[0]
                
                # Determine if files are in a subdirectory
                if '/' in manifest_path:
                    # Files are in subdirectory - get the prefix
                    subdir_prefix = manifest_path.rsplit('/', 1)[0] + '/'
                    signature_path = subdir_prefix + 'SIGNATURE'
                else:
                    # Files are at root level
                    subdir_prefix = ''
                    signature_path = 'SIGNATURE'
                
                # Check for SIGNATURE file
                if signature_path not in zip_contents:
                    return jsonify({'success': False, 'error': 'Missing SIGNATURE file'}), 400
                
                # Read and verify manifest signature
                manifest_json = zipf.read(manifest_path).decode('utf-8')
                manifest = json.loads(manifest_json)
                
                expected_signature = zipf.read(signature_path).decode('utf-8').strip()
                actual_signature = hashlib.sha256(manifest_json.encode()).hexdigest()
                
                if expected_signature != actual_signature:
                    return jsonify({'success': False, 'error': 'Manifest signature verification failed'}), 400
                
                # Verify all files match their hashes
                validation_errors = []
                for file_name, file_info in manifest['files'].items():
                    # Construct full path in zip (with subdirectory prefix if present)
                    zip_file_path = subdir_prefix + file_name
                    
                    if zip_file_path not in zip_contents:
                        validation_errors.append(f"Missing file: {file_name}")
                        continue
                    
                    # Calculate hash of file in zip
                    file_data = zipf.read(zip_file_path)
                    actual_hash = hashlib.sha256(file_data).hexdigest()
                    
                    if actual_hash != file_info['sha256']:
                        validation_errors.append(f"Hash mismatch: {file_name}")
                    
                    if len(file_data) != file_info['size']:
                        validation_errors.append(f"Size mismatch: {file_name}")
                
                if validation_errors:
                    return jsonify({
                        'success': False,
                        'error': 'File validation failed',
                        'details': validation_errors
                    }), 400
                
                # Determine archive name (use original or generate new if conflict)
                archive_name = manifest.get('archive_name', 'imported_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
                target_path = ARCHIVE_DIR / archive_name
                
                # Check for duplicate - if exists, return conflict status
                overwrite = request.form.get('overwrite', 'false') == 'true'
                
                if target_path.exists() and not overwrite:
                    return jsonify({
                        'success': False,
                        'conflict': True,
                        'archive_name': archive_name,
                        'message': f'Archive "{archive_name}" already exists.'
                    }), 409  # 409 Conflict
                
                # If overwriting, delete existing
                if target_path.exists() and overwrite:
                    shutil.rmtree(target_path)
                
                # Extract files (excluding MANIFEST.json and SIGNATURE)
                target_path.mkdir(parents=True, exist_ok=True)
                for file_name in manifest['files'].keys():
                    zip_file_path = subdir_prefix + file_name
                    file_data = zipf.read(zip_file_path)
                    file_path = target_path / file_name
                    file_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(file_path, 'wb') as f:
                        f.write(file_data)
                
                return jsonify({
                    'success': True,
                    'archive_name': archive_name,
                    'message': f'Successfully imported as {archive_name}'
                })
        
        except json.JSONDecodeError:
            return jsonify({'success': False, 'error': 'Invalid manifest format'}), 400
        except Exception as e:
            return jsonify({'success': False, 'error': f'Import failed: {str(e)}'}), 500
        finally:
            # Clean up temp file
            if zip_path.exists():
                zip_path.unlink()

@app.route('/plot-stream/<dirname>')
def plot_stream(dirname):
    """SSE endpoint: runs Sk_comparison.py, streams each stdout line to the browser,
    caches the finished plot HTML, then fires event: done so the loading page redirects."""
    import subprocess
    import tempfile
    import shutil

    theme = request.args.get('theme', 'dark')

    # If already cached (e.g. EventSource reconnect), immediately signal done.
    if (dirname, theme) in _plot_cache:
        def already_done():
            yield "event: done\ndata: ok\n\n"
        return Response(stream_with_context(already_done()),
                        mimetype='text/event-stream',
                        headers={'X-Accel-Buffering': 'no', 'Cache-Control': 'no-cache'})

    if dirname == 'current':
        data_path = DATA_DIR
    else:
        data_path = ARCHIVE_DIR / dirname

    plotly_template = 'plotly_dark' if theme == 'dark' else 'plotly_white'
    plot_script = REPO_ROOT / 'creutz-sim' / 'Sk_comparison.py'

    @stream_with_context
    def generate():
        if not data_path.exists():
            yield "data: ERROR: data path not found\n\n"
            yield "event: error\ndata: data path not found\n\n"
            return
        if not plot_script.exists():
            yield "data: ERROR: plotting script not found\n\n"
            yield "event: error\ndata: plotting script not found\n\n"
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            with open(plot_script, 'r') as f:
                script_content = f.read()

            # Inject the actual data path via args so no copytree is needed.
            # Pass an absolute path so Path(repo_root) / abs_path resolves to abs_path.
            modified_script = script_content.replace(
                "args = parser.parse_args()",
                f"args = parser.parse_args(['--data-dir', r'{data_path}'])"
            ).replace(
                'pio.templates.default = "plotly_white"',
                f'pio.templates.default = "{plotly_template}"'
            ).replace(
                "fig.show()",
                f"fig.write_html(r'{tmpdir_path}/plot1.html')"
            ).replace(
                "fig2.show()",
                f"fig2.write_html(r'{tmpdir_path}/plot2.html')"
            )

            temp_script = tmpdir_path / 'plot_script.py'
            with open(temp_script, 'w') as f:
                f.write(modified_script)

            python_bin = REPO_ROOT / 'venv' / 'bin' / 'python'
            proc = subprocess.Popen(
                [python_bin, str(temp_script)],
                cwd=tmpdir_path,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            for line in proc.stdout:
                line = line.rstrip()
                if line:
                    yield f"data: {line}\n\n"

            proc.wait()

            if proc.returncode != 0:
                yield "event: error\ndata: Script exited with non-zero status\n\n"
                return

            plot1 = tmpdir_path / 'plot1.html'
            plot2 = tmpdir_path / 'plot2.html'

            if not plot1.exists() and not plot2.exists():
                yield "event: error\ndata: No plot HTML files generated\n\n"
                return

            _plot_cache[(dirname, theme)] = {
                'plot1': plot1.read_text() if plot1.exists() else None,
                'plot2': plot2.read_text() if plot2.exists() else None,
            }

        yield "event: done\ndata: ok\n\n"

    return Response(generate(),
                    mimetype='text/event-stream',
                    headers={'X-Accel-Buffering': 'no', 'Cache-Control': 'no-cache'})


@app.route('/plot-loading/<dirname>')
def plot_loading(dirname):
    """Show loading page while plots are being generated."""
    from flask import request

    # Get theme from query parameter, default to dark
    theme = request.args.get('theme', 'dark')
    
    # Determine run type for spinner color
    if dirname == 'current':
        data_path = DATA_DIR
        title = 'Current Run'
    else:
        data_path = ARCHIVE_DIR / dirname
        try:
            dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
            title = dt.strftime('%b %d, %Y %H:%M:%S')
        except ValueError:
            title = dirname
    
    # Determine spinner color based on run type
    spinner_color = '#ab63fa'  # Default purple for reversible
    is_combined = False
    gradient_style = ""
    
    if data_path.exists():
        # Check if this is a combined archive
        is_combined = (data_path / 'rev').exists() and (data_path / 'irr').exists()
        
        if is_combined:
            # Use gradient for combined - apply to atom spinner
            spinner_color = '#ab63fa'
            gradient_style = """
            .nucleus {
                background: linear-gradient(135deg, #ab63fa 0%, #00b8d4 100%);
                box-shadow: 0 0 10px #ab63fa, 0 0 10px #00b8d4;
            }
            
            .orbit-1 {
                border-color: #ab63fa;
            }
            
            .orbit-2 {
                border-color: #8b7dd8;
            }
            
            .orbit-3 {
                border-color: #00b8d4;
            }
            
            .orbit-1 .electron {
                background: #ab63fa;
                box-shadow: 0 0 5px #ab63fa;
            }
            
            .orbit-2 .electron {
                background: #8b7dd8;
                box-shadow: 0 0 5px #8b7dd8;
            }
            
            .orbit-3 .electron {
                background: #00b8d4;
                box-shadow: 0 0 5px #00b8d4;
            }
            """
        else:
            # Check if irreversible
            params = parse_start_file(data_path)
            if params and params.get('flag') == 'i':
                spinner_color = '#00b8d4'  # Teal for irreversible
            else:
                spinner_color = '#ab63fa'  # Purple for reversible
    
    loading_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Generating Plots...</title>
        <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={FAVICON_VERSION}">
        <script>
            // Set theme immediately to prevent flash
            (function() {{
                const urlParams = new URLSearchParams(window.location.search);
                const urlTheme = urlParams.get('theme');
                const savedTheme = urlTheme || localStorage.getItem('theme') || 'dark';
                document.documentElement.setAttribute('data-theme', savedTheme);
            }})();
        </script>
        <style>
            :root {{
                --bg-primary: #1e1e1e;
                --bg-secondary: #2d2d2d;
                --text-primary: #e0e0e0;
                --text-secondary: #b0b0b0;
                --border-color: #404040;
                --log-bg: rgba(0,0,0,0.35);
                --log-border: rgba(255,255,255,0.12);
                --log-text: #c8e6c9;
                --log-error: #ff8a80;
            }}
            
            [data-theme="light"] {{
                --bg-primary: #ffffff;
                --bg-secondary: #f8f9fa;
                --text-primary: #333333;
                --text-secondary: #666666;
                --border-color: #dddddd;
                --log-bg: #f8f9fa;
                --log-border: #d8dee4;
                --log-text: #1f2937;
                --log-error: #b42318;
            }}
            
            body {{
                margin: 0;
                padding: 20px;
                font-family: Arial, sans-serif;
                background: var(--bg-primary);
                color: var(--text-primary);
                display: flex;
                flex-direction: column;
                align-items: center;
                justify-content: center;
                min-height: 100vh;
            }}
            
            .loading-container {{
                text-align: center;
                max-width: 600px;
            }}
            
            h1 {{
                margin-bottom: 30px;
                color: var(--text-primary);
            }}
            
            .atom-spinner {{
                position: relative;
                width: 80px;
                height: 80px;
                margin: 30px auto;
            }}
            
            .nucleus {{
                position: absolute;
                top: 50%;
                left: 50%;
                width: 12px;
                height: 12px;
                background: {spinner_color};
                border-radius: 50%;
                transform: translate(-50%, -50%);
                margin: 2px 0 0 2px;
                box-shadow: 0 0 10px {spinner_color}, 0 0 20px {spinner_color};
            }}
            
            .orbit {{
                position: absolute;
                top: 50%;
                left: 50%;
                border-radius: 50%;
                opacity: 0.4;
            }}
            
            .orbit-1 {{
                width: 40px;
                height: 40px;
                margin: -20px 0 0 -20px;
                border: 2px solid {spinner_color};
                animation: rotate 1.5s linear infinite;
                box-shadow: 0 0 8px {spinner_color};
            }}
            
            .orbit-2 {{
                width: 60px;
                height: 60px;
                margin: -30px 0 0 -30px;
                border: 2px solid {spinner_color};
                animation: rotate 2s linear infinite;
                animation-delay: -0.66s;
                box-shadow: 0 0 8px {spinner_color};
            }}
            
            .orbit-3 {{
                width: 80px;
                height: 80px;
                margin: -40px 0 0 -40px;
                border: 2px solid {spinner_color};
                animation: rotate 2.5s linear infinite;
                box-shadow: 0 0 8px {spinner_color};
            }}
            
            .electron {{
                position: absolute;
                width: 6px;
                height: 6px;
                background: {spinner_color};
                border-radius: 50%;
                box-shadow: 0 0 5px {spinner_color}, 0 0 10px {spinner_color};
            }}
            
            .orbit-1 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
            }}
            
            .orbit-2 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
            }}
            
            .orbit-3 .electron {{
                top: -4px;
                left: 50%;
                margin-left: -3px;
            }}
            
            @keyframes rotate {{
                0% {{ transform: rotate(0deg); }}
                100% {{ transform: rotate(360deg); }}
            }}
            
            {gradient_style}
            
            .status {{
                font-size: 1.2em;
                color: var(--text-secondary);
                margin-top: 20px;
            }}
            
            .back-link {{
                position: absolute;
                top: 20px;
                left: 20px;
            }}
            
            .back-link a {{
                color: #007bff;
                font-size: 34px;
                font-weight: 800;
                line-height: 1;
                text-decoration: none;
                padding: 0px 14px 4px 14px;
                border-radius: 4px;
                transition: all 0.2s;
            }}
            
            .back-link a:hover {{
                background: #007bff;
                color: white;
            }}
        </style>
        <style>
            #log-box {{
                margin-top: 24px;
                width: 640px;
                max-width: 92vw;
                max-height: 220px;
                overflow-y: auto;
                background: var(--log-bg);
                border: 1px solid var(--log-border);
                border-radius: 6px;
                padding: 10px 14px;
                font-family: 'Menlo', 'Consolas', monospace;
                font-size: 12px;
                line-height: 1.6;
                color: var(--log-text);
                text-align: left;
            }}
            #log-box .log-line {{ margin: 0; white-space: pre-wrap; word-break: break-all; }}
            #log-box .log-error {{ color: var(--log-error); }}
            .loading-help {{
                color: var(--text-secondary);
                margin-top: 16px;
                font-size: 0.85em;
            }}
        </style>
    </head>
    <body>
        <div class="back-link"><a href="/" title="Back to simulation browser">‹</a></div>
        <div class="loading-container">
            <h1>Simulation Plots - {title}</h1>
            <div class="atom-spinner">
                <div class="nucleus"></div>
                <div class="orbit orbit-1">
                    <div class="electron"></div>
                </div>
                <div class="orbit orbit-2">
                    <div class="electron"></div>
                </div>
                <div class="orbit orbit-3">
                    <div class="electron"></div>
                </div>
            </div>
            <div class="status" id="status-text">Connecting...</div>
            <div id="log-box"></div>
            <p class="loading-help">This may take a moment for larger datasets.</p>
        </div>
        <script>
            (function() {{
                const urlParams = new URLSearchParams(window.location.search);
                const theme = urlParams.get('theme') || localStorage.getItem('theme') || 'dark';
                const streamUrl = '/plot-stream/{dirname}?theme=' + theme;
                const destUrl   = '/plot/{dirname}?theme=' + theme;

                const statusEl = document.getElementById('status-text');
                const logBox   = document.getElementById('log-box');

                function appendLog(text, isError) {{
                    const p = document.createElement('p');
                    p.className = 'log-line' + (isError ? ' log-error' : '');
                    p.textContent = text;
                    logBox.appendChild(p);
                    logBox.scrollTop = logBox.scrollHeight;
                    // Keep last 200 lines to avoid DOM bloat
                    while (logBox.children.length > 200) logBox.removeChild(logBox.firstChild);
                }}

                const es = new EventSource(streamUrl);

                es.onopen = function() {{
                    statusEl.textContent = 'Reading CSV data...';
                }};

                es.onmessage = function(e) {{
                    const line = e.data;
                    appendLog(line, false);
                    // Update the headline status from key phrases in the output
                    if      (line.includes('Serializing fig2')) statusEl.textContent = 'Serializing summary plot...';
                    else if (line.includes('Serializing fig1')) statusEl.textContent = 'Serializing main plot...';
                    else if (line.includes('[rev]'))            statusEl.textContent = 'Processing reversible data...';
                    else if (line.includes('[irr]'))            statusEl.textContent = 'Processing irreversible data...';
                }};

                es.addEventListener('done', function() {{
                    es.close();
                    statusEl.textContent = 'Done! Preparing browser rendering...';
                    appendLog('Server finished. Browser is now loading and rendering plots...', false);
                    // Let the status/log paint before navigating to the heavy plot page.
                    requestAnimationFrame(() => window.location.replace(destUrl));
                }});

                es.addEventListener('error', function(e) {{
                    es.close();
                    statusEl.textContent = 'Error generating plots.';
                    appendLog('ERROR: ' + (e.data || 'stream error'), true);
                }});

                es.onerror = function() {{
                    // Connection dropped mid-stream — surface a message but don't redirect
                    statusEl.textContent = 'Connection lost. Check server logs.';
                }};
            }})();
        </script>
    </body>
    </html>
    """
    
    return loading_html

@app.route('/plot/<dirname>')
def plot_data(dirname):
    """Assemble and return the full plot page from the SSE-populated cache.
    If the cache is empty (direct navigation / refresh), redirect to the loading page
    which will re-run the stream."""
    from flask import request, redirect, url_for, make_response

    theme = request.args.get('theme', 'dark')
    cached = _plot_cache.get((dirname, theme))
    if not cached:
        # Not yet generated — send to loading page which drives /plot-stream
        return redirect(url_for('plot_loading', dirname=dirname, theme=theme))
    
    # Helper function to format numbers with commas
    def format_param(value):
        """Format a parameter value with commas if it's a number, otherwise return as-is."""
        if value is None or value == 'N/A':
            return 'N/A'
        try:
            return f"{int(value):,}"
        except (ValueError, TypeError):
            return value

    def status_chip_class(status_value):
        status_classes = {
            'COMPLETED': 'status-completed',
            'RUNNING': 'status-running',
            'INTERRUPTED': 'status-interrupted',
            'ERROR': 'status-error',
            'UNKNOWN': 'status-unknown'
        }
        return status_classes.get(status_value, 'status-unknown')
    
    # Pull cached plot HTML produced by /plot-stream
    plot1_html = cached.get('plot1')
    plot2_html = cached.get('plot2')

    if dirname == 'current':
        data_path = DATA_DIR
        if not data_path.exists():
            return "Data directory not found", 404
    else:
        data_path = ARCHIVE_DIR / dirname
        if not data_path.exists():
            return "Archive not found", 404

    try:
            # Combine plots into a single page
            if dirname == 'current':
                title = 'Current Run'
            else:
                # Format the title the same way as the cards
                try:
                    dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
                    title = dt.strftime('%b %d, %Y %H:%M:%S')
                except ValueError:
                    title = dirname
            
            # Check if this is a combined archive
            is_combined = (data_path / 'rev').exists() and (data_path / 'irr').exists()
            
            # Get status information
            status_info = parse_status_file(data_path)
            if status_info:
                status = status_info.get('Status', 'UNKNOWN')
            else:
                # If no status file, check if simulation has started but not completed
                start_file = find_status_file(data_path, 'sim_started.txt')
                completion_file = find_status_file(data_path, 'sim_completed.txt')
                
                if start_file and not completion_file:
                    status = 'RUNNING'
                else:
                    status = 'UNKNOWN'
            
            # Get parameters for this run
            params_html = ""
            if is_combined:
                # Check root first for combined archives
                root_params = parse_start_file_root_first(data_path)
                
                if root_params:
                    # Root has params - use them for both rev and irr
                    rev_params = root_params
                    irr_params = root_params
                else:
                    # Root doesn't have params - check subdirectories
                    rev_params = parse_start_file(data_path / 'rev')
                    irr_params = parse_start_file(data_path / 'irr')
                
                # Check if parameters match
                params_mismatch = False
                if rev_params and irr_params:
                    for key in ['n', 'sweeps', 'radius', 'runs']:
                        if rev_params.get(key) != irr_params.get(key):
                            params_mismatch = True
                            break
                
                # Only show params section if we actually have parameters
                if rev_params or irr_params:
                    params_html = '<div class="params-panel">'
                    
                    if params_mismatch:
                        params_html += '<div class="chip-combined">COMBINED</div>'
                        # Add status chip
                        status_class = status_chip_class(status)
                        params_html += f'<div class="chip chip-status chip-spaced {status_class}">{status}</div>'
                        params_html += '<div class="param-warning">⚠️ Warning: Reversible and irreversible parameters don\'t match</div>'
                        
                        # Show separate parameters when they don't match
                        if rev_params:
                            params_html += f"""
                            <div class="params-subpanel">
                                <strong>Reversible:</strong>
                                <strong>Lattice:</strong> n={format_param(rev_params.get('n', 'N/A'))} | 
                                <strong>Sweeps:</strong> s={format_param(rev_params.get('sweeps', 'N/A'))} | 
                                <strong>Radius:</strong> r={format_param(rev_params.get('radius', 'N/A'))} | 
                                <strong>Runs:</strong> m={format_param(rev_params.get('runs', 'N/A'))}
                            </div>
                            """
                        
                        if irr_params:
                            params_html += f"""
                            <div class="params-subpanel">
                                <strong>Irreversible:</strong>
                                <strong>Lattice:</strong> n={format_param(irr_params.get('n', 'N/A'))} | 
                                <strong>Sweeps:</strong> s={format_param(irr_params.get('sweeps', 'N/A'))} | 
                                <strong>Radius:</strong> r={format_param(irr_params.get('radius', 'N/A'))} | 
                                <strong>Runs:</strong> m={format_param(irr_params.get('runs', 'N/A'))}
                            </div>
                            """
                    else:
                        # Parameters match - show chip and params on same row
                        if rev_params:
                            # Build status chip HTML
                            status_class = status_chip_class(status)
                            status_chip_html = f'<div class="chip chip-status {status_class}">{status}</div>'
                            
                            params_html += f"""
                            <div class="params-row">
                                <div class="chip-combined">COMBINED</div>
                                {status_chip_html}
                                <div class="params-values">
                                    <strong>Lattice:</strong> n={format_param(rev_params.get('n', 'N/A'))} | 
                                    <strong>Sweeps:</strong> s={format_param(rev_params.get('sweeps', 'N/A'))} | 
                                    <strong>Radius:</strong> r={format_param(rev_params.get('radius', 'N/A'))} | 
                                    <strong>Runs:</strong> m={format_param(rev_params.get('runs', 'N/A'))}
                                </div>
                            </div>
                            """
                        else:
                            # No parameters available
                            params_html += '<div class="inline-muted">Parameters not available</div>'
                    
                    params_html += '</div>'
            else:
                # Single archive - show parameters normally
                params = parse_start_file(data_path)
                if params:
                    is_irreversible = params.get('flag') == 'i'
                    dynamics_label = "IRREVERSIBLE" if is_irreversible else "REVERSIBLE"
                    dynamics_class = "chip-dynamics-irr" if is_irreversible else "chip-dynamics-rev"
                    
                    # Build status chip HTML
                    status_class = status_chip_class(status)
                    status_chip_html = f'<div class="chip chip-status {status_class}">{status}</div>'
                    
                    params_html = f"""
                    <div class="params-panel">
                        <div class="params-row">
                            <div class="chip {dynamics_class}">{dynamics_label}</div>
                            {status_chip_html}
                            <div class="params-values">
                                <strong>Lattice:</strong> n={format_param(params.get('n', 'N/A'))} | 
                                <strong>Sweeps:</strong> s={format_param(params.get('sweeps', 'N/A'))} | 
                                <strong>Radius:</strong> r={format_param(params.get('radius', 'N/A'))} | 
                                <strong>Runs:</strong> m={format_param(params.get('runs', 'N/A'))}
                            </div>
                        </div>
                    </div>
                    """
            
            # Get notes for this run
            notes = read_notes(data_path)
            import html
            notes_escaped = html.escape(notes) if notes else ""
            
            notes_html = f"""
            <div class="notes-panel">
                <div class="notes-header">
                    <strong>Notes:</strong>
                    <div class="notes-actions">
                        <button id="editNotesBtn" class="notes-edit-btn" onclick="toggleEditMode()">Edit</button>
                        <span id="saveStatus" class="save-status"></span>
                        <button id="saveNotesBtn" class="notes-save-btn" onclick="saveNotes()">Save</button>
                    </div>
                </div>
                <div id="notesDisplay" class="notes-display" onclick="if(!isEditMode) toggleEditMode()">{notes_escaped if notes else '<span class="inline-muted">No notes yet. Click to add notes.</span>'}</div>
                <textarea id="notesTextarea" class="notes-textarea" data-dirname="{dirname}" placeholder="Add notes about this simulation run..." rows="10">{notes_escaped}</textarea>
            </div>
            """
            
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>Simulation Plots</title>
                <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={FAVICON_VERSION}">
                <script>
                    // Set theme immediately to prevent flash
                    (function() {{
                        const urlParams = new URLSearchParams(window.location.search);
                        const urlTheme = urlParams.get('theme');
                        const savedTheme = urlTheme || localStorage.getItem('theme') || 'dark';
                        document.documentElement.setAttribute('data-theme', savedTheme);
                    }})();
                </script>
                <style>
                    :root {{
                        --bg-primary: #1e1e1e;
                        --bg-secondary: #2d2d2d;
                        --text-primary: #e0e0e0;
                        --text-secondary: #b0b0b0;
                        --border-color: #404040;
                        --param-bg: #383838;
                        --msg-muted: #9aa0a6;
                        --msg-success: #28a745;
                    }}
                    
                    [data-theme="light"] {{
                        --bg-primary: #ffffff;
                        --bg-secondary: #f8f9fa;
                        --text-primary: #333333;
                        --text-secondary: #666666;
                        --border-color: #dddddd;
                        --param-bg: #f0f0f0;
                        --msg-muted: #6c757d;
                        --msg-success: #1e7e34;
                    }}
                    
                    body {{ 
                        margin: 0; 
                        padding: 20px; 
                        font-family: Arial, sans-serif;
                        background: var(--bg-primary);
                        color: var(--text-primary);
                        transition: background-color 0.3s ease, color 0.3s ease;
                    }}
                    h1 {{ 
                        text-align: center;
                        color: var(--text-primary);
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        gap: 15px;
                    }}
                    .plot-container {{ margin: 20px 0; }}
                    .back-link {{ 
                        margin-bottom: 20px;
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                    }}
                    .back-link a {{ 
                        color: #007bff;
                        font-size: 34px;
                        font-weight: 800;
                        line-height: 1;
                        text-decoration: none;
                        padding: 0px 14px 4px 14px;
                        border-radius: 4px;
                        align-items: center;
                        justify-content: center;
                        transition: all 0.2s;
                    }}
                    .back-link a:hover {{ 
                        background: #007bff;
                        color: white;
                        text-decoration: none;
                    }}
                    .refresh-btn {{
                        background: var(--bg-secondary);
                        color: var(--text-primary);
                        border: 1px solid var(--border-color);
                        padding: 8px;
                        border-radius: 4px;
                        cursor: pointer;
                        font-size: 18px;
                        transition: all 0.3s ease;
                        display: flex;
                        align-items: center;
                        height: 38px;
                    }}
                    .refresh-btn:hover {{
                        opacity: 0.8;
                        transform: scale(1.05);
                    }}
                    .theme-toggle {{
                        background: var(--bg-secondary);
                        color: var(--text-primary);
                        border: 1px solid var(--border-color);
                        padding: 8px;
                        border-radius: 4px;
                        cursor: pointer;
                        font-size: 18px;
                        transition: all 0.3s ease;
                        display: flex;
                        align-items: center;
                        height: 38px;
                    }}
                    .theme-toggle:hover {{
                        opacity: 0.8;
                        transform: scale(1.05);
                    }}
                    .export-btn {{
                        background: #28a745;
                        color: white;
                        border: none;
                        padding: 8px 16px;
                        border-radius: 4px;
                        cursor: pointer;
                        font-size: 14px;
                        font-weight: 500;
                        transition: all 0.3s ease;
                        display: flex;
                        align-items: center;
                        height: 38px;
                    }}
                    .export-btn:hover {{
                        background: #218838;
                        transform: scale(1.05);
                    }}
                    .params-panel {{
                        text-align: center;
                        margin: 20px 0;
                        padding: 15px;
                        background: var(--bg-secondary);
                        color: var(--text-primary);
                        border-radius: 5px;
                    }}
                    .params-row {{
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        gap: 15px;
                        flex-wrap: wrap;
                    }}
                    .params-values {{
                        padding: 10px;
                        background: var(--param-bg);
                        border-radius: 4px;
                    }}
                    .params-subpanel {{
                        margin-top: 10px;
                        padding: 10px;
                        background: var(--param-bg);
                        border-radius: 4px;
                    }}
                    .chip-combined {{
                        display: inline-block;
                        padding: 4px 12px;
                        border-radius: 12px;
                        font-size: 0.85em;
                        font-weight: bold;
                        background: linear-gradient(90deg, #ab63fa 0%, #00b8d4 100%);
                        color: white;
                    }}
                    .chip {{
                        display: inline-block;
                        padding: 4px 12px;
                        border-radius: 12px;
                        font-size: 0.85em;
                        font-weight: bold;
                    }}
                    .chip-status.status-completed {{ background: #28a745; color: white; }}
                    .chip-status.status-running {{ background: #ff9800; color: white; }}
                    .chip-status.status-interrupted {{ background: #ffc107; color: black; }}
                    .chip-status.status-error {{ background: #dc3545; color: white; }}
                    .chip-status.status-unknown {{ background: #e2e3e5; color: #383d41; }}
                    .chip-dynamics-rev {{ background: #ab63fa; color: #ffffff; }}
                    .chip-dynamics-irr {{ background: #00b8d4; color: #000000; }}
                    .chip-spaced {{
                        margin-left: 8px;
                    }}
                    .inline-muted {{
                        color: var(--text-secondary);
                    }}
                    .param-warning {{
                        background: #ff9800;
                        color: white;
                        padding: 8px 12px;
                        border-radius: 4px;
                        margin-top: 8px;
                        font-size: 0.9em;
                    }}
                    .save-status {{
                        display: none;
                        margin-left: 8px;
                        font-size: 14px;
                    }}
                    .save-status-muted {{
                        color: var(--msg-muted);
                    }}
                    .save-status-success {{
                        color: var(--msg-success);
                    }}
                    .notes-panel {{
                        margin: 20px 0;
                        padding: 15px;
                        background: var(--bg-secondary);
                        border-radius: 5px;
                    }}
                    .notes-header {{
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        margin-bottom: 10px;
                    }}
                    .notes-actions {{
                        display: flex;
                        align-items: center;
                    }}
                    .top-controls {{
                        display: flex;
                        gap: 8px;
                        align-items: center;
                    }}
                    .notes-edit-btn {{
                        padding: 4px 10px;
                        background: transparent;
                        color: #007bff;
                        border: 2px solid #007bff;
                        border-radius: 3px;
                        cursor: pointer;
                        font-size: 0.85em;
                        font-weight: bold;
                        transition: all 0.2s;
                    }}
                    .notes-edit-btn:hover {{
                        background: #007bff;
                        color: white;
                    }}
                    .notes-save-btn {{
                        padding: 6px 12px;
                        background: #28a745;
                        color: white;
                        border: none;
                        border-radius: 4px;
                        cursor: pointer;
                        font-size: 14px;
                        display: none;
                        margin-left: 8px;
                    }}
                    .notes-display {{
                        white-space: pre-wrap;
                        min-height: 50px;
                        padding: 10px;
                        background: var(--bg-primary);
                        color: var(--text-primary);
                        border-radius: 4px;
                        font-size: 11pt;
                        cursor: pointer;
                    }}
                    .notes-textarea {{
                        display: none;
                        width: 100%;
                        padding: 10px;
                        border: 1px solid var(--border-color);
                        border-radius: 4px;
                        font-family: inherit;
                        font-size: 14px;
                        resize: vertical;
                        box-sizing: border-box;
                        background: var(--bg-primary);
                        color: var(--text-primary);
                    }}
                </style>
                <script>
                    // Theme management
                    function initTheme() {{
                        // Check URL parameter first, then localStorage, default to dark
                        const urlParams = new URLSearchParams(window.location.search);
                        const urlTheme = urlParams.get('theme');
                        const savedTheme = urlTheme || localStorage.getItem('theme') || 'dark';
                        document.documentElement.setAttribute('data-theme', savedTheme);
                        updateThemeIcon();
                    }}
                    
                    function toggleTheme() {{
                        const currentTheme = document.documentElement.getAttribute('data-theme');
                        const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
                        localStorage.setItem('theme', newTheme);
                        
                        // Reload the page with the new theme parameter
                        const url = new URL(window.location);
                        url.searchParams.set('theme', newTheme);
                        window.location.href = url.toString();
                    }}
                    
                    function updateThemeIcon() {{
                        const theme = document.documentElement.getAttribute('data-theme');
                        const btn = document.getElementById('themeToggle');
                        if (btn) {{
                            btn.textContent = theme === 'dark' ? '☀️' : '🌙';
                            btn.title = theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';
                        }}
                    }}
                    
                    let originalNotes = '';
                    let isEditMode = false;
                    
                    document.addEventListener('DOMContentLoaded', function() {{
                        initTheme();  // Initialize theme on page load
                        
                        const textarea = document.getElementById('notesTextarea');
                        const saveBtn = document.getElementById('saveNotesBtn');
                        const saveStatus = document.getElementById('saveStatus');
                        originalNotes = textarea.value;
                        
                        // Show save button when content changes
                        textarea.addEventListener('input', function() {{
                            if (isEditMode && textarea.value !== originalNotes) {{
                                saveBtn.style.display = 'inline-block';
                                saveStatus.style.display = 'none';
                            }} else {{
                                saveBtn.style.display = 'none';
                            }}
                        }});
                        
                        // Save with Ctrl+S or Cmd+S
                        textarea.addEventListener('keydown', function(e) {{
                            if ((e.ctrlKey || e.metaKey) && e.key === 's') {{
                                e.preventDefault();
                                if (isEditMode) {{
                                    saveNotes();
                                }}
                            }}
                        }});
                        
                        // Auto-save on blur (when clicking away from textarea)
                        textarea.addEventListener('blur', function() {{
                            if (isEditMode && textarea.value !== originalNotes) {{
                                saveNotes();
                            }}
                        }});
                    }});
                    
                    function toggleEditMode() {{
                        isEditMode = true;
                        const saveStatus = document.getElementById('saveStatus');
                        document.getElementById('notesDisplay').style.display = 'none';
                        document.getElementById('notesTextarea').style.display = 'block';
                        document.getElementById('editNotesBtn').style.display = 'none';
                        saveStatus.style.display = 'none';
                        saveStatus.className = 'save-status';
                        document.getElementById('notesTextarea').focus();
                    }}
                    
                    function saveNotes() {{
                        const textarea = document.getElementById('notesTextarea');
                        const saveBtn = document.getElementById('saveNotesBtn');
                        const saveStatus = document.getElementById('saveStatus');
                        const dirname = textarea.getAttribute('data-dirname');
                        const notes = textarea.value;
                        
                        saveBtn.disabled = true;
                        saveBtn.style.display = 'none';
                        saveStatus.style.display = 'inline-block';
                        saveStatus.textContent = 'Saving...';
                        saveStatus.className = 'save-status save-status-muted';
                        
                        fetch(window.location.origin + '/notes/' + dirname, {{
                            method: 'POST',
                            headers: {{
                                'Content-Type': 'application/json'
                            }},
                            body: JSON.stringify({{ notes: notes }})
                        }})
                        .then(response => response.json())
                        .then(data => {{
                            if (data.success) {{
                                originalNotes = notes;
                                saveBtn.disabled = false;
                                saveStatus.textContent = '✓ Saved';
                                saveStatus.className = 'save-status save-status-success';
                                
                                setTimeout(() => {{
                                    saveStatus.style.display = 'none';
                                    saveStatus.className = 'save-status';
                                }}, 2000);
                            }} else {{
                                saveBtn.disabled = false;
                                saveBtn.style.display = 'inline-block';
                                saveStatus.style.display = 'none';
                                alert('Error saving notes: ' + data.error);
                            }}
                        }})
                        .catch(error => {{
                            saveBtn.disabled = false;
                            saveBtn.style.display = 'inline-block';
                            saveStatus.style.display = 'none';
                            alert('Error saving notes');
                        }});
                    }}
                </script>
            </head>
            <body>
                <div class="back-link">
                    <a href="/" title="Back to simulation browser">‹</a>
                    <div class="top-controls">
                        <button class="export-btn" onclick="window.location.href='/export-loading/{dirname}?theme=' + (document.documentElement.getAttribute('data-theme') || localStorage.getItem('theme') || 'dark')" title="Export this simulation archive">Export</button>
                        <button class="refresh-btn" onclick="location.reload()" title="Refresh the page to see latest data">🔄</button>
                        <button id="themeToggle" class="theme-toggle" onclick="toggleTheme()" title="Switch between dark and light themes">🌙</button>
                    </div>
                </div>
                <h1>Simulation Plots - {title}</h1>
                {params_html}
                {notes_html}
            """
            
            if plot1_html:
                html_content += f'<div class="plot-container">{plot1_html}</div>'
            
            if plot2_html:
                html_content += f'<div class="plot-container">{plot2_html}</div>'
            
            html_content += """
            </body>
            </html>
            """
            
            # Add cache control headers to prevent browser caching
            response = make_response(html_content)
            response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
            response.headers['Pragma'] = 'no-cache'
            response.headers['Expires'] = '0'
            return response

    except Exception as e:
        import traceback
        return f"<h1>Error generating plots</h1><pre>{traceback.format_exc()}</pre>", 500

if __name__ == '__main__':
    print("\n" + "="*60)
    print("  Simulation Browser")
    print("="*60)
    print(f"\n  📂 Data directory: {DATA_DIR}")
    print(f"  📦 Archive directory: {ARCHIVE_DIR}")
    print(f"\n  🌐 Open in browser: http://127.0.0.1:5001")
    print("\n  Press Ctrl+C to stop\n")
    
    app.run(debug=True, port=5001, host='127.0.0.1')
