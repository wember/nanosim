#!/usr/bin/env venv/bin/python
"""Simple web interface to browse current and archived simulation runs."""

from flask import Flask, render_template_string, send_file, request, jsonify, session, Response, stream_with_context
from pathlib import Path
from datetime import datetime
from urllib.parse import quote
import hashlib
import json
import logging
import os
import re
import zipfile
import tempfile
import shutil
import secrets

app = Flask(__name__)
app.secret_key = secrets.token_hex(16)  # For session management
app.logger.setLevel(logging.INFO)

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
DELETED_DIR = REPO_ROOT / 'archive'  # Non-destructive archive destination
EXPORT_EXTENSION = '.nanosim'
EXPORT_MODE_FULL = 'full'
EXPORT_MODE_PLOTS_ONLY = 'plots-only'
PLOTS_ONLY_ROOT_FILES = {
    'sim_started.txt',
    'sim_status.txt',
    'sim_completed.txt',
    'sim_notes.txt',
}
PLOT_CACHE_VERSION = 1
PLOT_CACHE_THEMES = ('dark', 'light')
PLOT_CACHE_STATIC_PREFIX = os.getenv('NANOSIM_PLOT_CACHE_STATIC_PREFIX', '/plot-cache').rstrip('/')
WEB_PLOT_BUILD_POLICY = os.getenv('NANOSIM_WEB_PLOT_BUILD', 'localhost-only').strip().lower()

# Cache for completed exports, keyed by token.
# Populated by /export-stream; consumed by /export-download.
_export_cache: dict = {}

FAVICON_SVG = "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'><rect width='32' height='32' rx='6' fill='%23ab63fa'/><polyline points='8,26 8,6 24,26 24,6' fill='none' stroke='white' stroke-width='4' stroke-linecap='round' stroke-linejoin='round'/></svg>"
FAVICON_VERSION = hashlib.md5(FAVICON_SVG.encode()).hexdigest()[:8]


def is_loopback_host(host_value: str) -> bool:
    """Return True if a Host header resolves to localhost/loopback."""
    host = (host_value or '').strip().lower()
    if not host:
        return False

    if host.startswith('localhost') or host.startswith('127.'):
        return True
    if host == '::1' or host.startswith('[::1]'):
        return True

    return False


def can_run_web_plot_build(req) -> bool:
    """Policy gate for expensive plot builds triggered by web requests."""
    if WEB_PLOT_BUILD_POLICY in ('enabled', 'allow-all', 'true', '1'):
        return True
    if WEB_PLOT_BUILD_POLICY in ('disabled', 'cache-only', 'false', '0'):
        return False

    # Default: localhost-only
    host = req.host or ''
    forwarded_host = req.headers.get('X-Forwarded-Host', '')
    return is_loopback_host(host) or is_loopback_host(forwarded_host)


def can_show_data_transfer_controls(req) -> bool:
    """Show import/export controls only for localhost sessions."""
    host = req.host or ''
    forwarded_host = req.headers.get('X-Forwarded-Host', '')
    return is_loopback_host(host) or is_loopback_host(forwarded_host)


app.logger.info(
    'plot_cache web_build_policy=%s',
    WEB_PLOT_BUILD_POLICY,
)


def plot_payload_total_bytes(plot1_html: str | None, plot2_html: str | None) -> int:
    """Return total UTF-8 payload size for two optional plot HTML strings."""
    total = 0
    if plot1_html is not None:
        total += len(plot1_html.encode('utf-8'))
    if plot2_html is not None:
        total += len(plot2_html.encode('utf-8'))
    return total


def to_mb(size_bytes: int) -> float:
    """Convert bytes to MB for compact logging."""
    return size_bytes / (1024 * 1024)


def normalize_export_mode(mode_value: str | None) -> str:
    """Normalize export mode to one of the supported values."""
    if (mode_value or '').strip().lower() == EXPORT_MODE_PLOTS_ONLY:
        return EXPORT_MODE_PLOTS_ONLY
    return EXPORT_MODE_FULL


def normalize_theme(raw_theme: str | None) -> str:
    """Normalize theme values from URLs and avoid duplicate '?theme=' fragments."""
    theme = (raw_theme or 'dark').strip()
    if not theme:
        return 'dark'
    if '?theme=' in theme:
        theme = theme.split('?theme=', 1)[0]
    if theme.endswith('?'):
        theme = theme[:-1]
    return theme if theme in {'dark', 'light'} else 'dark'


def export_file_list(archive_path: Path, export_mode: str) -> list[Path]:
    """Return the files to include for the selected export mode."""
    all_files = sorted([p for p in archive_path.rglob('*') if p.is_file()])
    if export_mode == EXPORT_MODE_FULL:
        return all_files

    selected_files = []
    for file_path in all_files:
        rel_path = file_path.relative_to(archive_path)
        rel_parts = rel_path.parts
        if not rel_parts:
            continue

        if rel_parts[0] == 'plot_cache':
            selected_files.append(file_path)
            continue

        if len(rel_parts) == 1 and rel_parts[0] in PLOTS_ONLY_ROOT_FILES:
            selected_files.append(file_path)

    return selected_files


def export_filename(dirname: str, export_mode: str) -> str:
    """Return download filename for an export mode."""
    if export_mode == EXPORT_MODE_PLOTS_ONLY:
        return f'ember_nanosim_plots_{dirname}{EXPORT_EXTENSION}'
    return f'ember_nanosim_{dirname}{EXPORT_EXTENSION}'


def export_placeholder_dirs(archive_path: Path, export_mode: str) -> list[str]:
    """Return empty directory placeholders to preserve for plots-only exports."""
    if export_mode != EXPORT_MODE_PLOTS_ONLY:
        return []

    placeholders = []
    if (archive_path / 'rev').is_dir():
        placeholders.append('rev')
    if (archive_path / 'irr').is_dir():
        placeholders.append('irr')
    return placeholders


def write_zip_dir_entry(zipf: zipfile.ZipFile, dirname: str) -> None:
    """Write an explicit empty directory entry into a zip file."""
    zipf.writestr(f"{dirname.rstrip('/')}/", b'')


def resolve_plot_data_path(dirname: str) -> Path:
    if dirname == 'current':
        return DATA_DIR
    return ARCHIVE_DIR / dirname


def plot_cache_paths(data_path: Path, theme: str, plot_options: dict[str, bool] | None = None) -> dict[str, Path]:
    safe_theme = normalize_theme(theme)
    cache_dir = data_path / 'plot_cache'
    return {
        'dir': cache_dir,
        'plot1': cache_dir / f'plot1_{safe_theme}.html',
        'plot2': cache_dir / f'plot2_{safe_theme}.html',
        'meta': cache_dir / f'manifest_{safe_theme}.json',
    }


def plot_options_match_manifest(manifest: dict | None, plot_options: dict[str, bool] | None = None) -> bool:
    if not isinstance(manifest, dict):
        return False
    stored = manifest.get('plot_options')
    if not isinstance(stored, dict):
        return False
    expected = parse_plot_options(fallback=plot_options)
    for key in default_plot_options():
        if stored.get(key, True) != expected.get(key, True):
            return False
    return True


def plot_cache_public_url(dirname: str, filename: str, cache_bust: int | None = None) -> str:
    """Return the public URL for a cached plot artifact.

    In production, nginx should serve this path directly from disk.
    In local development, Flask provides a fallback route on the same path.
    """
    safe_dir = quote(dirname, safe='')
    safe_file = quote(filename, safe='')
    if cache_bust is None:
        cache_bust = int((resolve_plot_data_path(dirname) / 'plot_cache' / filename).stat().st_mtime_ns) if (resolve_plot_data_path(dirname) / 'plot_cache' / filename).exists() else 0
    return f"{PLOT_CACHE_STATIC_PREFIX}/{safe_dir}/plot_cache/{safe_file}?v={cache_bust}"


def default_plot_options() -> dict[str, bool]:
    return {
        'include_entropy': True,
        'include_zoom': True,
        'include_psd': True,
        'include_summary': True,
    }


def parse_plot_options(request_args=None, fallback: dict[str, bool] | None = None) -> dict[str, bool]:
    options = dict(fallback or default_plot_options())
    if request_args is None:
        return options

    for key, default_value in default_plot_options().items():
        raw_value = request_args.get(key)
        if raw_value is None:
            options[key] = default_value
            continue
        if raw_value in ('0', 'false', 'False', 'no', 'No'):
            options[key] = False
        elif raw_value in ('1', 'true', 'True', 'yes', 'Yes'):
            options[key] = True
        else:
            options[key] = bool(raw_value) if raw_value not in ('', 'None') else default_value
    return options


def plot_options_query_string(options: dict[str, bool] | None = None) -> str:
    opts = parse_plot_options(fallback=options)
    return '&'.join(f"{key}={'1' if value else '0'}" for key, value in (
        ('include_entropy', opts.get('include_entropy', True)),
        ('include_zoom', opts.get('include_zoom', True)),
        ('include_psd', opts.get('include_psd', True)),
        ('include_summary', opts.get('include_summary', True)),
    ))


def plot_options_checkbox_html(options: dict[str, bool] | None = None) -> str:
    opts = parse_plot_options(fallback=options)
    label_map = {
        'include_entropy': 'Entropy',
        'include_zoom': 'Zoom',
        'include_psd': 'PSD',
        'include_summary': 'Summary',
    }
    checkbox_html = []
    for key, label in label_map.items():
        checked = 'checked' if opts.get(key, True) else ''
        checkbox_html.append(
            f'<label class="plot-option-picker"><input type="checkbox" data-key="{key}" {checked}> {label}</label>'
        )
    return '<div class="plot-option-set">' + ' '.join(checkbox_html) + '</div>'


def plot_options_status_html(options: dict[str, bool] | None = None) -> str:
    opts = parse_plot_options(fallback=options)
    label_map = {
        'include_entropy': 'Entropy',
        'include_zoom': 'Zoom',
        'include_psd': 'PSD',
        'include_summary': 'Summary',
    }
    chips = []
    for key, label in label_map.items():
        included = opts.get(key, True)
        chip_class = 'included' if included else 'hidden'
        chips.append(
            f'<span class="plot-option-chip {chip_class}" data-key="{key}">{label}</span>'
        )
    return '<div class="plot-option-set plot-option-readout">' + ' '.join(chips) + '</div>'


def write_plot_cache_artifacts(
    data_path: Path,
    theme: str,
    plot1_html: str | None,
    plot2_html: str | None,
    source: str,
    plot_options: dict[str, bool] | None = None,
) -> dict:
    options = parse_plot_options(fallback=plot_options)
    paths = plot_cache_paths(data_path, theme, options)
    tmp_root = Path(tempfile.mkdtemp(prefix='plot_cache_tmp_', dir=str(data_path)))
    tmp_cache_dir = tmp_root / 'plot_cache'
    tmp_cache_dir.mkdir(parents=True, exist_ok=True)

    if plot1_html is not None:
        (tmp_cache_dir / paths['plot1'].name).write_text(plot1_html)
    if plot2_html is not None:
        (tmp_cache_dir / paths['plot2'].name).write_text(plot2_html)

    manifest = {
        'cache_version': PLOT_CACHE_VERSION,
        'theme': theme,
        'source': source,
        'built_at': datetime.now().isoformat(),
        'plot_options': options,
    }
    (tmp_cache_dir / paths['meta'].name).write_text(json.dumps(manifest, indent=2))

    paths['dir'].mkdir(parents=True, exist_ok=True)

    for stale_path in paths['dir'].glob(f'plot1_{normalize_theme(theme)}*.html'):
        if stale_path.name != paths['plot1'].name:
            stale_path.unlink(missing_ok=True)
    for stale_path in paths['dir'].glob(f'plot2_{normalize_theme(theme)}*.html'):
        if stale_path.name != paths['plot2'].name:
            stale_path.unlink(missing_ok=True)
    for stale_path in paths['dir'].glob(f'manifest_{normalize_theme(theme)}*.json'):
        if stale_path.name != paths['meta'].name:
            stale_path.unlink(missing_ok=True)

    for entry in tmp_cache_dir.iterdir():
        shutil.move(str(entry), str(paths['dir'] / entry.name))

    shutil.rmtree(tmp_root, ignore_errors=True)
    return manifest


def build_plot_cache(dirname: str, theme: str, source: str = 'local-precompute', log_fn=None,
                    include_entropy: bool = True, include_zoom: bool = True,
                    include_psd: bool = True, include_summary: bool = True) -> tuple[dict, float]:
    import subprocess
    import time

    if theme not in PLOT_CACHE_THEMES:
        raise ValueError(f"Unsupported theme '{theme}'")

    data_path = resolve_plot_data_path(dirname)
    plot_script = REPO_ROOT / 'creutz-sim' / 'Sk_comparison.py'
    plotly_template = 'plotly_dark' if theme == 'dark' else 'plotly_white'
    emit = log_fn or (lambda _message: None)
    started_at = time.time()

    if not data_path.exists():
        raise FileNotFoundError(f"data path not found: {data_path}")
    if not plot_script.exists():
        raise FileNotFoundError(f"plotting script not found: {plot_script}")

    emit(f"cache_build_scan: preparing plot build inputs for {theme}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        script_content = plot_script.read_text()
        cli_args = [
            '--data-dir', str(data_path),
            '--include-entropy' if include_entropy else '--no-entropy',
            '--include-zoom' if include_zoom else '--no-zoom',
            '--include-psd' if include_psd else '--no-psd',
            '--include-summary' if include_summary else '--no-summary',
        ]
        modified_script = script_content.replace(
            "args = parser.parse_args()",
            f"args = parser.parse_args({cli_args!r})"
        ).replace(
            'pio.templates.default = "plotly_white"',
            f'pio.templates.default = "{plotly_template}"'
        ).replace(
            'show_plotly_figure(fig, "fig1")',
            f'fig.write_html(r\'{tmpdir_path}/plot1.html\')\n'
            'print("[web-build] wrote plot1.html", flush=True)'
        ).replace(
            'show_plotly_figure(fig2, "fig2")',
            f'fig2.write_html(r\'{tmpdir_path}/plot2.html\')\n'
            'print("[web-build] wrote plot2.html", flush=True)'
        )

        temp_script = tmpdir_path / 'plot_script.py'
        temp_script.write_text(modified_script)

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
                emit(line)

        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError("plot build subprocess exited with non-zero status")

        plot1 = tmpdir_path / 'plot1.html'
        plot2 = tmpdir_path / 'plot2.html'
        if not plot1.exists() and not plot2.exists():
            raise RuntimeError("no plot HTML files generated")

        emit(f"cache_build_write: writing cache artifacts for {theme}")
        plot1_html = plot1.read_text() if plot1.exists() else None
        plot2_html = plot2.read_text() if plot2.exists() else None
        cache_meta = write_plot_cache_artifacts(
            data_path,
            theme,
            plot1_html,
            plot2_html,
            source,
            plot_options={
                'include_entropy': include_entropy,
                'include_zoom': include_zoom,
                'include_psd': include_psd,
                'include_summary': include_summary,
            },
        )

    elapsed = time.time() - started_at
    emit(f"cache_build_done: {theme} cache ready in {elapsed:.2f}s")
    return cache_meta, elapsed


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
        .replot-link {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            padding: 2px 8px;
            border: 2px solid #f0b400;
            border-radius: 3px;
            font-size: 1em;
            line-height: 1.1;
            color: #b88400;
            background: rgba(240, 180, 0, 0.08);
            transition: all 0.2s;
            min-width: 80px;
            height: 28px;
            text-align: center;
            vertical-align: middle;
            appearance: none;
            -webkit-appearance: none;
            cursor: pointer;
            font: inherit;
            box-sizing: border-box;
        }
        .replot-link:hover {
            background: #f0b400;
            color: white;
            text-decoration: none;
        }
        .archive-link {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #fd7e14;
            border-radius: 3px;
            font-size: 1em;
            color: #fd7e14;
            transition: all 0.2s;
            min-width: 80px;
            text-align: center;
        }
        .archive-link:hover {
            background: #fd7e14;
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

        .export-menu {
            position: relative;
            display: inline-block;
        }

        .export-menu .export-link {
            cursor: pointer;
            appearance: none;
            -webkit-appearance: none;
            background: transparent;
            font: inherit;
            line-height: inherit;
        }

        .export-menu .export-link:focus-visible {
            outline: 2px solid #28a745;
            outline-offset: 2px;
        }

        .export-menu-items {
            position: absolute;
            right: 0;
            top: calc(100% + 6px);
            min-width: 170px;
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            box-shadow: 0 4px 10px var(--shadow);
            z-index: 30;
            padding: 6px;
            display: none;
        }

        .export-menu.open .export-menu-items {
            display: block;
        }

        .export-menu-item {
            display: block;
            padding: 7px 9px;
            border-radius: 4px;
            color: var(--text-primary);
            text-decoration: none;
            font-size: 0.9em;
        }

        .export-menu-item:hover {
            background: var(--bg-hover);
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
        #replotModal {
            display: none;
            position: fixed;
            z-index: 1200;
            width: min(260px, calc(100vw - 24px));
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            box-shadow: 0 6px 18px rgba(0, 0, 0, 0.16);
            padding: 8px 8px 6px;
            max-height: none;
            height: auto;
        }
        #replotModal::before {
            content: "";
            position: absolute;
            left: 18px;
            top: -8px;
            width: 12px;
            height: 12px;
            background: var(--bg-secondary);
            border-left: 1px solid var(--border-color);
            border-top: 1px solid var(--border-color);
            transform: rotate(45deg);
        }
        #replotModal .modal-content {
            margin: 0;
            padding: 0;
            width: 100%;
            max-width: none;
            max-height: none;
            overflow: visible;
            border: none;
            background: transparent;
        }
        #replotModal .modal-header {
            margin-bottom: 8px;
            padding-bottom: 6px;
            border-bottom: 1px solid var(--border-color);
        }
        #replotModal .modal-header h2 {
            font-size: 0.95rem;
            font-weight: 600;
        }
        #replotModal .close-btn {
            font-size: 20px;
            line-height: 1;
        }
        #replotOptionsForm {
            display: flex;
            flex-direction: column;
            gap: 4px;
            margin: 0 0 8px;
            padding: 0;
        }
        #replotOptionsForm .plot-option-picker {
            justify-content: flex-start;
            padding: 5px 7px;
            border-radius: 4px;
            font-size: 12px;
            background: var(--bg-primary);
            border-color: var(--border-color);
        }
        #replotModal .modal-footer {
            gap: 6px;
            margin-top: 8px;
            justify-content: flex-end;
        }
        #replotModal .modal-footer .btn {
            padding: 5px 10px;
            font-size: 12px;
        }
        .replot-link::after {
            content: '▾';
            margin-left: 6px;
            font-size: 0.8em;
            opacity: 0.8;
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
        .plot-actions-scope {
            display: none;
        }
        html[data-theme="dark"] .plot-actions-scope[data-theme-scope="dark"],
        html[data-theme="light"] .plot-actions-scope[data-theme-scope="light"] {
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
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
            .archive-link,
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
            .replot-link,
            .archive-link,
            .export-link {
                font-size: 0.8em;
                padding: 4px 10px;
            }
        }
    </style>
    <script>
        let currentDirname = null;
        let currentExportDirname = null;
        let currentExportMode = 'full';
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
            const theme = document.documentElement.getAttribute('data-theme') || localStorage.getItem('theme') || 'dark';
            const url = new URL(link.href, window.location.origin);
            url.searchParams.set('theme', theme);
            window.location.href = url.toString();
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

        function openReplotModal(dirname, buttonEl) {
            const modal = document.getElementById('replotModal');
            const form = document.getElementById('replotOptionsForm');
            if (modal.style.display === 'block' && modal.dataset.dirname === dirname) {
                closeReplotModal();
                return;
            }
            modal.dataset.dirname = dirname;
            form.querySelectorAll('input[type="checkbox"]').forEach(function(cb) {
                cb.checked = true;
            });
            const rect = buttonEl.getBoundingClientRect();
            const modalWidth = Math.min(260, window.innerWidth - 24);
            const left = Math.min(Math.max(8, rect.left), window.innerWidth - modalWidth - 8);
            const fitsBelow = (rect.bottom + 220) <= window.innerHeight - 12;
            const top = fitsBelow ? rect.bottom + 8 : rect.top - 8 - 220;
            modal.style.left = left + 'px';
            modal.style.top = Math.max(8, top) + 'px';
            modal.style.display = 'block';
            modal.style.height = 'auto';
        }

        function closeReplotModal() {
            document.getElementById('replotModal').style.display = 'none';
        }

        function startReplot() {
            const modal = document.getElementById('replotModal');
            const dirname = modal.dataset.dirname;
            const form = document.getElementById('replotOptionsForm');
            const query = [];
            form.querySelectorAll('input[type="checkbox"]').forEach(function(cb) {
                query.push(cb.dataset.key + '=' + (cb.checked ? '1' : '0'));
            });
            const theme = (document.documentElement.getAttribute('data-theme') || localStorage.getItem('theme') || 'dark');
            const url = '/plot-loading/' + dirname + '?theme=' + theme + '&' + query.join('&') + '&force_replot=1';
            closeReplotModal();
            window.location.href = url;
        }

        function exportModeLabel(exportMode) {
            return exportMode === 'plots-only' ? 'Plots Only' : 'Full Data Set';
        }

        function openExportModal(dirname, exportMode = 'full') {
            currentExportDirname = dirname;
            currentExportMode = exportMode;
            closeAllExportMenus();
            const modal = document.getElementById('exportModal');
            document.getElementById('exportStatus').textContent = 'Preparing export (' + exportModeLabel(exportMode) + ')...';
            document.getElementById('exportProgress').textContent = '';
            document.getElementById('exportLog').innerHTML = '';
            document.getElementById('exportHelp').style.display = 'block';
            document.getElementById('exportCloseBtn').style.display = 'none';
            document.getElementById('exportSpinner').style.display = 'block';
            modal.style.display = 'block';
            startExportStream(dirname, exportMode);
        }

        function closeExportModal() {
            if (exportEventSource) {
                exportEventSource.close();
                exportEventSource = null;
            }
            document.getElementById('exportModal').style.display = 'none';
        }

        function closeAllExportMenus(exceptMenu = null) {
            document.querySelectorAll('.export-menu').forEach((menu) => {
                if (menu === exceptMenu) return;
                menu.classList.remove('open');
                const trigger = menu.querySelector('.export-link');
                if (trigger) trigger.setAttribute('aria-expanded', 'false');
            });
        }

        function toggleExportMenu(event, triggerButton) {
            event.preventDefault();
            event.stopPropagation();
            const menu = triggerButton.closest('.export-menu');
            if (!menu) return;
            const isOpen = menu.classList.contains('open');
            closeAllExportMenus(menu);
            menu.classList.toggle('open', !isOpen);
            triggerButton.setAttribute('aria-expanded', String(!isOpen));
        }

        document.addEventListener('click', function(event) {
            if (!event.target.closest('.export-menu')) {
                closeAllExportMenus();
            }
        });

        document.addEventListener('keydown', function(event) {
            if (event.key === 'Escape') {
                closeAllExportMenus();
            }
        });

        function appendExportLog(text, isError) {
            const logBox = document.getElementById('exportLog');
            const p = document.createElement('p');
            p.className = 'log-line' + (isError ? ' log-error' : '');
            p.textContent = text;
            logBox.appendChild(p);
            logBox.scrollTop = logBox.scrollHeight;
            while (logBox.children.length > 300) logBox.removeChild(logBox.firstChild);
        }

        function startExportStream(dirname, exportMode = 'full') {
            if (exportEventSource) {
                exportEventSource.close();
            }

            const statusEl = document.getElementById('exportStatus');
            const progressEl = document.getElementById('exportProgress');
            const helpTextEl = document.getElementById('exportHelp');
            const closeBtnEl = document.getElementById('exportCloseBtn');
            const spinnerEl = document.getElementById('exportSpinner');

            const theme = document.documentElement.getAttribute('data-theme') || localStorage.getItem('theme') || 'dark';
            const streamUrl = '/export-stream/' + dirname + '?theme=' + theme + '&export_mode=' + encodeURIComponent(exportMode);
            exportEventSource = new EventSource(streamUrl);

            exportEventSource.onopen = function() {
                statusEl.textContent = 'Preparing export (' + exportModeLabel(exportMode) + ')...';
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
            const replotModal = document.getElementById('replotModal');
            const importModal = document.getElementById('importModal');
            const combineModal = document.getElementById('combineModal');
            const exportModal = document.getElementById('exportModal');
            if (event.target == notesModal) closeNotesModal();
            if (event.target == replotModal) closeReplotModal();
            if (replotModal && replotModal.style.display !== 'none' && !event.target.closest('#replotModal') && !event.target.closest('.replot-link')) {
                closeReplotModal();
            }
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
        
        function archiveRun(dirname) {
            if (confirm('Move this run to the archive folder? It will no longer appear in the list but can be recovered from the archive/ directory.\\n\\nRun: ' + dirname)) {
                fetch('/archive-run/' + dirname, {
                    method: 'POST'
                })
                .then(response => response.json())
                .then(data => {
                    if (data.success) {
                        sessionStorage.setItem('toastMessage', 'Run archived successfully');
                        window.location.reload();
                    } else {
                        alert('Error archiving run: ' + data.error);
                    }
                })
                .catch(error => {
                    alert('Error archiving run: ' + error);
                });
            }
        }
        
        function exportArchive(dirname, exportMode = 'full') {
            openExportModal(dirname, exportMode);
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

    <div id="replotModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h2>Replot Options</h2>
                <button class="close-btn" onclick="closeReplotModal()">&times;</button>
            </div>
            <div id="replotOptionsForm" class="plot-option-set" style="display:flex; justify-content:flex-start; margin:0 0 16px;">
                <label class="plot-option-picker"><input type="checkbox" data-key="include_entropy" checked> Entropy</label>
                <label class="plot-option-picker"><input type="checkbox" data-key="include_zoom" checked> Zoom</label>
                <label class="plot-option-picker"><input type="checkbox" data-key="include_psd" checked> PSD</label>
                <label class="plot-option-picker"><input type="checkbox" data-key="include_summary" checked> Summary</label>
            </div>
            <div class="modal-footer">
                <button class="btn btn-secondary" onclick="closeReplotModal()">Cancel</button>
                <button class="btn btn-primary" onclick="startReplot()">Start</button>
            </div>
        </div>
    </div>
    
    {% if show_data_transfer_controls %}
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
    {% endif %}
    
    {% if show_data_transfer_controls %}
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
    {% endif %}

    {% if show_data_transfer_controls %}
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
    {% endif %}
    
    <h1>
        <span>📊 Simulation Browser</span>
        <div class="header-controls">
            {% if show_data_transfer_controls %}
            <button class="import-btn" onclick="openImportModal()" title="Import a previously exported simulation archive">Import</button>
            <button class="combine-btn" onclick="openCombineModal()" title="Combine reversible and irreversible simulations into one archive">Combine</button>
            {% endif %}
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
                    {% for scope_theme in ['dark', 'light'] %}
                    <div class="plot-actions-scope" data-theme-scope="{{ scope_theme }}">
                        {% if archive.has_plot_cache_by_theme[scope_theme] %}
                        <a href="/plot/{{ archive.dirname }}" class="plot-link" data-base-url="/plot/{{ archive.dirname }}" onclick="this.href = this.dataset.baseUrl; addThemeToUrl(event, this); return true;" title="Open the cached plot if available, or build it on demand">View plots</a>
                        {% if show_web_plot_build_controls %}
                        <button type="button" class="replot-link" onclick="openReplotModal('{{ archive.dirname }}', this); return false;" title="Rebuild the plot cache from the source data using the selected panels">Replot</button>
                        {% endif %}
                        {% elif show_web_plot_build_controls %}
                        <button type="button" class="replot-link" onclick="openReplotModal('{{ archive.dirname }}', this); return false;" title="Build the plot cache from the source data using the selected panels">Plot</button>
                        {% endif %}
                    </div>
                    {% endfor %}
                    <a href="/view/{{ archive.dirname }}" class="view-link" title="Browse all files in this archive">View files</a>
                    {% if not archive.notes %}
                    <a href="#" class="notes-link" data-dirname="{{ archive.dirname }}" data-notes="" onclick="openNotesModalFromLink(this); return false;" title="Add notes about this simulation">Add notes</a>
                    {% endif %}
                    {% if show_data_transfer_controls %}
                    <div class="export-menu">
                        <button type="button" class="export-link" title="Export this simulation archive" aria-haspopup="true" aria-expanded="false" onclick="toggleExportMenu(event, this)">Export</button>
                        <div class="export-menu-items">
                            <a href="#" class="export-menu-item" onclick="exportArchive('{{ archive.dirname }}', 'full'); return false;">Full Data Set</a>
                            <a href="#" class="export-menu-item" onclick="exportArchive('{{ archive.dirname }}', 'plots-only'); return false;">Plots Only</a>
                        </div>
                    </div>
                    {% endif %}
                    {% if show_data_transfer_controls %}
                    <a href="#" class="archive-link" onclick="archiveRun('{{ archive.dirname }}'); return false;" title="Move this run to the archive folder">Archive</a>
                    {% endif %}
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

RESTORE_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Restore Archived Run</title>
    <link rel="icon" type="image/svg+xml" href="/favicon.svg?v={{ favicon_version }}">
    <script>
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
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            max-width: 900px;
            margin: 0 auto;
            padding: 30px 20px;
        }
        h1 { color: var(--text-primary); margin-bottom: 20px; }
        .page-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 24px;
        }
        .back-link a {
            color: #007bff;
            font-size: 34px;
            font-weight: 800;
            line-height: 1;
            text-decoration: none;
            padding: 0px 14px 4px 14px;
            border-radius: 4px;
            transition: all 0.2s;
        }
        .back-link a:hover { background: #007bff; color: white; text-decoration: none; }
        .theme-toggle {
            background: var(--bg-secondary);
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 18px;
            display: flex;
            align-items: center;
            height: 38px;
        }
        .archive-list { background: transparent; box-shadow: none; }
        .archive-item {
            padding: 15px;
            margin-bottom: 12px;
            border-radius: 8px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 4px var(--shadow);
        }
        .archive-item:hover { background: var(--bg-hover); }
        .timestamp {
            font-size: 1.2em;
            font-weight: bold;
            color: var(--text-primary);
            margin-bottom: 8px;
        }
        .archive-chip-row {
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
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
        .dynamics-chip.reversible { background: #ab63fa; color: white; }
        .dynamics-chip.irreversible { background: #00b8d4; color: black; }
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
        .meta-block { margin-top: 8px; }
        .notes-preview {
            white-space: pre-wrap;
            margin-top: 4px;
            max-height: 60px;
            overflow: hidden;
        }
        .action-buttons {
            margin-top: 10px;
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }
        .restore-btn {
            display: inline-block;
            padding: 2px 8px;
            border: 2px solid #28a745;
            border-radius: 3px;
            font-size: 1em;
            color: #28a745;
            background: transparent;
            cursor: pointer;
            min-width: 80px;
            text-align: center;
            transition: all 0.2s;
        }
        .restore-btn:hover { background: #28a745; color: white; }
        .restore-btn:disabled { border-color: #6c757d; color: #6c757d; cursor: not-allowed; }
        .no-archives {
            text-align: center;
            padding: 60px 20px;
            color: var(--text-secondary);
        }
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
            display: none;
        }
    </style>
</head>
<body>
    <div class="page-header">
        <div class="back-link"><a href="/" title="Back to simulation browser">‹</a></div>
        <h1>Restore Archived Run</h1>
    </div>

    {% if archived %}
    <div class="archive-list">
        {% for archive in archived %}
        <div class="archive-item" id="row-{{ archive.dirname }}">
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
                {% if archive.rev_params %}
                <div class="details-grid">
                    <div><strong>Lattice:</strong> <span class="param">n={{ archive.rev_params.n }}</span></div>
                    <div><strong>Sweeps:</strong> <span class="param">s={{ archive.rev_params.sweeps }}</span></div>
                    <div><strong>Radius:</strong> <span class="param">r={{ archive.rev_params.radius }}</span></div>
                    <div><strong>Runs:</strong> <span class="param">m={{ archive.rev_params.runs }}</span></div>
                </div>
                {% endif %}
                {% elif archive.params %}
                <div class="details-grid">
                    <div><strong>Lattice:</strong> <span class="param">n={{ archive.params.n }}</span></div>
                    <div><strong>Sweeps:</strong> <span class="param">s={{ archive.params.sweeps }}</span></div>
                    <div><strong>Radius:</strong> <span class="param">r={{ archive.params.radius }}</span></div>
                    <div><strong>Runs:</strong> <span class="param">m={{ archive.params.runs }}</span></div>
                    <div><strong>Total sims:</strong> <span class="param">{{ archive.params.total }}</span></div>
                </div>
                {% endif %}
                {% if archive.completion_info %}
                <div class="meta-block">
                    <strong>Runtime:</strong> {{ archive.completion_info.total_time }}<br>
                    <strong>Throughput:</strong> {{ archive.completion_info.throughput }}
                </div>
                {% endif %}
                {% if archive.progress %}
                <div class="meta-block"><strong>Progress:</strong> {{ archive.progress }}</div>
                {% endif %}
                {% if archive.notes %}
                <div class="meta-block">
                    <strong>Notes:</strong>
                    <div class="notes-preview">{{ archive.notes }}</div>
                </div>
                {% endif %}
                <div class="action-buttons">
                    <button class="restore-btn" onclick="restoreRun('{{ archive.dirname }}', this)">Restore</button>
                </div>
            </div>
        </div>
        {% endfor %}
    </div>
    {% else %}
    <div class="no-archives">
        <h2>No archived runs found</h2>
        <p>Runs moved to <code>archive/</code> will appear here.</p>
    </div>
    {% endif %}

    <div class="toast" id="toast"></div>

    <script>
        const savedTheme = localStorage.getItem('theme') || 'dark';
        document.documentElement.setAttribute('data-theme', savedTheme);

        function restoreRun(dirname, btn) {
            if (!confirm('Restore "' + dirname + '" back to the data directory?')) return;
            btn.disabled = true;
            btn.textContent = 'Restoring…';
            fetch('/restore-run/' + dirname, { method: 'POST' })
                .then(r => r.json())
                .then(data => {
                    if (data.success) {
                        document.getElementById('row-' + dirname).remove();
                        showToast('Restored successfully');
                    } else {
                        btn.disabled = false;
                        btn.textContent = 'Restore';
                        alert('Error: ' + data.error);
                    }
                })
                .catch(err => {
                    btn.disabled = false;
                    btn.textContent = 'Restore';
                    alert('Error: ' + err);
                });
        }

        function showToast(msg) {
            const t = document.getElementById('toast');
            t.textContent = msg;
            t.style.display = 'block';
            setTimeout(() => { t.style.display = 'none'; }, 3000);
        }
    </script>
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
    def format_runtime_display(raw_runtime: str) -> str:
        """Display long runtimes in days with hours in parentheses."""
        match = re.match(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*h\b", raw_runtime)
        if not match:
            return raw_runtime

        hours = float(match.group(1))
        total_minutes = int(round(hours * 60.0))
        total_hours = total_minutes // 60
        minutes = total_minutes % 60
        hours_minutes = f"{total_hours}h {minutes}m"

        if hours <= 48:
            return hours_minutes

        days = total_hours // 24
        rem_hours = total_hours % 24
        return f"{days}d {rem_hours}h ({hours_minutes})"

    completion_file = find_status_file(archive_path, 'sim_completed.txt')
    if not completion_file:
        return None
    
    info = {}
    with open(completion_file) as f:
        for line in f:
            if 'Total time:' in line:
                raw_total_time = line.split(':', 1)[1].strip()
                info['total_time'] = format_runtime_display(raw_total_time)
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
            plot_cache_dir = archive_dir / 'plot_cache'
            has_plot_cache_by_theme = {
                theme: (plot_cache_dir / f'plot1_{theme}.html').exists() or (plot_cache_dir / f'plot2_{theme}.html').exists()
                for theme in PLOT_CACHE_THEMES
            }
            
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
                'params_mismatch': params_mismatch if is_combined else False,
                'has_plot_cache_by_theme': has_plot_cache_by_theme,
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
    
    return render_template_string(
        HTML_TEMPLATE,
        archives=archives,
        rev_archives=rev_archives,
        irr_archives=irr_archives,
        favicon_version=FAVICON_VERSION,
        show_data_transfer_controls=can_show_data_transfer_controls(request),
        show_web_plot_build_controls=can_run_web_plot_build(request),
    )

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

@app.route('/archive-run/<dirname>', methods=['POST'])
def archive_run(dirname):
    """Move a run to the archive folder (non-destructive)."""
    from flask import jsonify
    import shutil

    def next_archive_destination(base_dir: Path, run_name: str) -> Path:
        """Return an available archive destination using +n suffix on collisions."""
        candidate = base_dir / run_name
        if not candidate.exists():
            return candidate

        suffix = 1
        while True:
            candidate = base_dir / f"{run_name}+{suffix}"
            if not candidate.exists():
                return candidate
            suffix += 1
    
    run_path = ARCHIVE_DIR / dirname
    if not run_path.exists():
        return jsonify({'success': False, 'error': 'Run not found'}), 404
    
    if not run_path.is_dir():
        return jsonify({'success': False, 'error': 'Not a directory'}), 400
    
    try:
        DELETED_DIR.mkdir(parents=True, exist_ok=True)
        dest = next_archive_destination(DELETED_DIR, dirname)
        shutil.move(str(run_path), str(dest))
        return jsonify({'success': True, 'archive_name': dest.name})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/restore')
def restore_page():
    """Hidden page listing archived runs available for restoration."""
    archived = []
    if DELETED_DIR.exists():
        for run_dir in sorted(DELETED_DIR.iterdir(), reverse=True):
            if not run_dir.is_dir():
                continue

            dirname = run_dir.name
            try:
                dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
                display_time = dt.strftime('%b %d, %Y %H:%M:%S')
            except ValueError:
                display_time = dirname

            status_info = parse_status_file(run_dir)
            if status_info:
                status = status_info.get('Status', 'UNKNOWN')
                progress = status_info.get('Completed', None)
            else:
                start_file = find_status_file(run_dir, 'sim_started.txt')
                completion_file = find_status_file(run_dir, 'sim_completed.txt')
                status = 'RUNNING' if (start_file and not completion_file) else 'UNKNOWN'
                progress = None

            status_class = status.lower()
            is_combined = (run_dir / 'rev').exists() and (run_dir / 'irr').exists()
            params = parse_start_file(run_dir)
            rev_params = None
            irr_params = None

            if is_combined:
                root_params = parse_start_file_root_first(run_dir)
                if root_params:
                    rev_params = root_params
                    irr_params = root_params
                else:
                    rev_params = parse_start_file(run_dir / 'rev')
                    irr_params = parse_start_file(run_dir / 'irr')

            completion_info = parse_completion_file(run_dir)
            notes = read_notes(run_dir)

            archived.append({
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
            })
    return render_template_string(RESTORE_TEMPLATE, archived=archived, favicon_version=FAVICON_VERSION)

@app.route('/restore-run/<dirname>', methods=['POST'])
def restore_run(dirname):
    """Move a run from the archive folder back to the data directory."""
    from flask import jsonify

    run_path = DELETED_DIR / dirname
    if not run_path.exists():
        return jsonify({'success': False, 'error': 'Run not found in archive'}), 404
    if not run_path.is_dir():
        return jsonify({'success': False, 'error': 'Not a directory'}), 400

    dest = DATA_DIR / dirname
    if dest.exists():
        return jsonify({'success': False, 'error': 'A run with this name already exists in data/'}), 409

    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        shutil.move(str(run_path), str(dest))
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

    theme = normalize_theme(request.args.get('theme'))
    export_mode = normalize_export_mode(request.args.get('export_mode'))
    export_mode_label = 'Plots Only' if export_mode == EXPORT_MODE_PLOTS_ONLY else 'Full Data Set'

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

            function buildPlotQueryString() {{
                const urlParams = new URLSearchParams(window.location.search);
                const keys = ['include_entropy', 'include_zoom', 'include_psd', 'include_summary'];
                const queryParts = [];
                for (const key of keys) {{
                    const value = urlParams.get(key);
                    if (value !== null) {{
                        queryParts.push(key + '=' + value);
                    }} else {{
                        queryParts.push(key + '=1');
                    }}
                }}
                return queryParts.join('&');
            }}
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
            <h1>Exporting Archive ({export_mode_label}) - {title}</h1>
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
                const streamUrl = '/export-stream/{dirname}?theme=' + theme + '&export_mode={export_mode}';

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
    export_mode = normalize_export_mode(request.args.get('export_mode'))

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

            file_list = export_file_list(archive_path, export_mode)
            placeholder_dirs = export_placeholder_dirs(archive_path, export_mode)
            total_files = len(file_list)
            mode_label = 'Plots Only' if export_mode == EXPORT_MODE_PLOTS_ONLY else 'Full Data Set'
            yield f"data: Export mode: {mode_label}\n\n"
            yield f"data: Found {total_files:,} files to export\n\n"

            manifest = {
                'export_version': '1.0',
                'archive_name': dirname,
                'export_mode': export_mode,
                'directories': placeholder_dirs,
                'export_date': datetime.now().isoformat(),
                'files': {}
            }

            # Use stored mode to reduce CPU/memory pressure on small instances.
            with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as zipf:
                for dir_name in placeholder_dirs:
                    write_zip_dir_entry(zipf, dir_name)

                for index, file_path in enumerate(file_list, start=1):
                    rel_path = file_path.relative_to(archive_path)
                    file_size = file_path.stat().st_size

                    sha256_hash = hashlib.sha256()
                    bytes_written = 0

                    # Single-pass read/write/hash so we avoid reading each file twice.
                    with open(file_path, 'rb') as src, zipf.open(str(rel_path), 'w') as dst:
                        for chunk in iter(lambda: src.read(4 * 1024 * 1024), b''):
                            sha256_hash.update(chunk)
                            dst.write(chunk)
                            bytes_written += len(chunk)

                            # Emit periodic keep-alive progress to prevent worker timeout
                            # when a single large file takes a long time to package.
                            if bytes_written % (128 * 1024 * 1024) < len(chunk):
                                pct = int((bytes_written / file_size) * 100) if file_size > 0 else 100
                                yield (
                                    f"data: Packing {index:,}/{total_files:,}: {rel_path} "
                                    f"{bytes_written:,}/{file_size:,} bytes ({pct}%)\n\n"
                                )

                    manifest['files'][str(rel_path)] = {
                        'sha256': sha256_hash.hexdigest(),
                        'size': file_size
                    }

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
                'filename': export_filename(dirname, export_mode)
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
    export_mode = normalize_export_mode(request.args.get('export_mode'))
    if not archive_path.exists():
        return "Archive not found", 404
    
    # Create a temporary zip file
    with tempfile.NamedTemporaryFile(mode='w+b', suffix='.zip', delete=False) as tmp_zip:
        zip_path = Path(tmp_zip.name)
    
    try:
        # Create manifest with file hashes for validation
        placeholder_dirs = export_placeholder_dirs(archive_path, export_mode)

        manifest = {
            'export_version': '1.0',
            'archive_name': dirname,
            'export_mode': export_mode,
            'directories': placeholder_dirs,
            'export_date': datetime.now().isoformat(),
            'files': {}
        }
        
        # Use stored mode to reduce CPU/memory pressure on small instances.
        # Single-pass write keeps IO predictable for large archives.
        with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as zipf:
            for dir_name in placeholder_dirs:
                write_zip_dir_entry(zipf, dir_name)

            for file_path in export_file_list(archive_path, export_mode):
                rel_path = file_path.relative_to(archive_path)
                file_size = file_path.stat().st_size

                # Calculate SHA256 hash of file
                sha256_hash = hashlib.sha256()
                with open(file_path, 'rb') as src, zipf.open(str(rel_path), 'w') as dst:
                    for chunk in iter(lambda: src.read(4 * 1024 * 1024), b''):
                        sha256_hash.update(chunk)
                        dst.write(chunk)

                # Store file and hash in manifest
                manifest['files'][str(rel_path)] = {
                    'sha256': sha256_hash.hexdigest(),
                    'size': file_size
                }
            
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
            download_name=export_filename(dirname, export_mode),
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
                for dir_name in manifest.get('directories', []):
                    dir_path = archive_dir / dir_name
                    if not dir_path.exists() or not dir_path.is_dir():
                        validation_errors.append(f"Missing directory: {dir_name}")

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
                for dir_name in manifest.get('directories', []):
                    (target_path / dir_name).mkdir(parents=True, exist_ok=True)

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
                for dir_name in manifest.get('directories', []):
                    zip_dir_path = subdir_prefix + dir_name.rstrip('/') + '/'
                    if zip_dir_path not in zip_contents:
                        validation_errors.append(f"Missing directory: {dir_name}")

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
                for dir_name in manifest.get('directories', []):
                    (target_path / dir_name).mkdir(parents=True, exist_ok=True)

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
    import time

    theme = normalize_theme(request.args.get('theme'))
    force_replot = request.args.get('force_replot', '0') == '1'

    def load_disk_cache(data_path: Path, active_theme: str, plot_options: dict[str, bool]) -> dict | None:
        paths = plot_cache_paths(data_path, active_theme, plot_options)
        if not paths['dir'].exists():
            return None

        plot1_exists = paths['plot1'].exists()
        plot2_exists = paths['plot2'].exists()
        manifest = None

        if paths['meta'].exists():
            try:
                manifest = json.loads(paths['meta'].read_text())
            except Exception:
                manifest = None

        if not plot1_exists and not plot2_exists:
            return None
        if not plot_options_match_manifest(manifest, plot_options):
            return None

        payload_bytes = 0
        if plot1_exists:
            payload_bytes += paths['plot1'].stat().st_size
        if plot2_exists:
            payload_bytes += paths['plot2'].stat().st_size

        return {
            'meta': manifest,
            'payload_bytes': payload_bytes,
        }

    def format_cache_meta(meta: dict | None) -> str:
        if not meta:
            return "source=unknown built_at=unknown version=unknown"

        return (
            f"source={meta.get('source', 'unknown')} "
            f"built_at={meta.get('built_at', 'unknown')} "
            f"version={meta.get('cache_version', 'unknown')}"
        )

    plot_options = parse_plot_options(request.args)

    if dirname == 'current':
        data_path = DATA_DIR
    else:
        data_path = ARCHIVE_DIR / dirname

    plotly_template = 'plotly_dark' if theme == 'dark' else 'plotly_white'
    plot_script = REPO_ROOT / 'creutz-sim' / 'Sk_comparison.py'

    @stream_with_context
    def generate():
        started_at = time.time()
        web_build_allowed = can_run_web_plot_build(request)

        if not data_path.exists():
            yield "data: ERROR: data path not found\n\n"
            yield "event: error\ndata: data path not found\n\n"
            return
        if not plot_script.exists():
            yield "data: ERROR: plotting script not found\n\n"
            yield "event: error\ndata: plotting script not found\n\n"
            return

        if not force_replot:
            disk_cache = load_disk_cache(data_path, theme, plot_options)
            if disk_cache:
                elapsed = time.time() - started_at
                cache_meta = disk_cache.get('meta')
                payload_bytes = int(disk_cache.get('payload_bytes') or 0)
                app.logger.info(
                    "plot_cache event=cache_hit_disk run=%s theme=%s options=%s elapsed=%.3fs cache_mode=%s payload_mb=%.2f",
                    dirname,
                    theme,
                    plot_options,
                    elapsed,
                    'disk-only',
                    to_mb(payload_bytes),
                )
                yield f"data: cache_hit: loaded cached plot artifacts from disk ({elapsed:.2f}s, disk-only)\n\n"
                yield f"data: cache_meta: {format_cache_meta(cache_meta)}\n\n"
                yield "event: done\ndata: ok\n\n"
                return
            app.logger.info("plot_cache event=cache_miss run=%s theme=%s options=%s", dirname, theme, plot_options)
            yield "data: cache_miss: no cached plot artifacts found for the selected panels, starting cache build\n\n"
        else:
            app.logger.info("plot_cache event=force_replot run=%s theme=%s", dirname, theme)
            yield "data: cache_rebuild: force replot requested, rebuilding cache\n\n"

        if not web_build_allowed:
            app.logger.warning(
                "plot_cache event=cache_build_blocked run=%s theme=%s policy=%s host=%s",
                dirname,
                theme,
                WEB_PLOT_BUILD_POLICY,
                request.host,
            )
            if force_replot:
                yield "data: cache_rebuild_blocked: web replot is disabled for this host\n\n"
            else:
                yield "data: cache_build_blocked: server plotting is disabled for this host\n\n"
            yield "data: hint: run make plot-cache locally, then export/import the run including plot_cache\n\n"
            yield "event: error\ndata: Server plotting disabled for non-localhost requests\n\n"
            return

        yield "data: cache_build_scan: preparing plot build inputs\n\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            with open(plot_script, 'r') as f:
                script_content = f.read()

            # Inject the actual data path via args so no copytree is needed.
            # Pass an absolute path so Path(repo_root) / abs_path resolves to abs_path.
            plot_cli = [
                '--data-dir', str(data_path),
                '--include-entropy' if plot_options['include_entropy'] else '--no-entropy',
                '--include-zoom' if plot_options['include_zoom'] else '--no-zoom',
                '--include-psd' if plot_options['include_psd'] else '--no-psd',
                '--include-summary' if plot_options['include_summary'] else '--no-summary',
            ]
            modified_script = script_content.replace(
                "args = parser.parse_args()",
                f"args = parser.parse_args({plot_cli!r})"
            ).replace(
                'pio.templates.default = "plotly_white"',
                f'pio.templates.default = "{plotly_template}"'
            ).replace(
                'show_plotly_figure(fig, "fig1")',
                f'fig.write_html(r\'{tmpdir_path}/plot1.html\')\n'
                'print("[web-build] wrote plot1.html", flush=True)'
            ).replace(
                'show_plotly_figure(fig2, "fig2")',
                f'fig2.write_html(r\'{tmpdir_path}/plot2.html\')\n'
                'print("[web-build] wrote plot2.html", flush=True)'
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
                app.logger.error("plot_cache event=cache_build_failed run=%s theme=%s reason=subprocess_nonzero", dirname, theme)
                yield "data: cache_build_failed: subprocess exited non-zero, using fallback path failed\n\n"
                yield "event: error\ndata: Script exited with non-zero status\n\n"
                return

            plot1 = tmpdir_path / 'plot1.html'
            plot2 = tmpdir_path / 'plot2.html'

            if not plot1.exists() and not plot2.exists():
                app.logger.error("plot_cache event=cache_build_failed run=%s theme=%s reason=no_plot_html", dirname, theme)
                yield "data: cache_build_failed: no plot html artifacts were produced\n\n"
                yield "event: error\ndata: No plot HTML files generated\n\n"
                return

            yield "data: cache_build_write: writing cache artifacts\n\n"

            plot1_html = plot1.read_text() if plot1.exists() else None
            plot2_html = plot2.read_text() if plot2.exists() else None
            payload_bytes = plot_payload_total_bytes(plot1_html, plot2_html)
            cache_meta = {
                'cache_version': 1,
                'theme': theme,
                'source': 'dynamic-build',
                'built_at': datetime.now().isoformat(),
                'plot_options': plot_options,
            }

            try:
                cache_meta = write_plot_cache_artifacts(
                    data_path,
                    theme,
                    plot1_html,
                    plot2_html,
                    source='dynamic-build',
                    plot_options=plot_options,
                )
            except Exception as e:
                app.logger.error("plot_cache event=cache_build_failed run=%s theme=%s reason=persist_failed error=%s", dirname, theme, str(e))
                yield f"data: cache_build_failed: unable to persist cache to disk ({e})\n\n"

        elapsed = time.time() - started_at
        app.logger.info(
            "plot_cache event=cache_build_done run=%s theme=%s elapsed=%.3fs payload_mb=%.2f",
            dirname,
            theme,
            elapsed,
            to_mb(payload_bytes),
        )
        yield f"data: cache_meta: {format_cache_meta(cache_meta)}\n\n"
        yield f"data: cache_build_done: cache ready and plots rendered ({elapsed:.2f}s)\n\n"

        yield "event: done\ndata: ok\n\n"

    return Response(generate(),
                    mimetype='text/event-stream',
                    headers={'X-Accel-Buffering': 'no', 'Cache-Control': 'no-cache'})


@app.route('/plot-loading/<dirname>')
def plot_loading(dirname):
    """Show loading page while plots are being generated."""
    from flask import request

    # Get theme from query parameter, default to dark
    theme = normalize_theme(request.args.get('theme'))
    plot_options = parse_plot_options(request.args)
    plot_options_bar = plot_options_status_html(plot_options)
    
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

            function buildPlotQueryString() {{
                const urlParams = new URLSearchParams(window.location.search);
                const keys = ['include_entropy', 'include_zoom', 'include_psd', 'include_summary'];
                const queryParts = [];
                for (const key of keys) {{
                    const value = urlParams.get(key);
                    if (value !== null) {{
                        queryParts.push(key + '=' + value);
                    }} else {{
                        queryParts.push(key + '=1');
                    }}
                }}
                return queryParts.join('&');
            }}
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
            
            .plot-option-set {{
                display: flex;
                flex-wrap: wrap;
                justify-content: center;
                gap: 8px 12px;
                margin: 18px auto 0;
                font-size: 0.9em;
                color: var(--text-secondary);
            }}
            .plot-option-chip {{
                display: inline-flex;
                align-items: center;
                gap: 6px;
                padding: 6px 10px;
                border-radius: 999px;
                border: 1px solid var(--border-color);
                background: var(--bg-secondary);
                color: var(--text-primary);
                font-size: 0.85em;
                font-weight: 600;
                letter-spacing: 0.01em;
                user-select: none;
            }}
            .plot-option-chip.included {{
                border-color: rgba(171, 99, 250, 0.65);
                background: rgba(171, 99, 250, 0.12);
                color: var(--text-primary);
            }}
            .plot-option-chip.hidden {{
                opacity: 0.75;
                color: var(--text-secondary);
            }}
            .plot-option-picker {{
                display: inline-flex;
                align-items: center;
                gap: 6px;
                padding: 6px 10px;
                border-radius: 999px;
                border: 1px solid var(--border-color);
                background: var(--bg-secondary);
                cursor: pointer;
            }}
            .plot-option-picker input {{
                accent-color: #ab63fa;
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
            {plot_options_bar}
            <div id="log-box"></div>
            <p class="loading-help">This may take a moment for larger datasets.</p>
        </div>
        <script>
            (function() {{
                const urlParams = new URLSearchParams(window.location.search);
                const theme = urlParams.get('theme') || localStorage.getItem('theme') || 'dark';
                const forceReplot = urlParams.get('force_replot') === '1';
                const baseQuery = buildPlotQueryString();
                const streamUrl = '/plot-stream/{dirname}?theme=' + theme + (baseQuery ? '&' + baseQuery : '') + (forceReplot ? '&force_replot=1' : '');
                const destUrl   = '/plot/{dirname}?theme=' + theme + (baseQuery ? '&' + baseQuery : '');

                const statusEl = document.getElementById('status-text');
                const logBox   = document.getElementById('log-box');
                let cacheMetaLine = '';

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
                    statusEl.textContent = forceReplot ? 'Force refresh requested...' : 'Checking cache...';
                }};

                es.onmessage = function(e) {{
                    const line = e.data;
                    appendLog(line, false);
                    // Update the headline status from key phrases in the output
                    if      (line.includes('cache_hit:'))       statusEl.textContent = 'Using cached plot data...';
                    else if (line.includes('cache_miss:'))      statusEl.textContent = 'Cache miss. Building cache on server...';
                    else if (line.includes('cache_rebuild:'))   statusEl.textContent = 'Rebuilding cache from source data...';
                    else if (line.includes('cache_build_blocked:')) statusEl.textContent = 'Server plotting disabled for this host.';
                    else if (line.includes('cache_rebuild_blocked:')) statusEl.textContent = 'Web replot disabled for this host.';
                    else if (line.includes('cache_build_scan:')) statusEl.textContent = 'Scanning run data...';
                    else if (line.includes('cache_build_write:')) statusEl.textContent = 'Writing cache artifacts...';
                    else if (line.includes('cache_build_done:')) statusEl.textContent = 'Cache ready. Preparing render...';
                    else if (line.includes('cache_build_failed:')) statusEl.textContent = 'Cache build warning. Continuing...';
                    else if (line.includes('cache_meta:'))      cacheMetaLine = line.replace('cache_meta:', '').trim();
                    else if (line.includes('Serializing fig2')) statusEl.textContent = 'Serializing summary plot...';
                    else if (line.includes('Serializing fig1')) statusEl.textContent = 'Serializing main plot...';
                    else if (line.includes('[rev]'))            statusEl.textContent = 'Processing reversible data...';
                    else if (line.includes('[irr]'))            statusEl.textContent = 'Processing irreversible data...';
                }};

                es.addEventListener('done', function() {{
                    es.close();
                    statusEl.textContent = 'Done! Preparing browser rendering...';
                    appendLog('Server finished. Browser is now loading and rendering plots...', false);
                    if (cacheMetaLine) appendLog('Cache metadata: ' + cacheMetaLine, false);
                    // Let the status/log paint before navigating to the heavy plot page.
                    setTimeout(() => window.location.replace(destUrl), 300);
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

    import time

    started_at = time.time()
    theme = normalize_theme(request.args.get('theme'))
    force_replot = request.args.get('force_replot', '0') == '1'
    plot_options = parse_plot_options(request.args)

    if dirname == 'current':
        data_path = DATA_DIR
        if not data_path.exists():
            return "Data directory not found", 404
    else:
        data_path = ARCHIVE_DIR / dirname
        if not data_path.exists():
            return "Archive not found", 404

    cache_dir = data_path / 'plot_cache'
    plot_paths = plot_cache_paths(data_path, theme, plot_options)
    plot1_filename = plot_paths['plot1'].name
    plot2_filename = plot_paths['plot2'].name
    plot1_path = plot_paths['plot1']
    plot2_path = plot_paths['plot2']
    meta_path = plot_paths['meta']
    plot1_exists = plot1_path.exists()
    plot2_exists = plot2_path.exists()
    cache_meta = None
    if meta_path.exists():
        try:
            cache_meta = json.loads(meta_path.read_text())
        except Exception:
            cache_meta = None

    if force_replot:
        if not can_run_web_plot_build(request):
            return (
                "Web replot is disabled for this host. "
                "Run make plot-cache locally and import plot_cache artifacts.",
                403,
            )
        return redirect(url_for('plot_loading', dirname=dirname, theme=theme, force_replot='1', **plot_options))

    if not plot1_exists and not plot2_exists:
        # Not yet generated — send to loading page which drives /plot-stream
        return redirect(url_for('plot_loading', dirname=dirname, theme=theme, **plot_options))
    
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
    
    plot_payload_bytes = 0
    if plot1_exists:
        plot_payload_bytes += plot1_path.stat().st_size
    if plot2_exists:
        plot_payload_bytes += plot2_path.stat().st_size

    cache_meta = cache_meta or {}
    plot1_url = plot_cache_public_url(dirname, plot1_filename, int(plot1_path.stat().st_mtime_ns) if plot1_exists else None) if plot1_exists else None
    plot2_url = plot_cache_public_url(dirname, plot2_filename, int(plot2_path.stat().st_mtime_ns) if plot2_exists else None) if plot2_exists else None

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

            cache_source = cache_meta.get('source', 'session-cache')
            cache_theme = cache_meta.get('theme', theme)
            cache_version = cache_meta.get('cache_version', 'N/A')
            cache_built_at = cache_meta.get('built_at')
            freshness_label = 'UNKNOWN'
            if cache_built_at:
                try:
                    built_dt = datetime.fromisoformat(cache_built_at)
                    age_seconds = max(0.0, (datetime.now() - built_dt).total_seconds())
                    if age_seconds < 86400:
                        freshness_label = 'FRESH'
                    elif age_seconds < 7 * 86400:
                        freshness_label = 'AGING'
                    else:
                        freshness_label = 'STALE'
                    cache_built_at = built_dt.strftime('%Y-%m-%d %H:%M:%S')
                except ValueError:
                    pass
            else:
                cache_built_at = 'Unknown'
            replot_title = html.escape(
                "Rebuild plot cache from source data"
                f" | Status: {freshness_label}"
                f" | Source: {cache_source}"
                f" | Built: {cache_built_at}"
                f" | Theme: {cache_theme}"
                f" | v{cache_version}"
            )
            # Replot is handled from the archive cards, not the plot page itself.
            replot_button_html = ""

            if can_show_data_transfer_controls(request):
                export_menu_html = (
                    '<div class="export-menu-inline">'
                    '<button type="button" class="export-btn" title="Export this simulation archive" aria-haspopup="true" aria-expanded="false" onclick="toggleExportInlineMenu(event, this)">Export</button>'
                    '<div class="export-menu-inline-items">'
                    f'<a href="#" class="export-menu-inline-item" onclick="window.location.href=\'/export-loading/{dirname}?theme=\' + (document.documentElement.getAttribute(\'data-theme\') || localStorage.getItem(\'theme\') || \'dark\') + \'&export_mode=full\'; return false;">Full Data Set</a>'
                    f'<a href="#" class="export-menu-inline-item" onclick="window.location.href=\'/export-loading/{dirname}?theme=\' + (document.documentElement.getAttribute(\'data-theme\') || localStorage.getItem(\'theme\') || \'dark\') + \'&export_mode=plots-only\'; return false;">Plots Only</a>'
                    '</div>'
                    '</div>'
                )
            else:
                export_menu_html = ""
            
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
                        --msg-warning: #ff9800;
                        --msg-warning-hover: #f57c00;
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
                        --msg-warning: #ff9800;
                        --msg-warning-hover: #ef6c00;
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
                    .plot-frame-wrap {{
                        position: relative;
                        min-height: 320px;
                        border: 1px solid var(--border-color);
                        border-radius: 8px;
                        overflow: hidden;
                        background: var(--bg-secondary);
                    }}
                    .plot-frame {{
                        width: 100%;
                        min-height: 320px;
                        border: 0;
                        border-radius: 6px;
                        background: var(--bg-secondary);
                        opacity: 0.08;
                        transition: opacity 0.25s ease, height 0.2s ease;
                    }}
                    .plot-frame-wrap.loaded .plot-frame {{
                        opacity: 1;
                    }}
                    .plot-progress-track {{
                        width: 100%;
                        height: 8px;
                        border-radius: 999px;
                        overflow: hidden;
                        background: rgba(128, 128, 128, 0.28);
                        border: 1px solid var(--border-color);
                    }}
                    .plot-progress-fill {{
                        height: 100%;
                        width: 0%;
                        border-radius: 999px;
                        background: linear-gradient(90deg, #00b8d4 0%, #ab63fa 100%);
                        transition: width 0.6s ease;
                    }}
                    .plot-progress-text {{
                        margin-top: 6px;
                        font-size: 0.82em;
                        color: var(--text-secondary);
                    }}
                    .page-load-curtain {{
                        position: fixed;
                        inset: 0;
                        z-index: 10000;
                        display: flex;
                        align-items: center;
                        justify-content: center;
                        background: rgba(0, 0, 0, 0.62);
                        backdrop-filter: blur(2px);
                        -webkit-backdrop-filter: blur(2px);
                        transition: opacity 0.24s ease;
                    }}
                    .page-load-curtain.hidden {{
                        opacity: 0;
                        pointer-events: none;
                    }}
                    .page-load-card {{
                        width: min(560px, 92vw);
                        background: var(--bg-primary);
                        color: var(--text-primary);
                        border: 1px solid var(--border-color);
                        border-radius: 10px;
                        box-shadow: 0 16px 40px rgba(0,0,0,0.35);
                        padding: 18px 18px 16px;
                    }}
                    .page-load-title {{
                        font-size: 1.06em;
                        font-weight: 700;
                        margin-bottom: 6px;
                    }}
                    .page-load-status {{
                        color: var(--text-secondary);
                        font-size: 0.93em;
                        margin-bottom: 10px;
                    }}
                    .page-load-hint {{
                        margin-top: 8px;
                        font-size: 0.84em;
                        color: var(--text-secondary);
                    }}
                    .back-link {{ 
                        margin-bottom: 20px;
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                    }}
                    .back-link > a {{
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
                    .back-link > a:hover {{ 
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
                    .export-menu-inline {{
                        position: relative;
                        display: inline-block;
                    }}
                    .export-menu-inline .export-btn {{
                        cursor: pointer;
                        appearance: none;
                        -webkit-appearance: none;
                    }}
                    .export-menu-inline .export-btn:focus-visible {{
                        outline: 2px solid #28a745;
                        outline-offset: 2px;
                    }}
                    .export-menu-inline-items {{
                        position: absolute;
                        right: 0;
                        top: calc(100% + 6px);
                        min-width: 170px;
                        background: var(--bg-secondary);
                        border: 1px solid var(--border-color);
                        border-radius: 6px;
                        box-shadow: 0 4px 10px rgba(0,0,0,0.25);
                        z-index: 30;
                        padding: 6px;
                        display: none;
                    }}
                    .export-menu-inline.open .export-menu-inline-items {{
                        display: block;
                    }}
                    .export-menu-inline-item {{
                        display: block;
                        padding: 7px 9px;
                        border-radius: 4px;
                        color: var(--text-primary);
                        text-decoration: none;
                        font-size: 0.9em;
                        font-weight: 500;
                        line-height: 1.3;
                    }}
                    .export-menu-inline-item:hover {{
                        background: var(--param-bg);
                        text-decoration: none;
                    }}
                    .replot-btn {{
                        background: var(--msg-warning);
                        color: white;
                        border: none;
                        padding: 8px 16px;
                        border-radius: 4px;
                        cursor: pointer;
                        font-size: 14px;
                        font-weight: 600;
                        transition: all 0.3s ease;
                        display: flex;
                        align-items: center;
                        height: 38px;
                    }}
                    .replot-btn:hover {{
                        background: var(--msg-warning-hover);
                        transform: scale(1.05);
                    }}
                    .replot-btn-disabled,
                    .replot-btn-disabled:hover {{
                        background: #888;
                        cursor: not-allowed;
                        transform: none;
                        opacity: 0.75;
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

                    function initPlotFrameLoading() {{
                        const wraps = document.querySelectorAll('.plot-frame-wrap');
                        const autoWrapsForCurtain = Array.from(document.querySelectorAll('.plot-frame-wrap[data-autoload="1"]'));
                        const curtainEl = document.getElementById('pageLoadCurtain');
                        const curtainStatusEl = document.getElementById('pageLoadStatus');

                        function updateCurtainState() {{
                            if (!curtainEl || autoWrapsForCurtain.length === 0) return;

                            let doneCount = 0;
                            autoWrapsForCurtain.forEach((wrap) => {{
                                if (wrap.classList.contains('loaded') || wrap.classList.contains('error')) doneCount += 1;
                            }});

                            if (curtainStatusEl) {{
                                if (doneCount < autoWrapsForCurtain.length) {{
                                    curtainStatusEl.textContent = `Loading plots: ${{doneCount}}/${{autoWrapsForCurtain.length}} complete`;
                                }} else {{
                                    const failed = autoWrapsForCurtain.filter(w => w.classList.contains('error')).length;
                                    curtainStatusEl.textContent = failed > 0
                                        ? `Loaded with warnings: ${{autoWrapsForCurtain.length - failed}} ok, ${{failed}} failed`
                                        : 'All plots loaded';
                                }}
                            }}

                            if (doneCount === autoWrapsForCurtain.length) {{
                                setTimeout(() => curtainEl.classList.add('hidden'), 380);
                            }}
                        }}

                        wraps.forEach((wrap, index) => {{
                            const iframe = wrap.querySelector('iframe');
                            if (!iframe) return;

                            function resizeToContent() {{
                                try {{
                                    const doc = iframe.contentWindow && iframe.contentWindow.document;
                                    const contentHeight = doc && (doc.body ? doc.body.scrollHeight : 0);
                                    if (contentHeight) {{
                                        iframe.style.height = contentHeight + 'px';
                                    }}
                                }} catch (e) {{
                                    // Cross-origin or not-yet-ready iframe content; keep the CSS fallback height.
                                }}
                            }}

                            iframe.addEventListener('load', () => {{
                                wrap.classList.remove('error');
                                wrap.classList.add('loaded');
                                resizeToContent();
                                updateCurtainState();
                            }});

                            iframe.addEventListener('error', () => {{
                                wrap.classList.add('error');
                                updateCurtainState();
                            }});
                        }});

                        updateCurtainState();
                    }}

                    function closeAllInlineExportMenus(exceptMenu = null) {{
                        document.querySelectorAll('.export-menu-inline').forEach((menu) => {{
                            if (menu === exceptMenu) return;
                            menu.classList.remove('open');
                            const trigger = menu.querySelector('.export-btn');
                            if (trigger) trigger.setAttribute('aria-expanded', 'false');
                        }});
                    }}

                    function toggleExportInlineMenu(event, triggerButton) {{
                        event.preventDefault();
                        event.stopPropagation();
                        const menu = triggerButton.closest('.export-menu-inline');
                        if (!menu) return;
                        const isOpen = menu.classList.contains('open');
                        closeAllInlineExportMenus(menu);
                        menu.classList.toggle('open', !isOpen);
                        triggerButton.setAttribute('aria-expanded', String(!isOpen));
                    }}
                    
                    let originalNotes = '';
                    let isEditMode = false;
                    
                    document.addEventListener('DOMContentLoaded', function() {{
                        initTheme();  // Initialize theme on page load
                        initPlotFrameLoading();
                        
                        const textarea = document.getElementById('notesTextarea');
                        const saveBtn = document.getElementById('saveNotesBtn');
                        const saveStatus = document.getElementById('saveStatus');
                        if (!textarea || !saveBtn || !saveStatus) return;
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

                    document.addEventListener('click', function(event) {{
                        if (!event.target.closest('.export-menu-inline')) {{
                            closeAllInlineExportMenus();
                        }}
                    }});

                    document.addEventListener('keydown', function(event) {{
                        if (event.key === 'Escape') {{
                            closeAllInlineExportMenus();
                        }}
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
                <div id="pageLoadCurtain" class="page-load-curtain" aria-live="polite">
                    <div class="page-load-card">
                        <div class="page-load-title">Loading plots...</div>
                        <div id="pageLoadStatus" class="page-load-status">Preparing...</div>
                        <div class="page-load-hint">Please wait while both plots render. Large data sets can take several minutes to load.</div>
                    </div>
                </div>
                <div class="back-link">
                    <a href="/" title="Back to simulation browser">‹</a>
                    <div class="top-controls">
                        {replot_button_html}
                        {export_menu_html}
                        <button class="refresh-btn" onclick="location.reload()" title="Refresh the page to see latest data">🔄</button>
                    </div>
                </div>
                <h1>Simulation Plots - {title}</h1>
                {params_html}
                {notes_html}
            """
            
            if plot1_url:
                html_content += (
                    '<div class="plot-container">'
                    '<div class="plot-frame-wrap" data-plot-name="Plot 1" data-autoload="1">'
                    f'<iframe class="plot-frame" src="{plot1_url}" loading="eager" title="Plot 1" '
                    'referrerpolicy="same-origin"></iframe>'
                    '</div>'
                    '</div>'
                )
            
            if plot2_url:
                html_content += (
                    '<div class="plot-container">'
                    '<div class="plot-frame-wrap" data-plot-name="Plot 2" data-autoload="1">'
                    f'<iframe class="plot-frame" src="{plot2_url}" loading="lazy" title="Plot 2" '
                    'referrerpolicy="same-origin"></iframe>'
                    '</div>'
                    '</div>'
                )
            
            html_content += """
            </body>
            </html>
            """
            
            # Add cache control headers to prevent browser caching
            response = make_response(html_content)
            response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
            response.headers['Pragma'] = 'no-cache'
            response.headers['Expires'] = '0'
            response_bytes = len(html_content.encode('utf-8'))
            elapsed_ms = (time.time() - started_at) * 1000.0
            app.logger.info(
                "plot_serve event=plot_page_served run=%s theme=%s source=%s payload_mb=%.2f response_mb=%.2f elapsed_ms=%.1f",
                dirname,
                theme,
                cache_meta.get('source', 'unknown'),
                to_mb(plot_payload_bytes),
                to_mb(response_bytes),
                elapsed_ms,
            )
            return response

    except Exception as e:
        import traceback
        return f"<h1>Error generating plots</h1><pre>{traceback.format_exc()}</pre>", 500


@app.route('/plot-cache/<dirname>/plot_cache/<plot_filename>')
def serve_plot_cache_file(dirname, plot_filename):
    """Flask fallback for cached plot artifacts.

    In production, nginx should serve this path directly from disk.
    """
    allowed_files = {f'plot1_{theme}.html' for theme in PLOT_CACHE_THEMES}
    allowed_files.update({f'plot2_{theme}.html' for theme in PLOT_CACHE_THEMES})

    if plot_filename not in allowed_files:
        return "Invalid plot cache file", 400

    data_path = resolve_plot_data_path(dirname)
    if not data_path.exists():
        return "Run not found", 404

    file_path = data_path / 'plot_cache' / plot_filename
    if not file_path.exists() or not file_path.is_file():
        return "Plot cache file not found", 404

    response = send_file(file_path, mimetype='text/html', conditional=True)
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    return response

if __name__ == '__main__':
    print("\n" + "="*60)
    print("  Simulation Browser")
    print("="*60)
    print(f"\n  📂 Data directory: {DATA_DIR}")
    print(f"  📦 Archive directory: {ARCHIVE_DIR}")
    print(f"\n  🌐 Open in browser: http://127.0.0.1:5001")
    print("\n  Press Ctrl+C to stop\n")
    
    app.run(debug=True, port=5001, host='127.0.0.1')
