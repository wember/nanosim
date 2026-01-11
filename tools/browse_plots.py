#!/usr/bin/env venv/bin/python
"""Simple web interface to browse current and archived simulation runs."""

from flask import Flask, render_template_string, send_file
from pathlib import Path
from datetime import datetime
import json

app = Flask(__name__)

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / 'data'
ARCHIVE_DIR = DATA_DIR / 'archive'

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Simulation Browser</title>
    <style>
        body { 
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            max-width: 1200px; 
            margin: 40px auto; 
            padding: 0 20px;
            background: #f5f5f5;
        }
        h1 { 
            color: #333;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .refresh-btn {
            background: #007bff;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            display: flex;
            align-items: center;
            gap: 6px;
        }
        .refresh-btn:hover {
            background: #0056b3;
        }
        .archive-list { 
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .archive-item { 
            padding: 20px;
            border-bottom: 1px solid #eee;
        }
        .archive-item:last-child { border-bottom: none; }
        .archive-item:hover { background: #fafafa; }
        .timestamp { 
            font-size: 1.2em; 
            font-weight: bold; 
            color: #333;
            margin-bottom: 8px;
        }
        .status { 
            display: inline-block;
            padding: 4px 12px;
            border-radius: 4px;
            font-weight: 600;
            font-size: 0.85em;
            margin-bottom: 8px;
        }
        .status-completed { background: #d4edda; color: #155724; }
        .status-interrupted { background: #fff3cd; color: #856404; }
        .status-error { background: #f8d7da; color: #721c24; }
        .status-unknown { background: #e2e3e5; color: #383d41; }
        .status-current { background: #cfe2ff; color: #084298; }
        .details { 
            color: #666; 
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
            background: #f8f9fa;
            padding: 4px 8px;
            border-radius: 4px;
        }
        .no-archives {
            text-align: center;
            padding: 60px 20px;
            color: #999;
        }
        a { color: #007bff; text-decoration: none; }
        a:hover { text-decoration: underline; }
        
        /* Modal styles */
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.4);
        }
        .modal-content {
            background-color: #fefefe;
            margin: 10% auto;
            padding: 20px;
            border: 1px solid #888;
            border-radius: 8px;
            width: 80%;
            max-width: 600px;
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
            border: 1px solid #ddd;
            border-radius: 4px;
            font-family: inherit;
            font-size: 14px;
            resize: vertical;
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
    </style>
    <script>
        let currentDirname = null;
        
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
            const modal = document.getElementById('notesModal');
            if (event.target == modal) {
                closeNotesModal();
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
                        alert('Archive deleted successfully');
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
    
    <h1>
        <span>📊 Simulation Browser</span>
        <button class="refresh-btn" onclick="window.location.reload()">🔄 Refresh</button>
    </h1>
    
    {% if archives %}
    <div class="archive-list">
        {% for archive in archives %}
        <div class="archive-item">
            <div class="timestamp">{{ archive.display_time }}</div>
            <span class="status status-{{ archive.status_class }}">{{ archive.status }}</span>
            <div class="details">
                {% if archive.params %}
                <div class="details-grid">
                    <div><strong>Lattice:</strong> <span class="param">n={{ archive.params.n }}</span></div>
                    <div><strong>Sweeps:</strong> <span class="param">s={{ archive.params.sweeps }}</span></div>
                    <div><strong>Flag:</strong> <span class="param">{{ archive.params.flag }}</span></div>
                    <div><strong>Radius:</strong> <span class="param">r={{ archive.params.radius }}</span></div>
                    <div><strong>Runs:</strong> <span class="param">m={{ archive.params.runs }}</span></div>
                    <div><strong>Total sims:</strong> <span class="param">{{ archive.params.total }}</span></div>
                </div>
                {% endif %}
                {% if archive.completion_info %}
                <div style="margin-top: 8px;">
                    <strong>Runtime:</strong> {{ archive.completion_info.total_time }}<br>
                    <strong>Throughput:</strong> {{ archive.completion_info.throughput }}
                </div>
                {% endif %}
                {% if archive.progress %}
                <div style="margin-top: 8px;">
                    <strong>Progress:</strong> {{ archive.progress }}
                </div>
                {% endif %}
                {% if archive.notes %}
                <div style="margin-top: 8px;">
                    <strong>Notes:</strong> 
                    <a href="#" data-dirname="{{ archive.dirname }}" data-notes="{{ archive.notes|escape }}" onclick="openNotesModalFromLink(this); return false;" style="font-size: 0.9em;">📝 Edit</a>
                    <div style="font-style: italic; white-space: pre-wrap; margin-top: 4px;">{{ archive.notes }}</div>
                </div>
                {% endif %}
                <div style="margin-top: 12px;">
                    <a href="/plot/{{ archive.dirname }}">📈 Plot data</a>
                    <span style="color: #ccc; margin: 0 8px;">|</span>
                    <a href="/view/{{ archive.dirname }}">📂 View files</a>
                    {% if not archive.notes %}
                    <span style="color: #ccc; margin: 0 8px;">|</span>
                    <a href="#" data-dirname="{{ archive.dirname }}" data-notes="" onclick="openNotesModalFromLink(this); return false;">📝 Add notes</a>
                    {% endif %}
                    {% if archive.dirname != 'current' %}
                    <span style="color: #ccc; margin: 0 8px;">|</span>
                    <a href="#" onclick="deleteArchive('{{ archive.dirname }}'); return false;" style="color: #dc3545;">🗑️ Delete</a>
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

FILE_LIST_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Archive Files - {{ dirname }}</title>
    <style>
        body { 
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
            max-width: 1200px; 
            margin: 40px auto; 
            padding: 0 20px;
        }
        h1 { color: #333; }
        .file-list { list-style: none; padding: 0; }
        .file-item { 
            padding: 12px;
            border-bottom: 1px solid #eee;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .file-item:hover { background: #f8f9fa; }
        .file-name { font-family: 'Courier New', monospace; }
        .file-size { color: #666; font-size: 0.9em; }
        a { color: #007bff; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .back-link { margin-bottom: 20px; }
    </style>
</head>
<body>
    <div class="back-link"><a href="/">← Back to browser</a></div>
    <h1>📂 Archive: {{ dirname }}</h1>
    <ul class="file-list">
        {% for file in files %}
        <li class="file-item">
            <span class="file-name">{{ file.path }}</span>
            <span class="file-size">{{ file.size }}</span>
        </li>
        {% endfor %}
    </ul>
</body>
</html>
"""

def parse_status_file(archive_path):
    """Parse sim_status.txt to get status info."""
    status_file = archive_path / 'sim_status.txt'
    if not status_file.exists():
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
    start_file = archive_path / 'sim_started.txt'
    if not start_file.exists():
        return None
    
    params = {}
    with open(start_file) as f:
        for line in f:
            if 'Parameters:' in line:
                # Parse: n=10, sweeps=100, flag=0, radius=11, runs=5
                param_str = line.split('Parameters:')[1].strip()
                for param in param_str.split(','):
                    key, value = param.strip().split('=')
                    params[key] = value
            elif 'Total simulations:' in line:
                params['total'] = line.split(':')[1].strip()
    return params

def parse_completion_file(archive_path):
    """Parse sim_completed.txt to get completion info."""
    completion_file = archive_path / 'sim_completed.txt'
    if not completion_file.exists():
        return None
    
    info = {}
    with open(completion_file) as f:
        for line in f:
            if 'Total time:' in line:
                info['total_time'] = line.split(':', 1)[1].strip()
            elif 'Throughput:' in line:
                info['throughput'] = line.split(':', 1)[1].strip()
    return info

def read_notes(archive_path):
    """Read sim_notes.txt to get user notes."""
    notes_file = archive_path / 'sim_notes.txt'
    if not notes_file.exists():
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
    
    # Add current data directory if it exists
    if DATA_DIR.exists() and any(DATA_DIR.iterdir()):
        # Get status
        status_info = parse_status_file(DATA_DIR)
        if status_info:
            status = status_info.get('Status', 'UNKNOWN')
            progress = status_info.get('Completed', None)
        else:
            status = 'IN PROGRESS' if (DATA_DIR / 'sim_started.txt').exists() else 'UNKNOWN'
            progress = None
        
        # Get parameters
        params = parse_start_file(DATA_DIR)
        
        # Get completion info
        completion_info = parse_completion_file(DATA_DIR)
        
        # Get notes
        notes = read_notes(DATA_DIR)
        
        # Get modification time as display time
        try:
            status_file = DATA_DIR / 'sim_status.txt'
            if status_file.exists():
                mtime = datetime.fromtimestamp(status_file.stat().st_mtime)
            elif (DATA_DIR / 'sim_started.txt').exists():
                mtime = datetime.fromtimestamp((DATA_DIR / 'sim_started.txt').stat().st_mtime)
            else:
                mtime = datetime.fromtimestamp(DATA_DIR.stat().st_mtime)
            display_time = f"Current Run (updated {mtime.strftime('%H:%M:%S')})"
        except:
            display_time = "Current Run"
        
        archives.append({
            'dirname': 'current',
            'display_time': display_time,
            'status': status,
            'status_class': 'current',
            'params': params,
            'progress': progress,
            'completion_info': completion_info,
            'notes': notes
        })
    
    # Add archived runs
    if ARCHIVE_DIR.exists():
        for archive_dir in sorted(ARCHIVE_DIR.iterdir(), reverse=True):
            if not archive_dir.is_dir():
                continue
            
            # Parse timestamp from directory name (YYYYMMDD_HHMMSS)
            dirname = archive_dir.name
            try:
                dt = datetime.strptime(dirname, '%Y%m%d_%H%M%S')
                display_time = dt.strftime('%Y-%m-%d %H:%M:%S')
            except ValueError:
                display_time = dirname
            
            # Get status
            status_info = parse_status_file(archive_dir)
            if status_info:
                status = status_info.get('Status', 'UNKNOWN')
                progress = status_info.get('Completed', None)
            else:
                status = 'UNKNOWN'
                progress = None
            
            # Map status to CSS class
            status_class = status.lower()
            
            # Get parameters
            params = parse_start_file(archive_dir)
            
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
                'notes': notes
            })
    
    return render_template_string(HTML_TEMPLATE, archives=archives)

@app.route('/notes/<dirname>', methods=['POST'])
def update_notes(dirname):
    """Update notes for a run."""
    from flask import jsonify, request
    
    if dirname == 'current':
        notes_path = DATA_DIR / 'sim_notes.txt'
    else:
        archive_path = ARCHIVE_DIR / dirname
        if not archive_path.exists():
            return jsonify({'success': False, 'error': 'Archive not found'}), 404
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

@app.route('/delete/<dirname>', methods=['POST'])
def delete_archive(dirname):
    """Delete an archived run."""
    from flask import jsonify
    import shutil
    
    # Don't allow deleting current data directory
    if dirname == 'current':
        return jsonify({'success': False, 'error': 'Cannot delete current data directory'}), 400
    
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
    
    return render_template_string(FILE_LIST_TEMPLATE, dirname=dirname, files=files)

@app.route('/plot/<dirname>')
def plot_data(dirname):
    """Generate plots for a simulation run."""
    import subprocess
    import tempfile
    import shutil
    
    if dirname == 'current':
        data_path = DATA_DIR
        if not data_path.exists():
            return "Data directory not found", 404
    else:
        data_path = ARCHIVE_DIR / dirname
        if not data_path.exists():
            return "Archive not found", 404
    
    # Create a temporary directory to work in
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Copy data to temp location
        temp_data = tmpdir_path / 'data'
        shutil.copytree(data_path, temp_data)
        
        # Copy the plotting script
        plot_script = REPO_ROOT / 'creutz-sim' / 'Sk_comparison.py'
        if not plot_script.exists():
            return "Plotting script not found", 404
        
        temp_script = tmpdir_path / 'plot_script.py'
        
        # Read the original script and modify it to save HTML instead of showing
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Replace the data path references and change .show() to .write_html()
        modified_script = script_content.replace(
            "repo_root = Path(__file__).parent.parent",
            f"repo_root = Path('{tmpdir_path}')"
        ).replace(
            "filepath = repo_root / 'data'",
            f"filepath = Path('{temp_data}')"
        ).replace(
            "fig.show()",
            f"fig.write_html('{tmpdir_path}/plot1.html')"
        ).replace(
            "fig2.show()",
            f"fig2.write_html('{tmpdir_path}/plot2.html')"
        )
        
        with open(temp_script, 'w') as f:
            f.write(modified_script)
        
        # Run the plotting script
        try:
            result = subprocess.run(
                [REPO_ROOT / 'venv' / 'bin' / 'python', str(temp_script)],
                cwd=tmpdir_path,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode != 0:
                return f"<h1>Plot Generation Failed</h1><pre>stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}</pre>", 500
            
            # Find the generated HTML files
            plot1 = tmpdir_path / 'plot1.html'
            plot2 = tmpdir_path / 'plot2.html'
            
            if not plot1.exists() and not plot2.exists():
                return f"<h1>No plots generated</h1><pre>stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}</pre>", 500
            
            # Combine plots into a single page
            title = 'Current Run' if dirname == 'current' else dirname
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>Simulation Plots</title>
                <style>
                    body {{ margin: 0; padding: 20px; font-family: Arial, sans-serif; }}
                    h1 {{ text-align: center; }}
                    .plot-container {{ margin: 20px 0; }}
                    .back-link {{ margin-bottom: 20px; }}
                    .back-link a {{ color: #007bff; text-decoration: none; font-size: 16px; }}
                    .back-link a:hover {{ text-decoration: underline; }}
                </style>
            </head>
            <body>
                <div class="back-link"><a href="/">← Back to browser</a></div>
                <h1>Simulation Plots - {title}</h1>
            """
            
            if plot1.exists():
                with open(plot1, 'r') as f:
                    html_content += f'<div class="plot-container">{f.read()}</div>'
            
            if plot2.exists():
                with open(plot2, 'r') as f:
                    html_content += f'<div class="plot-container">{f.read()}</div>'
            
            html_content += """
            </body>
            </html>
            """
            
            return html_content
                
        except subprocess.TimeoutExpired:
            return "Plot generation timed out (>30s)", 500
        except Exception as e:
            import traceback
            return f"<h1>Error generating plots</h1><pre>{traceback.format_exc()}</pre>", 500

if __name__ == '__main__':
    print("\n" + "="*60)
    print("  Simulation Browser")
    print("="*60)
    print(f"\n  📂 Data directory: {DATA_DIR}")
    print(f"  📦 Archive directory: {ARCHIVE_DIR}")
    print(f"\n  🌐 Open in browser: http://127.0.0.1:5000")
    print("\n  Press Ctrl+C to stop\n")
    
    app.run(debug=True, port=5000, host='127.0.0.1')
