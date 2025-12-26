"""Single Simulation Visualization

Generates a 1x3 subplot layout showing:
- Demon energy vs sweeps
- Lattice temperature (energy) vs sweeps  
- Total entropy (S/Nk) vs sweeps

Red dashed line marks the midpoint where forward phase ends and reverse begins.

Usage:
  python creutz-sim/sim_plot.py [path/to/sim_data.csv]
  
If no file specified, automatically finds the most recent simulation CSV in:
- data/r{0-10}/sim_data*.csv
- data/irr/r{0-10}/sim_data*.csv

Alternatively, run via Makefile:
  make plot
"""

import pandas as pd
import numpy as np
from scipy.special import loggamma as logg
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import sys
import os
import glob
pio.templates.default = "plotly_white"

# =============================================================================
# CSV File Discovery
# =============================================================================
# Get CSV file from command line argument or auto-discover most recent file
if len(sys.argv) > 1:
    csv_file = sys.argv[1]
else:
    # Auto-discovery: search common data directories for simulation output
    project_root = os.path.dirname(os.path.dirname(__file__))
    
    # Check common locations for simulation output
    search_paths = [
        os.path.join(project_root, 'data', 'r*', 'sim_data*.csv'),
        os.path.join(project_root, 'data', 'irr', 'r*', 'sim_data*.csv'),
    ]
    
    csv_files = []
    for pattern in search_paths:
        csv_files.extend(glob.glob(pattern))
    
    if csv_files:
        # Use most recently modified file
        csv_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        csv_file = csv_files[0]
        print(f"Using most recent simulation CSV: {csv_file}")
    else:
        print("Error: No simulation CSV files found")
        print("Run a simulation first to generate data:")
        print("  make run-sim-small")
        print("")
        print("Or specify a CSV file:")
        print("  python creutz-sim/sim_plot.py path/to/sim_data.csv")
        sys.exit(1)

# =============================================================================
# Load and Validate CSV Data
# =============================================================================
try:
    df = pd.read_csv(csv_file)
    
    # Validate CSV has required columns from simulation output
    # (distinguishes from analysis pipeline CSVs which have different structure)
    required_cols = ['t', 'K', 'U', 'N0', 'Nx', 'n']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: CSV file missing required columns: {missing_cols}")
        print(f"Found columns: {list(df.columns)}")
        print("")
        print("This script requires simulation output CSV, not analysis pipeline CSV.")
        print("Run: make run-sim-small")
        sys.exit(1)
        
except FileNotFoundError:
    print(f"Error: File not found: {csv_file}")
    sys.exit(1)

# =============================================================================
# Extract Data Columns
# =============================================================================
n = df['n'][0]  # Lattice size
t = df['t']     # Sweep number (time step)
U = df['U']     # Lattice energy per site
K = df['K']     # Demon energy per site (kinetic)
Nx = df['Nx']   # Anti-aligned neighbor pairs per site
N0 = df['N0']   # Broken bonds per site

# Handle N0=0 edge case: use 2^(N0+1) instead of 2^N0 to avoid log(0)
df['N0_exp'] = df['N0'].replace(0, 1)
N0_exp = df['N0_exp']

# =============================================================================
# Entropy Calculation (High-Precision)
# =============================================================================
# Sk: Kinetic entropy from demon energy distribution
# Uses loggamma for numerical stability with large factorials
Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N)

# Su: Configurational entropy from spin/bond states  
# Uses N0_exp * log(2) instead of log(2^N0_exp) to prevent overflow
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**N0_exp) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

# Total entropy per site (in units of Boltzmann constant)
total_entropy = (Sk(n, K) + Su(n, N0, Nx))/n

# =============================================================================
# Create 1x3 Subplot Layout
# =============================================================================
fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

# =============================================================================
# Plot 1: Demon Energy vs Time
# =============================================================================
fig.add_trace(go.Scatter(x=t, y=K, showlegend=False),row=1, col=1)
fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="Demon Energy", row=1, col=1)
# Red line marks transition from forward to reverse phase
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)

# =============================================================================
# Plot 2: Lattice Temperature (Energy) vs Time
# =============================================================================
fig.add_trace(go.Scatter(x=t, y=U, showlegend=False),row=1, col=2)
fig.update_xaxes(title_text="Sweeps", row=1, col=2)
fig.update_yaxes(title_text="Lattice Temp", row=1, col=2)
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)

# =============================================================================
# Plot 3: Total Entropy vs Time
# =============================================================================
fig.add_trace(go.Scatter(x=t, y=total_entropy, showlegend=False),row=1, col=3)
fig.update_xaxes(title_text="Sweeps", row=1, col=3)
fig.update_yaxes(title_text="S/Nk", row=1, col=3)
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=3)

# Display interactive plot
fig.update_layout(title_text=f"Single Simulation Results (Lattice Size: {n})")
fig.show()
