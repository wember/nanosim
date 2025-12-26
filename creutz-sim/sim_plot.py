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

# Get CSV file from command line or look for simulation output
if len(sys.argv) > 1:
    csv_file = sys.argv[1]
else:
    # Look for simulation CSV files (not example pipelines)
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

try:
    df = pd.read_csv(csv_file)
    
    # Validate required columns
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

n = df['n'][0]  # Lattice size
t = df['t']    # step number
U = df['U']    # potential (lattice) energy
K = df['K']    # kinetic (demon) energy
Nx = df['Nx']   # anti-aligned spins
N0 = df['N0']   # broken bonds

df['N0_exp'] = df['N0'].replace(0, 1) # If no broken bonds, 2^N0+1 in Su equation
N0_exp = df['N0_exp']

Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
Su = lambda N, N0, Nx: logg(N+1) + np.log(2**N0_exp) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins

total_entropy = (Sk(n, K) + Su(n, N0, Nx))/n

fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

### Demon Energy ###
fig.add_trace(go.Scatter(x=t, y=K, showlegend=False),row=1, col=1)
fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="Demon Energy", row=1, col=1)
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)

### Lattice Energy ###
fig.add_trace(go.Scatter(x=t, y=U, showlegend=False),row=1, col=2)
# fig.add_trace(go.Scatter(x=tR, y=KR, showlegend=False),row=1, col=2)
fig.update_xaxes(title_text="Sweeps", row=1, col=2)
fig.update_yaxes(title_text="Lattice Temp", row=1, col=2)
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)

### Entropy ###
fig.add_trace(go.Scatter(x=t, y=total_entropy, showlegend=False),row=1, col=3)
fig.update_xaxes(title_text="Sweeps", row=1, col=3)
fig.update_yaxes(title_text="S/Nk", row=1, col=3)
fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=3)

fig.update_layout(title_text=f"Lattice Size: {n}")
fig.show()
