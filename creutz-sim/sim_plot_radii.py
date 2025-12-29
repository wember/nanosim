"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Multi-Radius Comparison Visualization

Compares simulation results across different demon-coupling radii (R=0 to R=10)
for a single simulation type (reversible or irreversible).

Generates a 1x3 subplot layout showing all radii overlaid:
- Demon energy vs sweeps
- Lattice temperature (energy) vs sweeps
- Total entropy (S/Nk) vs sweeps

Red dashed line marks the midpoint where forward phase ends and reverse begins.

Automatically searches for CSV files in data/r{0-10}/ directories relative to
the project root.

Usage:
  python creutz-sim/sim_plot_radii.py
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.special import loggamma as logg

pio.templates.default = "plotly_white"

# =============================================================================
# Create 1x3 Subplot Layout
# =============================================================================
fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

# =============================================================================
# File Discovery
# =============================================================================
# Auto-discover CSV files from data/r{0-10}/ directories using project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Build list of CSV files and display names
file_names = []
data_names = []

for R in range(11):  # R=0 to R=10
    folder_path = os.path.join(project_root, "data", f"r{R}")
    # Look for any sim_data CSV files in this radius directory
    csv_pattern = os.path.join(folder_path, "sim_data*.csv")
    csv_files = glob.glob(csv_pattern)

    if csv_files:
        # Use the first CSV found (or most recent if multiple)
        csv_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        file_names.append(csv_files[0])
        data_names.append(f"radius = {R}")
    else:
        print(f"Warning: No CSV files found in {folder_path}")

# Check if we found any files
if not file_names:
    print("Error: No simulation CSV files found in data/r{0-10}/ directories")
    print("Run a simulation first to generate data:")
    print("  make run-sim")
    sys.exit(1)

print(f"Found {len(file_names)} CSV files to plot")

# =============================================================================
# Process Each Radius
# =============================================================================
for idx, file in enumerate(file_names):
    df = pd.read_csv(file)

    # Extract data columns
    n = df["n"][0]  # Lattice size
    t = df["t"]  # Sweep number (time step)
    U = df["U"]  # Lattice energy per site
    K = df["K"]  # Demon energy per site (kinetic)
    Nx = df["Nx"]  # Anti-aligned neighbor pairs per site
    N0 = df["N0"]  # Broken bonds per site

    # Handle N0=0 edge case: use 2^(N0+1) instead of 2^N0 to avoid log(0)
    df["N0_exp"] = df["N0"].replace(0, 1)
    N0_exp = df["N0_exp"]

    # Entropy calculation (high-precision using loggamma)
    Sk = lambda N, K: logg(K + N) - logg(K + 1) - logg(N)  # Kinetic entropy
    Su = (
        lambda N, N0, Nx: logg(N + 1)
        + np.log(2**N0_exp)
        - (logg(N - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))
    )  # Configurational entropy

    # Total entropy per site (in units of Boltzmann constant)
    total_entropy = (Sk(n, K) + Su(n, N0, Nx)) / n

    # Plot 1: Demon Energy vs Time (all radii overlaid)
    fig.add_trace(go.Scatter(x=t, y=K, name=data_names[idx]), row=1, col=1)
    fig.update_xaxes(title_text="Sweeps", row=1, col=1)
    fig.update_yaxes(title_text="Demon Energy", row=1, col=1)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1
    )

    # Plot 2: Lattice Temperature vs Time (all radii overlaid)
    fig.add_trace(go.Scatter(x=t, y=U, name=data_names[idx]), row=1, col=2)
    fig.update_xaxes(title_text="Sweeps", row=1, col=2)
    fig.update_yaxes(title_text="Lattice Temp", row=1, col=2)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2
    )

    # Plot 3: Total Entropy vs Time (all radii overlaid)
    fig.add_trace(go.Scatter(x=t, y=total_entropy, name=data_names[idx]), row=1, col=3)
    fig.update_xaxes(title_text="Sweeps", row=1, col=3)
    fig.update_yaxes(title_text="S/Nk", row=1, col=3)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=3
    )

# Display interactive plot
fig.update_layout(title_text=f"Multi-Radius Comparison (Lattice Size: {n})")
fig.show()
