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
from typing import List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.special import loggamma as logg

pio.templates.default = "plotly_white"


def find_radius_csv_files(
    project_root: str, max_radius: int = 11
) -> Tuple[List[str], List[str]]:
    """
    Find CSV files for each radius from R=0 to R=max_radius-1.

    Args:
        project_root: Root directory of the project
        max_radius: Maximum radius value (exclusive), default 11 searches R=0 to R=10

    Returns:
        Tuple of (file_paths, display_names)
    """
    file_names = []
    data_names = []

    for R in range(max_radius):
        folder_path = os.path.join(project_root, "data", f"r{R}")
        csv_pattern = os.path.join(folder_path, "sim_data*.csv")
        csv_files = glob.glob(csv_pattern)

        if csv_files:
            # Use most recent CSV if multiple found
            csv_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            file_names.append(csv_files[0])
            data_names.append(f"radius = {R}")
        else:
            print(f"Warning: No CSV files found in {folder_path}")

    return file_names, data_names


def calculate_entropy_with_n0_handling(df: pd.DataFrame) -> pd.Series:
    """
    Calculate total entropy per site handling N0=0 edge case.

    Args:
        df: DataFrame with columns n, K, N0, Nx

    Returns:
        Series with total entropy per site (S/Nk)
    """
    n = df["n"][0]
    K = df["K"]
    N0 = df["N0"]
    Nx = df["Nx"]

    # Handle N0=0 edge case: use 2^(N0+1) instead of 2^N0
    N0_exp = df["N0"].replace(0, 1)

    # Kinetic entropy
    Sk = logg(K + n) - logg(K + 1) - logg(n)

    # Configurational entropy
    Su = (
        logg(n + 1)
        + np.log(2**N0_exp)
        - (logg(n - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))
    )

    # Total entropy per site
    return (Sk + Su) / n


def add_radius_trace_to_subplot(
    fig: go.Figure,
    df: pd.DataFrame,
    data_name: str,
    row: int,
    col: int,
    y_column: str,
    y_label: str,
) -> None:
    """
    Add a single trace to a subplot with standard formatting.

    Args:
        fig: Plotly figure object
        df: DataFrame with data
        data_name: Name for legend
        row: Subplot row
        col: Subplot column
        y_column: Column name for y-axis data
        y_label: Label for y-axis
    """
    t = df["t"]
    y = df[y_column]

    fig.add_trace(go.Scatter(x=t, y=y, name=data_name), row=row, col=col)
    fig.update_xaxes(title_text="Sweeps", row=row, col=col)
    fig.update_yaxes(title_text=y_label, row=row, col=col)
    fig.add_vline(
        x=len(df) // 2,
        line_width=1,
        line_dash="dash",
        line_color="Red",
        row=row,
        col=col,
    )


def create_multi_radius_plot(file_names: List[str], data_names: List[str]) -> go.Figure:
    """
    Create 1x3 subplot figure comparing multiple radii.

    Args:
        file_names: List of CSV file paths
        data_names: List of display names for legend

    Returns:
        Plotly Figure object
    """
    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

    lattice_size = None

    for idx, file in enumerate(file_names):
        df = pd.read_csv(file)

        if lattice_size is None:
            lattice_size = df["n"][0]

        # Calculate entropy for plot 3
        df["entropy"] = calculate_entropy_with_n0_handling(df)

        # Add traces for all three subplots
        add_radius_trace_to_subplot(fig, df, data_names[idx], 1, 1, "K", "Demon Energy")
        add_radius_trace_to_subplot(fig, df, data_names[idx], 1, 2, "U", "Lattice Temp")
        add_radius_trace_to_subplot(fig, df, data_names[idx], 1, 3, "entropy", "S/Nk")

    fig.update_layout(
        title_text=f"Multi-Radius Comparison (Lattice Size: {lattice_size})"
    )
    return fig


def main():
    """Main entry point for multi-radius plotting script."""
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Find CSV files
    file_names, data_names = find_radius_csv_files(project_root)

    if not file_names:
        print("Error: No simulation CSV files found in data/r{0-10}/ directories")
        print("Run a simulation first to generate data:")
        print("  make run-sim")
        sys.exit(1)

    print(f"Found {len(file_names)} CSV files to plot")

    # Create and display plot
    fig = create_multi_radius_plot(file_names, data_names)
    fig.show()


if __name__ == "__main__":
    main()
