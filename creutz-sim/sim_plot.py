"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Single Simulation Visualization

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

import glob
import os
import sys
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.special import loggamma as logg

pio.templates.default = "plotly_white"


# =============================================================================
# Testable Functions
# =============================================================================


def find_most_recent_csv(
    project_root: str, search_patterns: Optional[List[str]] = None
) -> Optional[str]:
    """Find the most recent CSV file matching search patterns.

    Args:
        project_root: Root directory to search from
        search_patterns: List of glob patterns to search. If None, uses default
                         patterns.

    Returns:
        Path to most recent CSV file, or None if no files found
    """
    if search_patterns is None:
        search_patterns = [
            os.path.join(project_root, "data", "r*", "sim_data*.csv"),
            os.path.join(project_root, "data", "irr", "r*", "sim_data*.csv"),
        ]

    csv_files = []
    for pattern in search_patterns:
        csv_files.extend(glob.glob(pattern))

    if not csv_files:
        return None

    # Sort by modification time, most recent first
    csv_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return csv_files[0]


def validate_csv_columns(
    df: pd.DataFrame, required_cols: Optional[List[str]] = None
) -> Tuple[bool, List[str]]:
    """Validate that DataFrame has required columns.

    Args:
        df: DataFrame to validate
        required_cols: List of required column names. If None, uses default list.

    Returns:
        Tuple of (is_valid, missing_columns)
    """
    if required_cols is None:
        required_cols = ["t", "K", "U", "N0", "Nx", "n"]

    missing_cols = [col for col in required_cols if col not in df.columns]
    return len(missing_cols) == 0, missing_cols


def load_and_validate_csv(csv_file: str) -> pd.DataFrame:
    """Load CSV file and validate it has required columns.

    Args:
        csv_file: Path to CSV file

    Returns:
        Validated DataFrame

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns are missing
    """
    df = pd.read_csv(csv_file)

    is_valid, missing_cols = validate_csv_columns(df)
    if not is_valid:
        raise ValueError(
            f"CSV file missing required columns: {missing_cols}. "
            f"Found columns: {list(df.columns)}"
        )

    return df


def calculate_entropy(df: pd.DataFrame) -> pd.Series:
    """Calculate total entropy per site from simulation data.

    Args:
        df: DataFrame with columns t, K, U, N0, Nx, n

    Returns:
        Series containing entropy values (S/Nk)
    """
    n = df["n"].iloc[0]
    K = df["K"]
    N0 = df["N0"]
    Nx = df["Nx"]

    # Handle N0=0 edge case
    N0_exp = df["N0"].replace(0, 1)

    # Sk: Kinetic entropy from demon energy distribution
    Sk = logg(K + n) - logg(K + 1) - logg(n)

    # Su: Configurational entropy from spin/bond states
    Su = (
        logg(n + 1)
        + N0_exp * np.log(2)
        - (logg(n - N0 - Nx + 1) + logg(N0 + 1) + logg(Nx + 1))
    )

    # Total entropy per site
    return (Sk + Su) / n


def create_simulation_plot(df: pd.DataFrame, title: Optional[str] = None) -> go.Figure:
    """Create 1x3 subplot visualization of simulation results.

    Args:
        df: DataFrame with simulation data
        title: Optional plot title. If None, generates from lattice size.

    Returns:
        Plotly Figure object
    """
    n = df["n"].iloc[0]
    t = df["t"]
    U = df["U"]
    K = df["K"]

    # Calculate entropy
    total_entropy = calculate_entropy(df)

    # Create subplot layout
    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

    # Plot 1: Demon Energy
    fig.add_trace(go.Scatter(x=t, y=K, showlegend=False), row=1, col=1)
    fig.update_xaxes(title_text="Sweeps", row=1, col=1)
    fig.update_yaxes(title_text="Demon Energy", row=1, col=1)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1
    )

    # Plot 2: Lattice Temperature (Energy)
    fig.add_trace(go.Scatter(x=t, y=U, showlegend=False), row=1, col=2)
    fig.update_xaxes(title_text="Sweeps", row=1, col=2)
    fig.update_yaxes(title_text="Lattice Temp", row=1, col=2)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2
    )

    # Plot 3: Total Entropy
    fig.add_trace(go.Scatter(x=t, y=total_entropy, showlegend=False), row=1, col=3)
    fig.update_xaxes(title_text="Sweeps", row=1, col=3)
    fig.update_yaxes(title_text="S/Nk", row=1, col=3)
    fig.add_vline(
        x=len(df) // 2, line_width=1, line_dash="dash", line_color="Red", row=1, col=3
    )

    # Set title
    if title is None:
        title = f"Single Simulation Results (Lattice Size: {n})"
    fig.update_layout(title_text=title)

    return fig


# =============================================================================
# Main Script Execution
# =============================================================================


def main():
    """Main function for command-line execution."""
    # Get CSV file from command line or auto-discover
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        project_root = os.path.dirname(os.path.dirname(__file__))
        csv_file = find_most_recent_csv(project_root)

        if csv_file is None:
            print("Error: No simulation CSV files found")
            print("Run a simulation first to generate data:")
            print("  make run-sim-small")
            print("")
            print("Or specify a CSV file:")
            print("  python creutz-sim/sim_plot.py path/to/sim_data.csv")
            sys.exit(1)

        print(f"Using most recent simulation CSV: {csv_file}")

    # Load and validate CSV
    try:
        df = load_and_validate_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        print("")
        print("This script requires simulation output CSV, not analysis pipeline CSV.")
        print("Run: make run-sim-small")
        sys.exit(1)

    # Create and display plot
    fig = create_simulation_plot(df)
    fig.show()


if __name__ == "__main__":
    main()
