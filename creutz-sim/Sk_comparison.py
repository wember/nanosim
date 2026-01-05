"""
Copyright (c) 2026 Winry Ember
Licensed under the MIT License.
See LICENSE file in the project root for full license information.

Entropy Comparison Visualization

Compares entropy (S/Nk) evolution between reversible and irreversible
simulations across different demon-coupling radii (R=0 to R=10).

Generates a 2x2 subplot layout:
- Row 1, Col 1: Full time series (raw data)
- Row 1, Col 2: Smoothed full time series (rolling average)
- Row 2, Col 1: Zoomed view of middle 25% (around forward/reverse transition)
- Row 2, Col 2: Smoothed zoomed view

Requires data in:
- data/r{0-10}/ (reversible simulations)
- data/irr/r{0-10}/ (irreversible simulations)
"""

import glob
import os
from typing import List, Tuple

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

pio.templates.default = "plotly_white"


def get_color_palette() -> List[str]:
    """
    Get standard color palette for R=0 through R=10.

    Returns:
        List of 11 hex color strings (purple/pink gradient)
    """
    return [
        "#301934",
        "#702963",
        "#800020",
        "#AA336A",
        "#9F2B68",
        "#800080",
        "#BF40BF",
        "#DA70D6",
        "#CF9FFF",
        "#E0B0FF",
        "#CBC3E3",
    ]


def load_and_average_csv_files(folder_path: str) -> pd.DataFrame:
    """
    Load all CSV files from folder and calculate ensemble average.

    Args:
        folder_path: Path to directory containing CSV files

    Returns:
        DataFrame with averaged values across all CSV files
    """
    all_csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

    if not all_csv_files:
        return pd.DataFrame()

    list_of_dfs = []
    for file_path in all_csv_files:
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)

    combined_df = pd.concat(list_of_dfs)
    average_df = combined_df.groupby(combined_df.index).mean()

    return average_df


def calculate_zoom_indices(total_length: int) -> Tuple[int, int]:
    """
    Calculate start and end indices for middle 25% zoom region.

    Args:
        total_length: Total length of time series

    Returns:
        Tuple of (start_index, end_index)
    """
    num_elements = total_length // 4
    middle_index = total_length // 2
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    return start_index, end_index


def add_entropy_traces(
    fig: go.Figure,
    df: pd.DataFrame,
    radius: int,
    color: str,
    label_prefix: str,
    bin_size: int = 10,
) -> Tuple[int, int]:
    """
    Add entropy traces for one radius to all four subplots.

    Args:
        fig: Plotly figure object
        df: DataFrame with 't' and 'S/nk' columns
        radius: Radius value
        color: Hex color string
        label_prefix: Prefix for legend labels ('irr' or '')
        bin_size: Window size for rolling average

    Returns:
        Tuple of (start_index, end_index) for zoom region
    """
    # Full time series (raw)
    fig.add_trace(
        go.Scatter(
            x=df["t"],
            y=df["S/nk"],
            name=f"{label_prefix}radius {radius}",
            line=dict(color=color),
        ),
        row=1,
        col=1,
    )

    # Full time series (smoothed)
    fig.add_trace(
        go.Scatter(
            x=df["t"].rolling(window=bin_size).mean(),
            y=df["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius{radius}",
            line=dict(color=color),
        ),
        row=1,
        col=2,
    )

    # Zoom region
    start_index, end_index = calculate_zoom_indices(len(df["t"]))
    zoom = df.iloc[start_index:end_index]

    # Zoom (raw)
    fig.add_trace(
        go.Scatter(
            x=zoom["t"],
            y=zoom["S/nk"],
            name=f"radius {radius}",
            line=dict(color=color),
        ),
        row=2,
        col=1,
    )

    # Zoom (smoothed)
    fig.add_trace(
        go.Scatter(
            x=zoom["t"].rolling(window=bin_size).mean(),
            y=zoom["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius {radius}",
            line=dict(color=color),
        ),
        row=2,
        col=2,
    )

    return start_index, end_index


def create_entropy_comparison_plot(
    project_root: str, max_radius: int = 11, bin_size: int = 10
) -> go.Figure:
    """
    Create 2x2 subplot comparing reversible and irreversible entropy evolution.

    Args:
        project_root: Root directory of project
        max_radius: Maximum radius value (exclusive)
        bin_size: Window size for rolling average smoothing

    Returns:
        Plotly Figure object
    """
    # Create 2x2 subplot grid
    fig = make_subplots(
        rows=2,
        cols=2,
        horizontal_spacing=0.2,
        vertical_spacing=0.02,
        row_heights=[0.8, 0.2],
    )

    colors = get_color_palette()
    lattice_size = None
    zoom_start, zoom_end = None, None
    midpoint = None

    # Process irreversible simulations
    for R in range(max_radius):
        folder_path = os.path.join(project_root, "data", "irr", f"r{R}")
        average_df = load_and_average_csv_files(folder_path)

        if average_df.empty:
            continue

        if lattice_size is None:
            lattice_size = int(average_df["n"][0])
            midpoint = len(average_df["t"]) // 2

        zoom_start, zoom_end = add_entropy_traces(
            fig, average_df, R, colors[R], "irr ", bin_size
        )

    # Process reversible simulations
    for R in range(max_radius):
        folder_path = os.path.join(project_root, "data", f"r{R}")
        average_df = load_and_average_csv_files(folder_path)

        if average_df.empty:
            continue

        if lattice_size is None:
            lattice_size = int(average_df["n"][0])
            midpoint = len(average_df["t"]) // 2

        zoom_start, zoom_end = add_entropy_traces(
            fig, average_df, R, colors[R], "", bin_size
        )

    # Configure axes
    fig.update_xaxes(title_text="Sweeps", row=1, col=1)
    fig.update_yaxes(title_text="S/Nk", row=1, col=1)
    fig.update_yaxes(title_text="S/Nk", row=1, col=2)

    # Add midpoint lines
    if midpoint is not None:
        fig.add_vline(
            x=midpoint,
            line_width=1,
            line_dash="dash",
            line_color="Red",
            row=1,
            col=1,
        )
        fig.add_vline(
            x=midpoint,
            line_width=1,
            line_dash="dash",
            line_color="Red",
            row=1,
            col=2,
        )

    # Add zoom region rectangles
    if zoom_start is not None and zoom_end is not None:
        fig.add_vrect(
            x0=zoom_start,
            x1=zoom_end,
            line_width=0,
            fillcolor="blue",
            opacity=0.1,
            row=1,
            col=1,
        )
        fig.add_vrect(
            x0=zoom_start,
            x1=zoom_end,
            line_width=0,
            fillcolor="blue",
            opacity=0.1,
            row=1,
            col=2,
        )

    # Set title
    title = f"Entropy Comparison: Reversible vs Irreversible"
    if lattice_size is not None:
        title += f" (Lattice Size: {lattice_size})"
    fig.update_layout(title_text=title)

    return fig


def main():
    """Main entry point for entropy comparison visualization."""
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    fig = create_entropy_comparison_plot(project_root)
    fig.show()


if __name__ == "__main__":
    main()
