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

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

pio.templates.default = "plotly_white"

# Configuration parameters
r = 11  # Max radius (tests R=0 to R=10)
bin_size = 10  # Window size for rolling average smoothing

# Create 2x2 subplot grid:
# [full data, smoothed] on top row, [zoomed, zoomed smoothed] on bottom
fig = make_subplots(
    rows=2,
    cols=2,
    horizontal_spacing=0.2,
    vertical_spacing=0.02,
    row_heights=[0.8, 0.2],
)

# Color palette: 11 distinct colors for R=0 through R=10 (purple/pink gradient)
colors = [
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


# =============================================================================
# IRREVERSIBLE SIMULATION DATA
# =============================================================================
# Process irreversible simulation results from data/irr/r{R}/ directories
# Combines multiple runs per radius and calculates ensemble averages

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for R in range(r):
    folder_path = os.path.join(project_root, "data", "irr", f"r{R}")
    all_csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

    # Load all CSV files for this radius (multiple independent runs)
    list_of_dfs = []
    for file_path in all_csv_files:
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)

    # Combine all runs and calculate ensemble average for each sweep
    # This gives statistical averaging across M independent runs
    combined_df = pd.concat(list_of_dfs)
    average_df = combined_df.groupby(combined_df.index).mean()

    n = int(average_df["n"][0])
    # Add raw entropy data to left plots
    fig.add_trace(
        go.Scatter(
            x=average_df["t"],
            y=average_df["S/nk"],
            name=f"irr radius {R}",
            line=dict(color=colors[R]),
        ),
        row=1,
        col=1,
    )
    # Add smoothed entropy data to right plots (rolling average reduces noise)
    fig.add_trace(
        go.Scatter(
            x=average_df["t"].rolling(window=bin_size).mean(),
            y=average_df["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius{R}",
            line=dict(color=colors[R]),
        ),
        row=1,
        col=2,
    )

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df["t"]) // 4

    # Calculate the middle index
    middle_index = len(average_df["t"]) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    fig.add_trace(
        go.Scatter(
            x=zoom["t"], y=zoom["S/nk"], name=f"radius {R}", line=dict(color=colors[R])
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=zoom["t"].rolling(window=bin_size).mean(),
            y=zoom["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius {R}",
            line=dict(color=colors[R]),
        ),
        row=2,
        col=1,
    )

# =============================================================================
# REVERSIBLE SIMULATION DATA
# =============================================================================
# Process reversible simulation results from data/r{R}/ directories
# Same processing as irreversible data above

for R in range(r):
    folder_path = os.path.join(project_root, "data", f"r{R}")
    all_csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

    # Load all CSV files for this radius (multiple independent runs)
    list_of_dfs = []
    for file_path in all_csv_files:
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)

    # Combine all runs and calculate ensemble average for each sweep
    combined_df = pd.concat(list_of_dfs)
    average_df = combined_df.groupby(combined_df.index).mean()

    n = int(average_df["n"][0])
    fig.add_trace(
        go.Scatter(
            x=average_df["t"],
            y=average_df["S/nk"],
            name=f"radius {R}",
            line=dict(color=colors[R]),
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=average_df["t"].rolling(window=bin_size).mean(),
            y=average_df["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius{R}",
            line=dict(color=colors[R]),
        ),
        row=1,
        col=2,
    )

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df["t"]) // 4

    # Calculate the middle index
    middle_index = len(average_df["t"]) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    fig.add_trace(
        go.Scatter(
            x=zoom["t"], y=zoom["S/nk"], name=f"radius {R}", line=dict(color=colors[R])
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=zoom["t"].rolling(window=bin_size).mean(),
            y=zoom["S/nk"].rolling(window=bin_size).mean(),
            name=f"radius {R}",
            line=dict(color=colors[R]),
        ),
        row=2,
        col=2,
    )
    # fig.add_trace(go.Histogram(x=average_df['S/nk'], nbinsx=1),row=1, col=2)

# Configure axes labels
fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=2)

# Add red dashed line at midpoint (where forward phase ends, reverse begins)
fig.add_vline(
    x=len(average_df["t"]) // 2,
    line_width=1,
    line_dash="dash",
    line_color="Red",
    row=1,
    col=1,
)
fig.add_vline(
    x=len(average_df["t"]) // 2,
    line_width=1,
    line_dash="dash",
    line_color="Red",
    row=1,
    col=2,
)

# Add blue shaded rectangle showing zoom region on full plots
fig.add_vrect(
    x0=start_index,
    x1=end_index,
    line_width=0,
    fillcolor="blue",
    opacity=0.1,
    row=1,
    col=1,
)
fig.add_vrect(
    x0=start_index,
    x1=end_index,
    line_width=0,
    fillcolor="blue",
    opacity=0.1,
    row=1,
    col=2,
)

# Set title and display
fig.update_layout(
    title_text=f"Entropy Comparison: Reversible vs Irreversible (Lattice Size: {n})"
)
fig.show()
