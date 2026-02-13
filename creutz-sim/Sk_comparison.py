import pandas as pd
import glob
import os
import numpy as np
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from pathlib import Path
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate simulation comparison plots')
parser.add_argument('--data-dir', type=str, default='data',
                    help='Directory containing simulation data (default: data)')
args = parser.parse_args()

pio.templates.default = "plotly_white"

# Detect theme from template
is_dark_mode = pio.templates.default == "plotly_dark"

# Max radius
r = 11
# window size for rolling average
bin_size = 10

fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2, vertical_spacing=0.02, row_heights=[0.8, 0.2])
fig2 = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

colors = ['#5A3A5E',
          '#702963',
          '#800020',
          '#AA336A',
          '#9F2B68',
          '#800080',
          '#BF40BF',
          '#DA70D6',
          '#CF9FFF',
          '#E0B0FF',
          '#CBC3E3']

irr_colors = ['#005A5E',
          '#047A80',
          '#068D94',
          '#07A8B0',
          '#08BFC9',
          '#09D3DE',
          '#07E0ED',
          '#A4D3E0',
          '#03F1FF',
          '#A1F3FF',
          '#BCEEF7']

# Use relative path from repo root
repo_root = Path(__file__).parent.parent
filepath = repo_root / args.data_dir

# Check if data directory exists
if not filepath.exists():
    print(f"Error: Data directory not found at {filepath}")
    print("Please run 'make sim' first to generate simulation data.")
    exit(1)

avg_Sk = np.array([])
irr_avg_Sk = np.array([])
SEM = np.array([])
irr_SEM = np.array([])

# Detect if this is a combined run (has both rev and irr data)
is_combined = (filepath / 'irr').exists() and (filepath / 'rev').exists()

# Track which radii already have legend entries
rev_legend_shown = set()
irr_legend_shown = set()

# Track if we have any data to plot
average_df = None
start_index = 0
end_index = 0
n = 0

######### Plot irreversible data (if available)
for R in range(r):
    folder_path = filepath / 'irr' / f'r{R}'
    if not folder_path.exists():
        continue
    
    all_csv_files = glob.glob(str(folder_path / '*.csv'))
    # Filter out _0.csv files as they have too much noise
    all_csv_files = [f for f in all_csv_files if not f.endswith('_0.csv')]
    if not all_csv_files:
        continue

    # Create an empty list to store individual DataFrames
    list_of_dfs = []

    # Loop through each CSV file, read it into a DataFrame, and append to the list
    for file_path in all_csv_files:
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)

    # Concatenate all DataFrames in the list into a single DataFrame
    # ignore_index=True resets the index of the combined DataFrame
    combined_df = pd.concat(list_of_dfs)
    average_df = combined_df.groupby(combined_df.index).mean()

    n = int(average_df['n'][0])
    show_legend = R not in irr_legend_shown
    if show_legend:
        irr_legend_shown.add(R)
    
    # For single runs, show "Radius X"; for combined runs, show empty string (group title shows radius)
    if is_combined:
        trace_name = ""
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=trace_name, 
                                 line=dict(color=irr_colors[R]), legendgroup=f"r{R}", legendgrouptitle_text=f"Radius {R}", 
                                 legendrank=R*2+1, showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=irr_colors[R]), legendgroup=f"r{R}", showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    if is_combined:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=trace_name, line=dict(color=irr_colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=irr_colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df['t']) // 4

    # Calculate the middle index
    middle_index = len(average_df['t']) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    if is_combined:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name="", 
                                 line=dict(color=irr_colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name="", line=dict(color=irr_colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    else:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=irr_colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=irr_colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    irr_avg_Sk = np.append(irr_avg_Sk, np.mean(zoom['S/nk']))
    irr_SEM = np.append(irr_SEM, np.std(zoom['S/nk']/math.sqrt(len(zoom['S/nk']))))

######### Plot reversible data
# Check if data is in rev subdirectory (when run with -f 0)
rev_filepath = filepath / 'rev' if (filepath / 'rev').exists() and any((filepath / 'rev').iterdir()) else filepath

for R in range(r):
    folder_path = rev_filepath / f'r{R}'
    all_csv_files = glob.glob(str(folder_path / '*.csv'))
    # Filter out _0.csv files as they have too much noise
    all_csv_files = [f for f in all_csv_files if not f.endswith('_0.csv')]
    
    if not all_csv_files:
        continue

    # Create an empty list to store individual DataFrames
    list_of_dfs = []

    # Loop through each CSV file, read it into a DataFrame, and append to the list
    for file_path in all_csv_files:
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)

    # Concatenate all DataFrames in the list into a single DataFrame
    # ignore_index=True resets the index of the combined DataFrame
    combined_df = pd.concat(list_of_dfs)
    average_df = combined_df.groupby(combined_df.index).mean()

    n = int(average_df['n'][0])
    show_legend = R not in rev_legend_shown
    if show_legend:
        rev_legend_shown.add(R)
    
    # For single runs, show "Radius X"; for combined runs, show empty string (group title shows radius)
    if is_combined:
        trace_name = ""
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=trace_name, 
                                 line=dict(color=colors[R]), legendgroup=f"r{R}", legendgrouptitle_text=f"Radius {R}", 
                                 legendrank=R*2, showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=colors[R]), legendgroup=f"r{R}", showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    if is_combined:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=trace_name, line=dict(color=colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df['t']) // 4

    # Calculate the middle index
    middle_index = len(average_df['t']) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = middle_index + (num_elements // 2)

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    if is_combined:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name="", 
                                 line=dict(color=colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name="", line=dict(color=colors[R]), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    else:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=colors[R]), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    avg_Sk = np.append(avg_Sk, np.mean(zoom['S/nk']))
    SEM = np.append(SEM, np.std(zoom['S/nk']/math.sqrt(len(zoom['S/nk']))))
    # fig.add_trace(go.Histogram(x=average_df['S/nk'], nbinsx=1),row=1, col=2)

fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=2)

# Only add vlines and vrects if we have data
if average_df is not None:
    fig.add_vline(x=len(average_df['t'])//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)
    fig.add_vline(x=len(average_df['t'])//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)
    fig.add_vrect(x0=start_index, x1=end_index, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=1)
    fig.add_vrect(x0=start_index, x1=end_index, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=2)

# Add toggle buttons to show/hide reversible and irreversible traces
num_traces = len(fig.data)

# Determine which traces are rev vs irr based on line color
# Rev traces use colors list (purple shades), Irr traces use irr_colors list (cyan shades)
rev_visibility = []
irr_visibility = []
for trace in fig.data:
    # Check if trace has a line color property
    if hasattr(trace, 'line') and hasattr(trace.line, 'color'):
        trace_color = trace.line.color
        # Check if color is in the reversible colors (purple shades)
        if trace_color in colors:
            rev_visibility.append(True)
            irr_visibility.append('legendonly')
        # Check if color is in the irreversible colors (cyan shades)
        elif trace_color in irr_colors:
            rev_visibility.append('legendonly')
            irr_visibility.append(True)
        else:
            # Unknown color, show in both
            rev_visibility.append(True)
            irr_visibility.append(True)
    else:
        # No color info, show in both
        rev_visibility.append(True)
        irr_visibility.append(True)

fig.update_layout(
    title_text=f"Lattice Size: {n}",
    height=500,
    legend=dict(
        orientation="v",
        yanchor="top",
        y=1.0,
        xanchor="left",
        x=1.01,
        bgcolor="rgba(30,30,30,0.9)" if is_dark_mode else "rgba(255,255,255,0.9)",
        bordercolor="rgba(100,100,100,0.5)" if is_dark_mode else "rgba(0,0,0,0.3)",
        borderwidth=1,
        font=dict(size=9, color="rgba(255,255,255,0.9)" if is_dark_mode else "rgba(0,0,0,0.9)"),
        grouptitlefont=dict(size=10, family="Arial Black", color="rgba(255,255,255,1)" if is_dark_mode else "rgba(0,0,0,1)"),
        itemsizing='constant',
        itemwidth=30,
        tracegroupgap=0,
        groupclick="togglegroup"
    ),
    updatemenus=[
        dict(
            type="buttons",
            direction="right",
            buttons=[
                dict(
                    args=[{"visible": True}],
                    label="Show All",
                    method="restyle"
                )
            ] + ([
                dict(
                    args=[{"visible": rev_visibility}],
                    label="Show Rev",
                    method="restyle"
                ),
                dict(
                    args=[{"visible": irr_visibility}],
                    label="Show Irr",
                    method="restyle"
                )
            ] if is_combined else []) + [
                dict(
                    args=[{"visible": "legendonly"}],
                    label="Hide All",
                    method="restyle"
                )
            ],
            pad={"r": 0, "t": 0, "b": 0, "l": 0},
            showactive=False,
            x=1,
            xanchor="right",
            y=1.12,
            yanchor="top",
            bgcolor="rgba(255,255,255,0.1)" if is_dark_mode else "rgba(0,0,0,0.05)",
            bordercolor="rgba(100,100,100,0.3)" if is_dark_mode else "rgba(0,0,0,0.2)",
            borderwidth=1,
            font=dict(size=11, color="#1f77b4")
        )
    ],
    annotations=[
        dict(
            text="💡 Double-click legend items to isolate traces",
            xref="paper",
            yref="paper",
            x=1,
            y=-0.05,
            xanchor="right",
            yanchor="top",
            showarrow=False,
            font=dict(size=10, color="rgba(100,100,100,0.7)" if is_dark_mode else "rgba(150,150,150,0.8)"),
            bgcolor="rgba(255,255,255,0.05)" if is_dark_mode else "rgba(255,255,255,0.7)",
            borderpad=4
        )
    ]
)
fig.show()

####################################################################################
###################### Avg Entropy v Radius & Entropy difference ###################
####################################################################################

### Avg Entropy v Radius
fig2.add_trace(go.Scatter(x=np.arange(r), y=avg_Sk, error_y=dict(type='data', array=SEM), name="Reversible", line=dict(color='purple'),
                          hovertemplate='<b>Reversible</b><br>Radius: %{x}<br>Avg S/Nk: %{y:.4f}<extra></extra>'),row=1, col=1)
fig2.add_trace(go.Scatter(x=np.arange(r), y=irr_avg_Sk, error_y=dict(type='data', array=irr_SEM), name="Irreversible", line=dict(color='#00C4C4'),
                          hovertemplate='<b>Irreversible</b><br>Radius: %{x}<br>Avg S/Nk: %{y:.4f}<extra></extra>'),row=1, col=1)

fig2.update_xaxes(title_text="Radius", row=1, col=1)
fig2.update_yaxes(title_text="Avg S/Nk", row=1, col=1)

### Entropy Difference

fig2.update_xaxes(title_text="Sweeps", row=1, col=2)
fig2.update_yaxes(title_text="(S<sub>irr</sub>-S<sub>rev</sub>)/Nk", row=1, col=2)

fig2.update_layout(title_text=f"Lattice Size: {n}")
fig2.show()
