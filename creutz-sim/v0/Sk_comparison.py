import pandas as pd
import glob
import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"

# Max radius
r = 11
# window size for rolling average
bin_size = 10

fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2, vertical_spacing=0.02, row_heights=[0.8, 0.2])

colors = ['#301934',
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


######### Plot irreversible sim
for R in range(r):
    folder_path = f'/Users/winry/Documents/ASU/thesis/dev/data/irr/r{R}'
    all_csv_files = glob.glob(os.path.join(folder_path, '*.csv'))

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
    fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"irr radius {R}", line=dict(color=colors[R])),row=1, col=1)
    fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), name=f"radius{R}", line=dict(color=colors[R])),row=1, col=2)

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df['t']) // 4

    # Calculate the middle index
    middle_index = len(average_df['t']) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"radius {R}", line=dict(color=colors[R])),row=2, col=1)
    fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), name=f"radius {R}", line=dict(color=colors[R])),row=2, col=1)
########### Reversible sims with radius 0 to r
for R in range(r):
    folder_path = f'/Users/winry/Documents/ASU/thesis/dev/data/r{R}'
    all_csv_files = glob.glob(os.path.join(folder_path, '*.csv'))

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
    fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"radius {R}", line=dict(color=colors[R])),row=1, col=1)
    fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), name=f"radius{R}", line=dict(color=colors[R])),row=1, col=2)

    # Zoomed in portion about center of dataframe
    num_elements = len(average_df['t']) // 4

    # Calculate the middle index
    middle_index = len(average_df['t']) // 2

    # Calculate the start and end indices for the middle elements
    start_index = middle_index - (num_elements // 2)
    end_index = start_index + num_elements

    # Extract the middle elements
    zoom = average_df.iloc[start_index:end_index]
    fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"radius {R}", line=dict(color=colors[R])),row=2, col=1)
    fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), name=f"radius {R}", line=dict(color=colors[R])),row=2, col=2)
    # fig.add_trace(go.Histogram(x=average_df['S/nk'], nbinsx=1),row=1, col=2)

fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=2)
fig.add_vline(x=len(average_df['t'])//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)
fig.add_vline(x=len(average_df['t'])//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)
fig.add_vrect(x0=start_index, x1=end_index, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=1)
fig.add_vrect(x0=start_index, x1=end_index, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=2)

fig.update_layout(title_text=f"Lattice Size: {n}")
fig.show()
