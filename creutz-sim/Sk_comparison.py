import pandas as pd
import glob
import os
import numpy as np

from scipy.special import loggamma as logg
from scipy.optimize import curve_fit

def Sk(N, K):
    return logg(K + N) - logg(K + 1) - logg(N)

def Su(N, N0, Nx):
    return logg(N+1) + N0 * np.log(2) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

def Su0(N, N0, Nx):
    return logg(N+1) + (N0+1) * np.log(2) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1))

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

# window size for rolling average
bin_size = 10

# Target number of block-averaged symbols per trace in the ln(p_n/p_n+1) /
# kT/J plot (Figure 2D recreation). Each symbol is a genuine average over a
# block of ~(total rows / this number) consecutive recorded rows, positioned
# at the block's median 't' - matching the paper's own methodology of
# averaging over large sweep windows (e.g. 10^4 sweeps) rather than plotting
# every recorded row or applying a small rolling-window smoothing. Also keeps
# the HTML payload for that figure a manageable size on long runs.
RATIO_PLOT_N_POINTS = 15

fig = make_subplots(rows=3, cols=2, horizontal_spacing=0.2, vertical_spacing=0.08,
                     row_heights=[0.5, 0.2, 0.3],
                     specs=[[{}, {}], [{}, {}], [{"colspan": 2}, None]])

def log_bin_average(freqs, psd, n_bins=120):
    """Average PSD values into bins that are evenly spaced in log10(f), so that
    the plotted curve doesn't get visually 'crowded' at high frequencies just
    because linearly-spaced FFT/Welch bins pile up when viewed on a log axis.
    Averaging is done in linear power (not dB) within each bin, which is the
    physically correct way to average power spectral densities."""
    logf = np.log10(freqs)
    edges = np.linspace(logf.min(), logf.max(), n_bins + 1)
    bin_idx = np.clip(np.digitize(logf, edges) - 1, 0, n_bins - 1)

    binned_f, binned_psd = [], []
    for i in range(n_bins):
        sel = bin_idx == i
        if np.any(sel):
            binned_f.append(np.mean(freqs[sel]))
            binned_psd.append(np.mean(psd[sel]))
    return np.array(binned_f), np.array(binned_psd)


def add_psd_trace(fig, list_of_dfs, color, R, is_irreversible, row=3, col=1):
    """Compute the power spectral density of sqrt(K) (the demon/heat-bath energy),
    matching the paper's inset in Figure 1 of Chamberlin (2024) as closely as
    possible:

    - PSD = |FFT(sqrt(K))|^2, taken over the FULL simulation record of each
      individual run (not a short window), since the paper's x-axis range
      (10*log10(f) up to ~50, i.e. f up to ~10^5) only makes sense as the raw
      FFT bin index over their full ~131,072-point records — a short window
      can't reach anywhere near that range.
    - "Relative frequency" is therefore the raw FFT bin index k = 1, 2, 3, ...
      (not k/N), so 10*log10(f) directly means 10*log10(k).
    - list_of_dfs holds the individual, un-averaged CSV runs for this radius
      (the _0.csv file is already excluded upstream as "too noisy"). The PSD
      is computed separately per run, then averaged across runs in linear
      power — this mirrors the paper averaging over repeated simulations.
    - The averaged PSD is then bin-averaged in log-frequency space purely for
      display smoothness (see log_bin_average docstring); this step is a
      plotting aid, not part of the paper's stated method.
    """
    all_psds = []
    for df in list_of_dfs:
        df_sorted = df.sort_values('t')
        K = df_sorted['E_demon'].values
        if len(K) < 10:
            continue
        signal = np.sqrt(np.clip(K, a_min=0, a_max=None))

        fft_vals = np.fft.rfft(signal)
        psd = np.abs(fft_vals) ** 2
        k = np.arange(len(psd))  # bin index, 0 = DC
        # Drop the DC bin: log10(0) is undefined and it's not a "frequency"
        all_psds.append((k[1:], psd[1:]))

    if not all_psds:
        return

    # Individual runs may differ slightly in length; align on the shortest
    # common bin range rather than assuming they all match exactly.
    min_len = min(len(k) for k, _ in all_psds)
    k_common = all_psds[0][0][:min_len].astype(float)
    psd_stack = np.array([psd[:min_len] for _, psd in all_psds])
    psd_avg = psd_stack.mean(axis=0)  # average across repeated runs, in linear power

    # Smooth out the log-axis crowding artifact (see log_bin_average docstring)
    freqs, psd = log_bin_average(k_common, psd_avg, n_bins=150)
    if len(freqs) == 0:
        return

    # Match the paper's axes: 10*log10(f) on x, PSD in dB on y, both on LINEAR
    # plotly axes (the dB conversion already does the "log" work).
    x_db = 10 * np.log10(freqs)
    y_db = 10 * np.log10(psd)

    label = "Irreversible" if is_irreversible else "Reversible"
    fig.add_trace(
        go.Scatter(
            x=x_db, y=y_db, mode='lines',
            line=dict(color=color), legendgroup=f"r{R}", showlegend=False,
            hovertemplate=f'<b>Radius {R} ({label})</b><br>10log(f): %{{x:.2f}}<br>PSD: %{{y:.2f}} dB<extra></extra>'
        ),
        row=row, col=col,
    )

fig2 = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

### Figure recreating Chamberlin (2024) Figure 2D: moving averages of
### ln(p_n / p_(n+1)) vs t (sweeps), for the ratio of occupation probabilities
### between adjacent demon energy levels n=0,1,2,3. Shown side by side with an
### equivalent temperature view: since ln(p_n/p_(n+1)) = J/kT for a Boltzmann
### distribution, kT/J = 1 / ln(p_n/p_(n+1)) is the same information rescaled
### as an effective temperature, which converges to a single value only when
### the dynamics is truly thermal.
fig3 = make_subplots(rows=1, cols=2, horizontal_spacing=0.15,
                      subplot_titles=("ln(p<sub>n</sub>/p<sub>n+1</sub>) (\u221d 1/T)",
                                       "kT/J (\u221d T)"))

# Legend-only entries mapping marker shape -> energy level n. These carry no
# data (x/y are empty) so they don't plot anything on either panel; they
# exist purely so the shape/level mapping shows up as real, always-visible
# legend entries instead of only in the caption below the figure.
_shape_legend_symbols = ['square', 'circle', 'triangle-up', 'triangle-down']
for _n_level, _symbol in enumerate(_shape_legend_symbols):
    fig3.add_trace(
        go.Scatter(
            x=[], y=[],
            mode='markers',
            marker=dict(symbol=_symbol, size=8, color='rgba(80,80,80,0.9)'),
            name=f"n={_n_level}",
            legendgroup="shape_legend",
            legendgrouptitle_text="Marker shape = energy level" if _n_level == 0 else None,
            showlegend=True,
        ),
        row=1, col=1,
    )

def add_ratio_traces(fig, average_df, color, R, is_irreversible):
    """
    Add ln(p_n / p_(n+1)) and kT/J moving-average traces for n = 0..3 to `fig`
    (col 1 and col 2 respectively), reproducing Figure 2D of Chamberlin (2024)
    and its equivalent effective-temperature view.

    Symbol shape identifies the energy level n (square/circle/up-triangle/
    down-triangle for n=0,1,2,3), matching the paper. Unlike the paper (which
    uses a single black/red color pair for reversible/irreversible), color here
    follows the same per-radius palette used elsewhere in this script, since
    this script compares multiple radii at once; dynamics type is instead
    distinguished via the legend group label and hover text.

    Averaging matches the paper's own methodology (see Appendix A / Figure 2
    caption): each plotted symbol is the average of the raw ratio over a
    block of consecutive recorded rows (not a rolling/decimated smoothing),
    positioned at the median 't' of that block. This also keeps the number of
    plotted points bounded (see RATIO_PLOT_N_POINTS) so the HTML stays a
    manageable size for browsers even on very long runs.
    """
    if average_df is None or len(average_df) == 0:
        return

    label = "Irreversible" if is_irreversible else "Reversible"
    ratio_cols = ['p0/p1', 'p1/p2', 'p2/p3', 'p3/p4']
    symbols = ['square', 'circle', 'triangle-up', 'triangle-down']

    n_rows = len(average_df)
    block_size = max(1, n_rows // RATIO_PLOT_N_POINTS)
    block_id = np.arange(n_rows) // block_size

    t_blocks = average_df['t'].groupby(block_id).median()

    for n_level, (col, symbol) in enumerate(zip(ratio_cols, symbols)):
        if col not in average_df:
            continue
        raw = average_df[col].replace(0, np.nan)
        ln_ratio = np.log(raw).groupby(block_id).mean()

        fig.add_trace(
            go.Scatter(
                x=t_blocks, y=ln_ratio,
                mode='lines+markers',
                marker=dict(symbol=symbol, size=6, color=color),
                line=dict(color=color, width=1.5),
                name=f"Radius {R} ({label})",
                legendgroup=f"r{R}_{label}",
                legendgrouptitle_text=f"Radius {R} ({label})" if n_level == 0 else None,
                showlegend=(n_level == 0),
                hovertemplate=(
                    f'<b>Radius {R} ({label}), n={n_level}</b><br>'
                    't: %{x}<br>ln(p<sub>n</sub>/p<sub>n+1</sub>): %{y:.4f}<extra></extra>'
                ),
            ),
            row=1, col=1,
        )

        # kT/J = 1 / ln(p_n/p_(n+1)); replace non-positive or exactly-zero
        # ratios (which give a negative, undefined, or infinite "temperature")
        # with NaN rather than plotting an unphysical value or breaking JSON
        # serialization with an inf.
        kT_over_J = ln_ratio.where(ln_ratio > 0, np.nan).rdiv(1).replace([np.inf, -np.inf], np.nan)
        fig.add_trace(
            go.Scatter(
                x=t_blocks, y=kT_over_J,
                mode='lines+markers',
                marker=dict(symbol=symbol, size=6, color=color),
                line=dict(color=color, width=1.5),
                name=f"Radius {R} ({label})",
                legendgroup=f"r{R}_{label}",
                showlegend=False,
                hovertemplate=(
                    f'<b>Radius {R} ({label}), n={n_level}</b><br>'
                    't: %{x}<br>kT/J: %{y:.4f}<extra></extra>'
                ),
            ),
            row=1, col=2,
        )

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


def get_available_radii(base_path):
    # Dynamic radius discovery: read existing r* folders instead of hard-coding max R.
    if not base_path.exists():
        return []
    radii = []
    for radius_dir in base_path.glob('r*'):
        if radius_dir.is_dir():
            suffix = radius_dir.name[1:]
            if suffix.isdigit():
                radii.append(int(suffix))
    return sorted(radii)


def is_run_directory(path):
    # A valid run can be either combined layout (rev/irr) or flat r* folders.
    if (path / 'rev').exists() or (path / 'irr').exists():
        return True
    return len(get_available_radii(path)) > 0


def resolve_data_directory(path):
    # If caller passes data/ (archive root), auto-select the newest valid run folder.
    # If caller passes a specific run folder, use it directly.
    if is_run_directory(path):
        return path

    candidate_runs = []
    for subdir in path.iterdir():
        if subdir.is_dir() and is_run_directory(subdir):
            candidate_runs.append(subdir)

    if not candidate_runs:
        return path

    latest_run = max(candidate_runs, key=lambda p: p.stat().st_mtime)
    print(f"Using latest run directory: {latest_run}")
    return latest_run


filepath = resolve_data_directory(filepath)

avg_Sk = np.array([])
irr_avg_Sk = np.array([])
SEM = np.array([])
irr_SEM = np.array([])

# Detect if this is a combined run (has both rev and irr data)
is_combined = (filepath / 'irr').exists() and (filepath / 'rev').exists()

irr_basepath = filepath / 'irr'
# Dynamic: plot whatever radii are actually present under irr/r*.
irr_radii = get_available_radii(irr_basepath)

# Track which radii already have legend entries
rev_legend_shown = set()
irr_legend_shown = set()

# Track if we have any data to plot
average_df = None
start_index = 0
end_index = 0
n = 0
rev_plotted_radii = []
irr_plotted_radii = []

######### Plot irreversible data (if available)
for R in irr_radii:
    folder_path = irr_basepath / f'r{R}'
    if not folder_path.exists():
        continue

    # Use modulo for safety when available radii exceed palette length.
    irr_color = irr_colors[R % len(irr_colors)]
    
    all_csv_files = glob.glob(str(folder_path / '*.csv'))
    # Filter out _0.csv files as they have too much noise
    all_csv_files = [f for f in all_csv_files if not f.endswith('_0.csv')]
    if not all_csv_files:
        continue

    print(f"[irr] R={R}: reading {len(all_csv_files)} CSV file(s)...", flush=True)

    # Create an empty list to store individual DataFrames
    list_of_dfs = []

    # Loop through each CSV file, read it into a DataFrame, and append to the list
    for idx, file_path in enumerate(all_csv_files):
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)
        pct = int((idx + 1) / len(all_csv_files) * 100)
        print(f"[irr] R={R}: {idx+1}/{len(all_csv_files)} ({pct}%) — {Path(file_path).name} ({len(df):,} rows)", flush=True)

    # Concatenate all DataFrames in the list into a single DataFrame
    # ignore_index=True resets the index of the combined DataFrame
    combined_df = pd.concat(list_of_dfs, ignore_index=True)
    average_df = combined_df.groupby('t', as_index=False).mean()
    # Recompute total entropy from stored observables
    N = int(average_df['n'][0])
    N0 = average_df['N0 (%)'] * N / 100
    Nx = average_df['Nx (%)'] * N / 100
    K  = average_df['E_demon']

    S_conf = np.where(
        N0 == 0,
        Su0(N, N0, Nx),
        Su(N, N0, Nx)
    )

    S_total = (Sk(N, K) + S_conf) / N
    average_df['S/nk'] = S_total


    n = int(average_df['n'][0])
    show_legend = R not in irr_legend_shown
    if show_legend:
        irr_legend_shown.add(R)
    
    # For single runs, show "Radius X"; for combined runs, show empty string (group title shows radius)
    if is_combined:
        trace_name = ""
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=trace_name, 
                                 line=dict(color=irr_color), legendgroup=f"r{R}", legendgrouptitle_text=f"Radius {R}", 
                                 legendrank=R*2+1, showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=irr_color), legendgroup=f"r{R}", showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    if is_combined:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=trace_name, line=dict(color=irr_color), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=irr_color), legendgroup=f"r{R}", showlegend=False,
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
                                 line=dict(color=irr_color), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name="", line=dict(color=irr_color), legendgroup=f"r{R}", 
                                 legendrank=R*2+1, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    else:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=irr_color), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=irr_color), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Irreversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    add_psd_trace(fig, list_of_dfs, irr_color, R, is_irreversible=True)
    add_ratio_traces(fig3, average_df, irr_color, R, is_irreversible=True)

    irr_avg_Sk = np.append(irr_avg_Sk, np.mean(average_df['S/nk']))
    irr_SEM = np.append(irr_SEM, np.std(average_df['S/nk']/math.sqrt(len(average_df['S/nk']))))
    # Track exact radii that produced points (after file filtering) for fig2 x-axis alignment.
    irr_plotted_radii.append(R)
    print(f"[irr] R={R}: done ({len(average_df):,} rows)", flush=True)

######### Plot reversible data
# Check if data is in rev subdirectory (when run with -f 0)
rev_filepath = filepath / 'rev' if (filepath / 'rev').exists() and any((filepath / 'rev').iterdir()) else filepath
# Dynamic: supports both rev/r* subfolders and flat r* folder layout.
rev_radii = get_available_radii(rev_filepath)

if not irr_radii and not rev_radii:
    print(f"Error: No simulation radius folders found under {filepath}")
    print("If your data is in a specific run folder, use: --data-dir data/<timestamp>")
    exit(1)

for R in rev_radii:
    folder_path = rev_filepath / f'r{R}'
    # Use modulo for safety when available radii exceed palette length.
    rev_color = colors[R % len(colors)]
    all_csv_files = glob.glob(str(folder_path / '*.csv'))
    # Filter out _0.csv files as they have too much noise
    all_csv_files = [f for f in all_csv_files if not f.endswith('_0.csv')]
    
    if not all_csv_files:
        continue

    print(f"[rev] R={R}: reading {len(all_csv_files)} CSV file(s)...", flush=True)

    # Create an empty list to store individual DataFrames
    list_of_dfs = []

    # Loop through each CSV file, read it into a DataFrame, and append to the list
    for idx, file_path in enumerate(all_csv_files):
        df = pd.read_csv(file_path)
        list_of_dfs.append(df)
        pct = int((idx + 1) / len(all_csv_files) * 100)
        print(f"[rev] R={R}: {idx+1}/{len(all_csv_files)} ({pct}%) — {Path(file_path).name} ({len(df):,} rows)", flush=True)

    # Concatenate all DataFrames in the list into a single DataFrame
    # ignore_index=True resets the index of the combined DataFrame
    combined_df = pd.concat(list_of_dfs, ignore_index=True)
    average_df = combined_df.groupby('t', as_index=False).mean()
    # Recompute total entropy from stored observables
    N = int(average_df['n'][0])
    N0 = average_df['N0 (%)'] * N / 100
    Nx = average_df['Nx (%)'] * N / 100
    K  = average_df['E_demon']

    S_conf = np.where(
        N0 == 0,
        Su0(N, N0, Nx),
        Su(N, N0, Nx)
    )

    S_total = (Sk(N, K) + S_conf) / N
    average_df['S/nk'] = S_total


    n = int(average_df['n'][0])
    show_legend = R not in rev_legend_shown
    if show_legend:
        rev_legend_shown.add(R)
    
    # For single runs, show "Radius X"; for combined runs, show empty string (group title shows radius)
    if is_combined:
        trace_name = ""
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=trace_name, 
                                 line=dict(color=rev_color), legendgroup=f"r{R}", legendgrouptitle_text=f"Radius {R}", 
                                 legendrank=R*2, showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'], y=average_df['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=rev_color), legendgroup=f"r{R}", showlegend=show_legend,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=1)
    if is_combined:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=trace_name, line=dict(color=rev_color), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=1, col=2)
    else:
        fig.add_trace(go.Scatter(x=average_df['t'].rolling(window=bin_size).mean(), y=average_df['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=rev_color), legendgroup=f"r{R}", showlegend=False,
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
                                 line=dict(color=rev_color), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name="", line=dict(color=rev_color), legendgroup=f"r{R}", 
                                 legendrank=R*2, showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    else:
        fig.add_trace(go.Scatter(x=zoom['t'], y=zoom['S/nk'], name=f"Radius {R}", 
                                 line=dict(color=rev_color), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=1)
        fig.add_trace(go.Scatter(x=zoom['t'].rolling(window=bin_size).mean(), y=zoom['S/nk'].rolling(window=bin_size).mean(), 
                                 name=f"Radius {R}", line=dict(color=rev_color), legendgroup=f"r{R}", showlegend=False,
                                 hovertemplate=f'<b>Radius {R} (Reversible)</b><br>Sweep: %{{x:.1f}}<br>S/Nk: %{{y:.4f}}<extra></extra>'),row=2, col=2)
    add_psd_trace(fig, list_of_dfs, rev_color, R, is_irreversible=False)
    add_ratio_traces(fig3, average_df, rev_color, R, is_irreversible=False)

    avg_Sk = np.append(avg_Sk, np.mean(average_df['S/nk']))
    SEM = np.append(SEM, np.std(average_df['S/nk']/math.sqrt(len(average_df['S/nk']))))
    # Track exact radii that produced points (after file filtering) for fig2 x-axis alignment.
    rev_plotted_radii.append(R)
    print(f"[rev] R={R}: done ({len(average_df):,} rows)", flush=True)
    # fig.add_trace(go.Histogram(x=average_df['S/nk'], nbinsx=1),row=1, col=2)

fig.update_xaxes(title_text="Sweeps", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=1)
fig.update_yaxes(title_text="S/Nk", row=1, col=2)
fig.update_xaxes(title_text="10 log(f)", row=3, col=1)
fig.update_yaxes(title_text="PSD (dB)", row=3, col=1)
print("Serializing fig1 to HTML...", flush=True)

# Only add vlines and vrects if we have data
if average_df is not None:
    t_values = average_df['t'].values
    mid_t = t_values[len(t_values) // 2]
    t_start = t_values[start_index]
    t_end = t_values[end_index - 1]

    fig.add_vline(x=mid_t, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)
    fig.add_vline(x=mid_t, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)
    fig.add_vrect(x0=t_start, x1=t_end, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=1)
    fig.add_vrect(x0=t_start, x1=t_end, line_width=0, fillcolor="blue", opacity=0.1, row=1, col=2)

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
    height=850,
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
fig2.add_trace(go.Scatter(x=rev_plotted_radii, y=avg_Sk, error_y=dict(type='data', array=SEM), name="Reversible", line=dict(color='purple'),
                          hovertemplate='<b>Reversible</b><br>Radius: %{x}<br>Avg S/Nk: %{y:.4f}<extra></extra>'),row=1, col=1)
fig2.add_trace(go.Scatter(x=irr_plotted_radii, y=irr_avg_Sk, error_y=dict(type='data', array=irr_SEM), name="Irreversible", line=dict(color='#00C4C4'),
                          hovertemplate='<b>Irreversible</b><br>Radius: %{x}<br>Avg S/Nk: %{y:.4f}<extra></extra>'),row=1, col=1)

fig2.update_xaxes(title_text="Local bath size", range=[0, 100], autorange=False, row=1, col=1)
fig2.update_yaxes(title_text="Avg S/Nk", row=1, col=1)

### Power-law fit of irreversible avg S/Nk vs radius (y = A x^B), excluding r=0
def power_law(x, A, B, C):
    return A * np.power(x, B) + C

def fit_limit_str(A, B, C, C_err):
    """
    Describe the x -> infinity limit of y = A*x^B + C, for display in a legend.
    If B < 0, the power-law term decays to 0 and the curve approaches C.
    If B >= 0, the power-law term does not vanish (or grows), so there is no
    finite limit; A's sign determines whether it diverges to +inf or -inf.
    """
    if B < 0:
        return f"lim={C:.5f}\u00b1{C_err:.5f}"
    elif B == 0:
        return "lim=undefined (B=0)"
    else:
        return "diverges (B>0)"

irr_r_arr_fit = np.array(irr_plotted_radii)
irr_avg_arr_fit = np.array(irr_avg_Sk)
fit_mask = irr_r_arr_fit > 0

if np.sum(fit_mask) >= 3:
    try:
        popt, pcov = curve_fit(
            power_law,
            irr_r_arr_fit[fit_mask],
            irr_avg_arr_fit[fit_mask],
            p0=[1.0, -1.0, float(np.mean(irr_avg_arr_fit[fit_mask]))],
            maxfev=10000,
        )
        A_fit, B_fit, C_fit = popt
        perr = np.sqrt(np.diag(pcov))

        r_fit_x = np.linspace(irr_r_arr_fit[fit_mask].min(), max(irr_r_arr_fit[fit_mask].max(), 100), 200)
        r_fit_y = power_law(r_fit_x, A_fit, B_fit, C_fit)

        fig2.add_trace(
            go.Scatter(
                x=r_fit_x,
                y=r_fit_y,
                mode='lines',
                name=f"Irr fit: y={A_fit:.4g}x^{B_fit:.4g}+{C_fit:.4g} ({fit_limit_str(A_fit, B_fit, C_fit, perr[2])})",
                line=dict(color='#00C4C4', dash='dash'),
                hovertemplate='<b>Irreversible Fit</b><br>Radius: %{x:.2f}<br>S/Nk: %{y:.6f}<extra></extra>',
            ),
            row=1, col=1,
        )

        print(
            f"Irreversible power-law fit (r>0): A={A_fit:.6g} (+/-{perr[0]:.2g}), "
            f"B={B_fit:.6g} (+/-{perr[1]:.2g}), C={C_fit:.6g} (+/-{perr[2]:.2g}), "
            f"{fit_limit_str(A_fit, B_fit, C_fit, perr[2])}",
            flush=True,
        )
    except RuntimeError as e:
        print(f"Warning: power-law fit to irreversible data failed to converge: {e}", flush=True)
else:
    print("Warning: not enough irreversible r>0 points to fit y=Ax^B+C (need at least 3).", flush=True)

### Power-law fit of reversible avg S/Nk vs radius (y = A x^B), excluding r=0
rev_r_arr_fit = np.array(rev_plotted_radii)
rev_avg_arr_fit = np.array(avg_Sk)
rev_fit_mask = rev_r_arr_fit > 0

if np.sum(rev_fit_mask) >= 3:
    try:
        popt_rev, pcov_rev = curve_fit(
            power_law,
            rev_r_arr_fit[rev_fit_mask],
            rev_avg_arr_fit[rev_fit_mask],
            p0=[1.0, -1.0, float(np.mean(rev_avg_arr_fit[rev_fit_mask]))],
            maxfev=10000,
        )
        A_fit_rev, B_fit_rev, C_fit_rev = popt_rev
        perr_rev = np.sqrt(np.diag(pcov_rev))

        r_fit_x_rev = np.linspace(rev_r_arr_fit[rev_fit_mask].min(), max(rev_r_arr_fit[rev_fit_mask].max(), 100), 200)
        r_fit_y_rev = power_law(r_fit_x_rev, A_fit_rev, B_fit_rev, C_fit_rev)

        fig2.add_trace(
            go.Scatter(
                x=r_fit_x_rev,
                y=r_fit_y_rev,
                mode='lines',
                name=f"Rev fit: y={A_fit_rev:.4g}x^{B_fit_rev:.4g}+{C_fit_rev:.4g} ({fit_limit_str(A_fit_rev, B_fit_rev, C_fit_rev, perr_rev[2])})",
                line=dict(color='purple', dash='dash'),
                hovertemplate='<b>Reversible Fit</b><br>Radius: %{x:.2f}<br>S/Nk: %{y:.6f}<extra></extra>',
            ),
            row=1, col=1,
        )

        print(
            f"Reversible power-law fit (r>0): A={A_fit_rev:.6g} (+/-{perr_rev[0]:.2g}), "
            f"B={B_fit_rev:.6g} (+/-{perr_rev[1]:.2g}), C={C_fit_rev:.6g} (+/-{perr_rev[2]:.2g}), "
            f"{fit_limit_str(A_fit_rev, B_fit_rev, C_fit_rev, perr_rev[2])}",
            flush=True,
        )
    except RuntimeError as e:
        print(f"Warning: power-law fit to reversible data failed to converge: {e}", flush=True)
else:
    print("Warning: not enough reversible r>0 points to fit y=Ax^B+C (need at least 3).", flush=True)

### Entropy Difference (ppm) vs Radius
# Align the two series on common radii, then express the gap in parts per million.
rev_r_arr = np.array(rev_plotted_radii)
irr_r_arr = np.array(irr_plotted_radii)
common_radii = np.intersect1d(rev_r_arr, irr_r_arr)

if len(common_radii) > 0:
    rev_idx = np.searchsorted(rev_r_arr, common_radii)
    irr_idx = np.searchsorted(irr_r_arr, common_radii)

    rev_vals  = avg_Sk[rev_idx]
    irr_vals  = irr_avg_Sk[irr_idx]
    rev_sems  = SEM[rev_idx]
    irr_sems  = irr_SEM[irr_idx]

    diff_ppm     = (irr_vals - rev_vals) * 1e6
    diff_sem_ppm = np.sqrt(irr_sems**2 + rev_sems**2) * 1e6

    fig2.add_trace(
        go.Scatter(
            x=common_radii,
            y=diff_ppm,
            error_y=dict(type='data', array=diff_sem_ppm),
            name="Irr − Rev",
            mode='lines+markers',
            line=dict(color='#E23CB1'),
            marker=dict(size=6),
            hovertemplate=(
                '<b>Irr − Rev</b><br>'
                'Radius: %{x}<br>'
                'Δ S/Nk: %{y:.2f} ppm<extra></extra>'
            ),
        ),
        row=1, col=2,
    )

    # Zero-reference line so the sign of the gap is obvious
    fig2.add_hline(y=0, line_width=1, line_dash="dot", line_color="gray", row=1, col=2)

fig2.update_xaxes(title_text="Local bath size", row=1, col=2)
fig2.update_yaxes(title_text="(S<sub>irr</sub>−S<sub>rev</sub>)/Nk [ppm]", row=1, col=2)

fig2.update_layout(title_text=f"Lattice Size: {n}")

# --- Stats box: total average entropy for r>0, irr vs rev, and their difference ---
# Filter to r>0 entries only (index 0 in the plotted_radii lists corresponds to r=0 if present)
def _total_avg_for_r_gt0(plotted_radii, avg_sk_array):
    """Return the mean of avg_sk_array entries where the corresponding radius > 0."""
    if len(plotted_radii) == 0 or len(avg_sk_array) == 0:
        return float('nan')
    mask = np.array(plotted_radii) > 0
    vals = avg_sk_array[mask]
    return float(np.mean(vals)) if len(vals) > 0 else float('nan')

total_avg_rev = _total_avg_for_r_gt0(rev_plotted_radii, avg_Sk)
total_avg_irr = _total_avg_for_r_gt0(irr_plotted_radii, irr_avg_Sk)
delta = total_avg_irr - total_avg_rev

def _fmt(val):
    return f"{val:.4f}" if not np.isnan(val) else "N/A"

stats_lines = [
    "<b>Total avg S/Nk (r > 0)</b>",
    f"  Rev:  {_fmt(total_avg_rev)}",
    f"  Irr:  {_fmt(total_avg_irr)}",
    f"  Δ (Irr − Rev):  {_fmt(delta)}",
]
stats_text = "<br>".join(stats_lines)

fig2.add_annotation(
    text=stats_text,
    xref="paper", yref="paper",
    x=1.01, y=0.45,
    xanchor="left", yanchor="top",
    showarrow=False,
    align="left",
    font=dict(size=12, family="monospace"),
    bgcolor="rgba(255,255,255,0.85)" if not is_dark_mode else "rgba(30,30,30,0.85)",
    bordercolor="rgba(0,0,0,0.3)" if not is_dark_mode else "rgba(200,200,200,0.3)",
    borderwidth=1,
    borderpad=8,
)

print("Serializing fig2 to HTML...", flush=True)
fig2.show()

### Figure 2D recreation: ln(p_n/p_(n+1)) and equivalent kT/J moving averages vs sweeps
fig3.update_xaxes(title_text="t (sweeps)", row=1, col=1)
fig3.update_xaxes(title_text="t (sweeps)", row=1, col=2)
fig3.update_yaxes(title_text="ln(p<sub>n</sub> / p<sub>n+1</sub>)", row=1, col=1)
fig3.update_yaxes(title_text="kT/J", row=1, col=2)
fig3.update_layout(
    title_text=f"Lattice Size: {n} — Ratio of Adjacent Demon-Level Occupation Probabilities (cf. Chamberlin 2024, Fig. 2D)",
)
fig3.add_annotation(
    text="Color: radius/dynamics (see legend) &nbsp;|&nbsp; Right panel: kT/J = 1/ln(p<sub>n</sub>/p<sub>n+1</sub>) \u2014 same data as an effective temperature",
    xref="paper", yref="paper",
    x=0.5, y=-0.15,
    xanchor="center", yanchor="top",
    showarrow=False,
    font=dict(size=11, color="rgba(100,100,100,0.8)" if is_dark_mode else "rgba(120,120,120,0.9)"),
)
print("Serializing fig3 to HTML...", flush=True)
fig3.show()