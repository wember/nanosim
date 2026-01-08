import pandas as pd
import numpy as np
from scipy.special import loggamma as logg
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"

fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.2)

file_names = ['sim_data.csv',
              'sim_data_r1.csv',
              'sim_data_r2.csv',
              'sim_data_r3.csv',
              'sim_data_r4.csv',
              'sim_data_r5.csv',
              'sim_data_r6.csv',
              'sim_data_r7.csv',
              'sim_data_r8.csv',
              'sim_data_r9.csv',
              'sim_data_r10.csv']

data_names = ['radius = 0',
              'radius = 1',
              'radius = 2',
              'radius = 3',
              'radius = 4',
              'radius = 5',
              'radius = 6',
              'radius = 7',
              'radius = 8',
              'radius = 9',
              'radius = 10']

j=0
for file in file_names:
    df = pd.read_csv(file)

    n = df['n'][0] # Lattice size
    t = df['t']    # step number
    U = df['U']    # potential (lattice) energy
    K = df['K']    # kinetic (demon) energy
    Nx = df['Nx']  # anti-aligned spins
    N0 = df['N0']  # broken bonds

    df['N0_exp'] = df['N0'].replace(0, 1) # If no broken bonds, 2^N0+1 in Su equation
    N0_exp = df['N0_exp']

    Sk = lambda N, K: logg(K + N) - logg(K+1) - logg(N) # N == lattice size, K == kinetic energy
    Su = lambda N, N0, Nx: logg(N+1) + np.log(2**N0_exp) - (logg(N-N0-Nx+1) + logg(N0+1) + logg(Nx+1)) # N == lattice size, N0 == broken bonds, Nx == bonds between anti-aligned spins

    total_entropy = (Sk(n, K) + Su(n, N0, Nx))/n

    ### Demon Energy ###
    fig.add_trace(go.Scatter(x=t, y=K, name=data_names[j]),row=1, col=1)
    fig.update_xaxes(title_text="Sweeps", row=1, col=1)
    fig.update_yaxes(title_text="Demon Energy", row=1, col=1)
    fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=1)

    ### Lattice Energy ###
    fig.add_trace(go.Scatter(x=t, y=U, name=data_names[j]),row=1, col=2)
    # fig.add_trace(go.Scatter(x=tR, y=KR, showlegend=False),row=1, col=2)
    fig.update_xaxes(title_text="Sweeps", row=1, col=2)
    fig.update_yaxes(title_text="Lattice Temp", row=1, col=2)
    fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=2)

    ### Entropy ###
    fig.add_trace(go.Scatter(x=t, y=total_entropy, name=data_names[j]),row=1, col=3)
    fig.update_xaxes(title_text="Sweeps", row=1, col=3)
    fig.update_yaxes(title_text="S/Nk", row=1, col=3)
    fig.add_vline(x=len(df)//2, line_width=1, line_dash="dash", line_color="Red", row=1, col=3)

    j+=1

fig.update_layout(title_text=f"Lattice Size: {n}")
fig.show()
