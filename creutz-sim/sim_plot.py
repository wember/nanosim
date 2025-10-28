import pandas as pd
import numpy as np
from scipy.special import loggamma as logg
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.templates.default = "plotly_white"

df = pd.read_csv('sim_data.csv')

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
