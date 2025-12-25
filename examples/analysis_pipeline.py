"""
End-to-end analysis pipeline example.

Shows the complete workflow:
1. Run simulation
2. Save data
3. Load and analyze results
4. Generate plots
"""
import sys
import os
import csv
import numpy as np
from scipy.special import loggamma as logg

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'creutz-sim'))

from inferno import Inferno

# Configuration
N = 100
R = 3
sweeps = 50

print(f"Running simulation: N={N}, R={R}, sweeps={sweeps}\n")

# Initialize
x = Inferno(N, R)

# Data collection
data = []

# Forward phase
print("Forward phase...")
for i in range(sweeps):
    for _ in range(N):
        x.demon_move()
    
    state = x.get_validated_state()
    
    # Calculate entropy
    Sk = logg(state['E_demon_sum'] + N) - logg(state['E_demon_sum'] + 1) - logg(N)
    N0_exp = max(int(state['bond_count'][1]), 1)
    Su = logg(N + 1) + N0_exp * np.log(2) - \
         (logg(N - int(state['bond_count'][1]) - int(state['bond_count'][2]) + 1) + 
          logg(int(state['bond_count'][1]) + 1) + 
          logg(int(state['bond_count'][2]) + 1))
    S = (Sk + Su) / N
    
    data.append({
        't': i,
        'K': state['K'],
        'U': state['U'],
        'S': S,
        'phase': 'forward'
    })

# Reverse phase
print("Reverse phase...")
for i in range(sweeps):
    for _ in range(N):
        x.demon_reverse()
    
    state = x.get_validated_state()
    
    Sk = logg(state['E_demon_sum'] + N) - logg(state['E_demon_sum'] + 1) - logg(N)
    N0_exp = max(int(state['bond_count'][1]), 1)
    Su = logg(N + 1) + N0_exp * np.log(2) - \
         (logg(N - int(state['bond_count'][1]) - int(state['bond_count'][2]) + 1) + 
          logg(int(state['bond_count'][1]) + 1) + 
          logg(int(state['bond_count'][2]) + 1))
    S = (Sk + Su) / N
    
    data.append({
        't': sweeps + i,
        'K': state['K'],
        'U': state['U'],
        'S': S,
        'phase': 'reverse'
    })

# Save to CSV
output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'examples')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'pipeline_example.csv')

with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['t', 'K', 'U', 'S', 'phase'])
    writer.writeheader()
    writer.writerows(data)

print(f"\nData saved to: {output_file}")

# Basic analysis
K_values = [d['K'] for d in data]
U_values = [d['U'] for d in data]
S_values = [d['S'] for d in data]

print("\nSummary statistics:")
print(f"  <K> = {np.mean(K_values):.4f} ± {np.std(K_values):.4f}")
print(f"  <U> = {np.mean(U_values):.4f} ± {np.std(U_values):.4f}")
print(f"  <S/k> = {np.mean(S_values):.4f} ± {np.std(S_values):.4f}")

# Check reversibility
forward_final = data[sweeps-1]
reverse_final = data[-1]
print(f"\nReversibility check:")
print(f"  Forward final:  K={forward_final['K']:.4f}, U={forward_final['U']:.4f}")
print(f"  Reverse final:  K={reverse_final['K']:.4f}, U={reverse_final['U']:.4f}")
print(f"  Difference:     ΔK={abs(forward_final['K']-reverse_final['K']):.6f}")

print("\n✓ Pipeline complete!")
