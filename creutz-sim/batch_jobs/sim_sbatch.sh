#!/bin/bash

# Usage: sbatch sim_sbatch.sh [args for sim.py]

#SBATCH -N 1  		# number of nodes
#SBATCH -c 1  		# number of "tasks" (cores)
#SBATCH --mem=1G        # GigaBytes of memory required (per node)
#SBATCH -t 7-00:00:00   # time in d-hh:mm:ss
#SBATCH -p public	# partition
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
##SBATCH --mail-user=wember@asu.edu # Mail-to address

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Activate virtual environment
source "$REPO_ROOT/venv/bin/activate"

# Run simulation, add any additional arguments passed to the script
python "$REPO_ROOT/creutz-sim/sim.py" $@
