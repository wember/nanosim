#!/bin/bash

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

source /home/wember/2025thesis/my_venv/bin/activate

python /home/wember/2025thesis/nanosim/creutz-sim/sim.py
python /home/wember/2025thesis/nanosim/creutz-sim/irr_sim.py

/home/wember/2025thesis/nanosim/data
