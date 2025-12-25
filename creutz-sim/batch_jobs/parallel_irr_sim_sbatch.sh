#!/bin/bash
#SBATCH --job-name=parallel_sim_irreversible
#SBATCH --output=logs/parallel_sim_irr_%j.out
#SBATCH --error=logs/parallel_sim_irr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8GB
#SBATCH --time=7-00:00:00

# Parallel irreversible Creutz demon simulation SLURM batch script
# Optimized for multi-core execution on HPC clusters

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load any required modules (uncomment and modify as needed)
# module load python/3.11
# module load numpy/2.0
# module load scipy/1.14

# Activate virtual environment
source venv/bin/activate

# Change to simulation directory
cd creutz-sim

# Run parallel simulation with explicit core count
# Uses SLURM_CPUS_PER_TASK to match allocated resources
python parallel_irr_sim.py \
    --n 1000000 \
    --s 5000 \
    --r 11 \
    --m 5 \
    --cores $SLURM_CPUS_PER_TASK \
    --validate off

echo "Job completed at: $(date)"
