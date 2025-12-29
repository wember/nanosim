#!/bin/bash

# ==============================================================================
# SLURM Job Configuration Directives
# ==============================================================================
# These lines starting with #SBATCH are NOT regular comments!
# They are read by SLURM workload manager to configure the job.
# MUST appear at the top of the file before any bash commands.
# ==============================================================================

#SBATCH --job-name=parallel_sim_irreversible  # Name displayed in job queue
#SBATCH --output=logs/parallel_sim_irr_%j.out  # stdout file (%j = SLURM job ID)
#SBATCH --error=logs/parallel_sim_irr_%j.err   # stderr file (%j = SLURM job ID)
#SBATCH --nodes=1                    # Request 1 compute node
#SBATCH --ntasks=1                   # Run 1 task (not using MPI)
#SBATCH --cpus-per-task=16           # Allocate 16 CPU cores per task
#SBATCH --mem=8GB                    # Request 8GB of RAM
#SBATCH --time=7-00:00:00            # Max runtime: 7 days (format: D-HH:MM:SS)

# ==============================================================================
# Parallel irreversible Creutz demon simulation SLURM batch script
# Optimized for multi-core execution on HPC clusters
# ==============================================================================

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load any required modules (uncomment and modify as needed)
# module load python/3.11
# module load numpy/2.0
# module load scipy/1.14


# Resolve project root (two levels up from this script)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Activate virtual environment
if [ -f "$PROJECT_ROOT/venv/bin/activate" ]; then
    source "$PROJECT_ROOT/venv/bin/activate"
else
    echo "[SBATCH WRAPPER] ERROR: venv/bin/activate not found at $PROJECT_ROOT/venv/bin/activate"
    exit 1
fi

# Change to simulation directory
cd "$PROJECT_ROOT/creutz-sim"

# Run parallel irreversible simulation (JIT enabled by default)
# $SLURM_CPUS_PER_TASK is automatically set by SLURM (16 cores from --cpus-per-task above)
python parallel_sim_irr.py \
    --n 1000000 \
    --s 5000 \
    --r 11 \
    --m 5 \
    --cores $SLURM_CPUS_PER_TASK \
    --validate off

echo "Job completed at: $(date)"
