#!/bin/bash

# ==============================================================================
# SLURM Job Configuration Directives
# ==============================================================================
# These lines starting with #SBATCH are NOT regular comments!
# They are read by SLURM workload manager to configure the job.
# MUST appear at the top of the file before any bash commands.
# ==============================================================================

#SBATCH --job-name=parallel_sim_reversible  # Name displayed in job queue
#SBATCH --output=logs/parallel_sim_rev_%j.out  # stdout file (%j = SLURM job ID)
#SBATCH --error=logs/parallel_sim_rev_%j.err   # stderr file (%j = SLURM job ID)
#SBATCH --nodes=1                    # Request 1 compute node
#SBATCH --ntasks=1                   # Run 1 task (not using MPI)
#SBATCH --cpus-per-task=16           # Allocate 16 CPU cores per task
#SBATCH --mem=8GB                    # Request 8GB of RAM
#SBATCH --time=7-00:00:00            # Max runtime: 7 days (format: D-HH:MM:SS)

# ==============================================================================
# Parallel reversible Creutz demon simulation SLURM batch script
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

# Set SLURM_CPUS_PER_TASK to number of physical cores if not set (for local testing)
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
    if command -v sysctl >/dev/null 2>&1; then
        export SLURM_CPUS_PER_TASK=$(sysctl -n hw.physicalcpu)
    else
        export SLURM_CPUS_PER_TASK=4
    fi
    echo "[Local] SLURM_CPUS_PER_TASK not set, defaulting to $SLURM_CPUS_PER_TASK cores."
fi


# Run parallel simulation with JIT optimization
python parallel_sim.py \
    --jit \
    --n 1000000 \
    --s 5000 \
    --r 11 \
    --m 5 \
    --cores $SLURM_CPUS_PER_TASK \
    --validate off

echo "Job completed at: $(date)"
