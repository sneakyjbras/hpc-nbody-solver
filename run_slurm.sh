#!/usr/bin/env bash
#SBATCH --job-name=bombardiro_crocodilo
#SBATCH --output=logs_%j.txt
#SBATCH --exclusive
# Note: Do not set --cpus-per-task here. Instead, pass it via the sbatch command line.
# Example submission: sbatch -N <node_count> --cpus-per-task=<OMP_THREADS> run_slurm.sh <OMP_THREADS>

# Check for the required argument for OpenMP threads.
if [ -z "$1" ]; then
    echo "Usage: sbatch -N <node_count> --cpus-per-task=<OMP_THREADS> run_slurm.sh <OMP_THREADS>"
    exit 1
fi

OMP_THREADS=$1
echo "Using OMP_THREADS=${OMP_THREADS}"

# Determine number of MPI processes (assuming one per node).
NP=${SLURM_JOB_NUM_NODES}

echo "======================================"
echo "Running with ${NP} MPI process(es) across ${NP} node(s)"
echo "OpenMP threads per process: ${OMP_THREADS}"
echo "======================================"

declare -a COMMANDS=(
    "./parsim/src/parsim 1 5000 100 1000000 4"
    "./parsim/src/parsim 1 5000 100 1000000 100"
    "./parsim/src/parsim 1 5000 20 1000000 10"
    "./parsim/src/parsim 1 1000 3 10000 10000"
    "./parsim/src/parsim 3 5000 50 1000000 300"
    "./parsim/src/parsim 3 5000 50 1000000 500"
    "./parsim/src/parsim -1 1000 30 100000 1000"
    "./parsim/src/parsim 12672 0.05 3 10 10"
    "./parsim/src/parsim 5893 0.05 3 10 10"
    "./parsim/src/parsim 8555 0.05 3 10 10"
    "./parsim/src/parsim 12 100 5 10000 10000"
    "./parsim/src/parsim -11 3500 20 500000 10"
)

# Loop over commands and execute each with the proper environment variable setting.
for CMD in "${COMMANDS[@]}"; do
    echo -e "\nCommand: OMP_NUM_THREADS=${OMP_THREADS} srun -n ${NP} ${CMD}"
    OMP_NUM_THREADS=${OMP_THREADS} srun -n "${NP}" $CMD
done

echo -e "\nFinished run for ${NP} MPI process(es) with ${OMP_THREADS} OpenMP thread(s) each."
echo "All tests complete."
