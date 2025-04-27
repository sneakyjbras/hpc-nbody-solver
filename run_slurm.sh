#!/usr/bin/env bash
#SBATCH --job-name=bombardiro_crocodilo
#SBATCH --output=logs_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# Number of nodes requested from sbatch
NP=${SLURM_JOB_NUM_NODES}

echo "======================================"
echo "Running with $NP MPI process(es) across $NP node(s)"
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

for CMD in "${COMMANDS[@]}"; do
    echo -e "\nCommand: srun -n $NP $CMD"
    srun -n "$NP" $CMD
done

echo -e "\nFinished run for $NP MPI process(es).\n"
echo "All tests complete."
