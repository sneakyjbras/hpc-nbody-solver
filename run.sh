#!/bin/bash
# Redirect all output (stdout and stderr) to the file "logs"
exec > logs 2>&1

# Define the MPI process counts to test
PROCS=(1 2 4 8 16 32 64)
# Define the OpenMP thread counts to test (e.g., 1 to 8)
THREADS=(1)

for NP in "${PROCS[@]}"; do
    echo "======================================"
    echo "Running with MPI processes: $NP"
    echo "======================================"

    # Change into the "src" folder
    if ! cd src; then
        echo "Error: Unable to change into src directory."
        exit 1
    fi

    # Clean parsim
    echo "Cleaning parsim..."
    make clean

    # Build parsim
    echo "Building parsim..."
    make
    if [ $? -ne 0 ]; then
        echo "Build failed. Please check the errors above."
        exit 1
    fi

    echo "Build succeeded. Running parsim commands..."

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 1 5000 100 1000000 4"
        mpirun -np $NP --oversubscribe ./parsim 1 5000 100 1000000 4

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 1 5000 100 1000000 100"
        mpirun -np $NP --oversubscribe ./parsim 1 5000 100 1000000 100

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 1 5000 20 1000000 10"
        mpirun -np $NP --oversubscribe ./parsim 1 5000 20 1000000 10

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 1 1000 3 10000 10000"
        mpirun -np $NP --oversubscribe ./parsim 1 1000 3 10000 10000

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 3 5000 50 1000000 300"
        mpirun -np $NP --oversubscribe ./parsim 3 5000 50 1000000 300

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 3 5000 50 1000000 500"
        mpirun -np $NP --oversubscribe ./parsim 3 5000 50 1000000 500

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim -1 1000 30 100000 1000"
        mpirun -np $NP --oversubscribe ./parsim -1 1000 30 100000 1000

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 12672 0.05 3 10 10"
        mpirun -np $NP --oversubscribe ./parsim 12672 0.05 3 10 10

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 5893 0.05 3 10 10"
        mpirun -np $NP --oversubscribe ./parsim 5893 0.05 3 10 10

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 8555 0.05 3 10 10"
        mpirun -np $NP --oversubscribe ./parsim 8555 0.05 3 10 10

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim 12 100 5 10000 10000"
        mpirun -np $NP --oversubscribe ./parsim 12 100 5 10000 10000

        echo -e "\nCommand: mpirun -np $NP --oversubscribe ./parsim -11 3500 20 500000 10"
        mpirun -np $NP --oversubscribe ./parsim -11 3500 20 500000 10

        echo -e "\nAll commands executed for OMP_NUM_THREADS: $NT with MPI processes: $NP."
    done

    # Return to the original directory for the next iteration
    cd ..
done

echo "Script finished. Check the logs file for details."
