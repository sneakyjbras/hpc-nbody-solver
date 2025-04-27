# N‑Body Solver – README

## How to Run the Program

To run the simulation, simply execute the provided run script:

```bash
./run.sh
```

This script compiles (if necessary) and then runs the simulation with the proper parameters.

## Folders

The following folders contain:
- `src/`: The source code of the project and its `Makefile`. The `main` function is found on the `main.cpp` file.
- `utils/`: Scripts that were used by us to measure the memory and performance of the program.
- `docs`: The project statement.
- `logs`: The examples given by the faculty and our best scenario on the `lab1p1` machine at **RNL**.

## Overview of Optimizations

This simulation has been optimized for high performance using several techniques:

### Sweep and Prune for Collision Detection

Instead of comparing every pair of particles in each cell, we sort the particles by their x-coordinate and then only compare particles until the x difference exceeds a fixed collision radius. This early-exit strategy—called sweep and prune—significantly reduces the number of comparisons when particles are spread out. More information about this can be found [here](https://en.wikipedia.org/wiki/Sweep_and_prune).

### BFS Implementation for Collision Grouping

When collisions are detected, we use a breadth-first search (BFS) to group colliding particles. BFS is chosen over recursion to avoid potential stack overflows with large collision groups and to provide more predictable, iterative performance.

### Computations Avoiding `sqrt`

In collision detection, we compare squared distances to a squared collision threshold. This avoids the expensive square root computation for every pair of particles. In gravitational force calculations, we use formulations that minimize `sqrt` calls—computing the inverse square root once per pair when needed.

## Order of Computations in Each Simulation Step

The simulation proceeds through the following steps in each time step:

1. **Compute Center of Mass**:
   For each grid cell, compute the center of mass of the particles in that cell.

2. **Compute In-Cell Gravitational Pull**:
   Calculate gravitational forces between particles within the same cell, using pairwise comparisons (optimized via squared distance checks).

3. **Compute Extra-Cell Gravitational Pull**:
   Use the computed center of mass from adjacent cells to compute gravitational pull on particles from outside their cell.

4. **Update Particle Positions**:
   Compute new positions for each particle based on their current velocity and gravitational forces.

5. **Update Particle Velocities**:
   Adjust particle velocities based on the computed forces.

6. **Check and Resolve Collisions**:
   Detect collisions using the optimized sweep and prune approach and group colliding particles using BFS. Colliding particles are then processed accordingly.

## Branches

The repository maintains dedicated branches for each parallelisation strategy:

### `main` (default ‑ serial)

The **main** branch contains the baseline, single‑threaded implementation. It is intended for correctness reference and for measuring raw speed‑ups of the parallel variants.

### `omp`

The **omp** branch parallelises the computation with OpenMP.

- **Parallel loops** on all compute‑intensive phases (forces, position & velocity updates, collision checks).
- **Scheduling**: static for regular loops, dynamic for collision resolution.
- **Thread‑local force accumulators** to avoid synchronisation overhead.
- **Build flags**: `-fopenmp` is added automatically by the branch‑local `Makefile`.
- **Performance**: Speed‑up is **≈ p** for uniform particle distributions and between **p ⁄ 2 – p** for Gaussian clusters (e.g. ~7.5× on 8 threads). A detailed technical report is available in `docs/omp_report.pdf`.

### `mpi`

The **mpi** branch distributes the simulation over multiple nodes using MPI.

- **Domain decomposition**: 2‑D block decomposition of the grid, halo‑exchange for border cells.
- **Non‑blocking communication** overlaps data exchange with computation.
- **Scalability**: Near‑linear speed‑up on up to 64 MPI ranks for uniform workloads; performance diminishes gracefully for highly clustered cases.

### `mpi‑omp` (hybrid)

The **mpi‑omp** branch combines MPI between nodes with OpenMP threading inside each rank.

- **Hybrid parallel regions**: Each rank launches an OpenMP team; MPI calls are confined to the master thread.
- **NUMA‑aware affinity** via `OMP_PLACES` and `OMP_PROC_BIND`.
- **Cluster job script**: A ready‑to‑use SLURM submission helper, `run_slurm.sh`, lives at the project root. Tweak the `#SBATCH` resource lines to match your cluster (nodes, tasks‑per‑node, OMP threads) and submit with `sbatch run_slurm.sh`.
- **Performance**:
  - **Super‑linear** for uniform distributions (e.g. **128×** on **32 MPI ranks × 8 threads** thanks to improved cache locality).
  - Between **p ⁄ 2 – p** for Gaussian distributions (e.g. ~55× on **8 MPI ranks × 8 threads**).
- A comprehensive technical report detailing experiments and analysis is provided in `docs/mpi_omp_report.pdf`.

---

## Conclusion

This project showcases how algorithmic optimisations coupled with multi‑level parallelism can scale an \(n\)-body particle simulation from a single core to large, multi‑node clusters.

* **Algorithmic strength first** – Sweep‑and‑prune collision detection, squared‑distance checks and BFS grouping deliver a solid serial baseline.
* **Parallel breadth next** – OpenMP pushes shared‑memory performance close to the theoretical limit \(≈ p\), MPI extends scaling across nodes, and the hybrid branch attains super‑linear gains when cache effects align.
* **Reproducibility** – Every branch builds and runs through the same `./run.sh` wrapper, while benchmark logs and two technical reports in `docs/` make it easy to verify and compare results.

Taken together, the code, documentation and benchmark artefacts form a compact playground for experimenting with high‑performance simulation techniques and exploring the trade‑offs between shared‑memory, distributed‑memory and hybrid approaches.


