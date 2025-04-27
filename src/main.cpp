#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "parsim.hpp"

int32_t main(int32_t argc, char *argv[]) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0]
              << " <seed> <side> <ncside> <nPart> <tsteps>" << std::endl;
    return 1;
  }

  // Initialize MPI communications.
  int32_t provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    std::cerr << "Error: MPI does not support MPI_THREAD_FUNNELED."
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Parse command-line arguments.
  const int32_t seed = static_cast<int32_t>(std::stol(argv[1], nullptr, 10));
  const double side = std::strtod(argv[2], nullptr);
  const int32_t ncside =
      static_cast<int32_t>(std::strtol(argv[3], nullptr, 10));
  const uint64_t nPart =
      static_cast<uint64_t>(std::strtoll(argv[4], nullptr, 10));
  const int32_t tsteps =
      static_cast<int32_t>(std::strtol(argv[5], nullptr, 10));
  const int32_t nThreads = static_cast<int32_t>(omp_get_max_threads());

  // Create simulation instance with parsed parameters.
  Simulation::Parsim parsim(seed, side, ncside, nPart, tsteps, nThreads,
                            MPI_COMM_WORLD);

  if (parsim.getWorldSize() < ncside) {
    // Only rank 0 prints the result.
  } else {
    std::cerr << "Total procs: " << parsim.getWorldSize()
              << "; should be less than ncside: " << ncside << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return -1;
  }

  // Compute the subdomain for each process and get its neighbors.
  parsim.computeSubdomain();

  // Initialize particles.
  parsim.initParticles();

  // Start timing the simulation.
  double startTime = omp_get_wtime();

  std::unordered_set<Particle *> visited; // Track visited particles.
  std::queue<Particle *> bfsQueue;        // Queue for BFS traversal.
  uint64_t totalProcCollisions = 0;

#pragma omp parallel private(visited, bfsQueue) num_threads(nThreads)          \
    reduction(+ : totalProcCollisions)
  {
    // Setup simulation.
    parsim.allocateGrid();
    parsim.populateGrid();
#pragma omp master
    {
      parsim.getAllNeighbors();
      parsim.createMPIType();
      parsim.initializeCommunication();
    }

    uint16_t tid = omp_get_thread_num();
    parsim.computeThreadIndexes(tid);
    totalProcCollisions += parsim.simulate(visited, bfsQueue);
  }
  parsim.gatherParticle0();
  uint64_t totalParsimCollisions =
      parsim.reduceCollisionCount(totalProcCollisions);

  if (!parsim.getWorldRank()) {
    std::cout << std::fixed << std::setprecision(3) << parsim.getParticle0().x
              << " " << parsim.getParticle0().y << std::endl;
    std::cout << totalParsimCollisions << std::endl;
    double execTime = omp_get_wtime() - startTime;
    std::cerr << std::fixed << std::setprecision(1) << execTime << "s"
              << std::endl;
  }

  MPI_Finalize();
  return 0;
}
