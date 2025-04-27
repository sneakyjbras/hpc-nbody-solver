#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "parsim.hpp"

int32_t main(int32_t argc, char *argv[]) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0]
              << " <seed> <side> <ncside> <nPart> <tsteps>" << std::endl;
    return 1;
  }

  // Enable nested parallelism if required by simulation routines.
  omp_set_nested(true);

  std::unordered_set<Simulation::Particle *>
      visited; // Set to track visited particles.
  std::queue<Simulation::Particle *> bfsQueue; // Queue for BFS traversal.

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

  Simulation::Parsim parsim(seed, side, ncside, nPart, tsteps);

  // Start timing the simulation.
  double startTime = omp_get_wtime();

  uint64_t totalCollisions = 0;
#pragma omp parallel private(visited, bfsQueue) num_threads(nThreads)          \
    reduction(+ : totalCollisions)
  {
    // Compute this thread's linear cell range
    int32_t tid = omp_get_thread_num();
    int32_t totalCells = ncside * ncside;
    int32_t cellsPerThread = (totalCells + nThreads - 1) / nThreads;
    int32_t idxStart = tid * cellsPerThread;
    int32_t idxEnd = std::min(totalCells, idxStart + cellsPerThread);

    parsim.allocateGrid();
    parsim.populateGrid();

    // Run the simulation for this subrange of cells
    totalCollisions += parsim.simulate(visited, bfsQueue, idxStart, idxEnd);
  }

  double execTime = omp_get_wtime() - startTime;
  std::cout << std::fixed << std::setprecision(3) << parsim.get_particle(0).x
            << " " << parsim.get_particle(0).y << std::endl;
  std::cout << totalCollisions << std::endl;
  std::cerr << std::fixed << std::setprecision(1) << execTime << "s"
            << std::endl;
  return 0;
}
