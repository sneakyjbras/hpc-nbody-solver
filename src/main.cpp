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

  Simulation::Parsim parsim(seed, side, ncside, nPart, tsteps);
  // Start timing the simulation.
  double startTime = omp_get_wtime();

  parsim.allocateGrid();
  parsim.populateGrid();

  uint64_t totalCollisions = 0;

  // Compute this thread's linear cell range
  int32_t total_cells = ncside * ncside;
  int32_t idxStart = 0, idxEnd = total_cells;

  // Run the simulation for this subrange of cells
  totalCollisions += parsim.simulate(visited, bfsQueue, idxStart, idxEnd);

  double exec_time = omp_get_wtime() - startTime;
  std::cout << std::fixed << std::setprecision(3) << parsim.get_particle(0).x
            << " " << parsim.get_particle(0).y << std::endl;
  std::cout << totalCollisions << std::endl;
  std::cerr << std::fixed << std::setprecision(1) << exec_time << "s"
            << std::endl;
  return 0;
}
