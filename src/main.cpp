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
  MPI_Init(&argc, &argv);

  // Parse command-line arguments.
  const int32_t seed = static_cast<int32_t>(std::stol(argv[1], nullptr, 10));
  const double side = std::strtod(argv[2], nullptr);
  const int32_t ncside =
      static_cast<int32_t>(std::strtol(argv[3], nullptr, 10));
  const uint64_t nPart =
      static_cast<uint64_t>(std::strtoll(argv[4], nullptr, 10));
  const int32_t tsteps =
      static_cast<int32_t>(std::strtol(argv[5], nullptr, 10));

  // Create simulation instance with parsed parameters.
  Simulation::Parsim parsim(seed, side, ncside, nPart, tsteps, MPI_COMM_WORLD);

  if (parsim.getWorldSize() < ncside) {
    // Only rank 0 prints the result.
  } else {
    std::cerr << "Total procs: " << parsim.getWorldSize()
              << "; should be less than ncside: " << ncside << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return -1;
  }

  // Start timing the simulation.
  double startTime = omp_get_wtime();

  // Setup simulation.
  parsim.getAllNeighbors();
  parsim.allocateGrid();
  parsim.populateGrid();
  parsim.createMPIType();
  parsim.initializeCommunication();

  std::unordered_set<Particle *> visited; // Track visited particles.
  std::queue<Particle *> bfsQueue;        // Queue for BFS traversal.
  uint64_t totalCollisions = 0;

  totalCollisions += parsim.simulate(visited, bfsQueue);

  if (parsim.getWorldSize() < ncside) {
    // Only rank 0 prints the result.
    if (!parsim.getWorldRank()) {
      std::cout << std::fixed << std::setprecision(3) << parsim.getParticle0().x
                << " " << parsim.getParticle0().y << std::endl;
      std::cout << totalCollisions << std::endl;
      double execTime = omp_get_wtime() - startTime;
      std::cerr << std::fixed << std::setprecision(1) << execTime << "s"
                << std::endl;
    }
  } else {
    std::cerr << "Total procs: " << parsim.getWorldSize()
              << "; should be less than ncside: " << ncside << std::endl;
  }

  MPI_Finalize();
  return 0;
}
