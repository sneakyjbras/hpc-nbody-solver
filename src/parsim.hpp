#ifndef PARSIM_HPP
#define PARSIM_HPP

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <memory>
#include <mpi.h>
#include <omp.h>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ds.hpp"
#include "flatGrid2D.hpp"
#include "parsim_utils.hpp"
#include "pointCalculations.hpp"

namespace Simulation {

/**
 * @brief The Parsim class encapsulates the functionality of a parallel
 * simulation. It manages particles, cells, and MPI communications,
 * performing the simulation through a series of timesteps that include grid
 * initialization, gravitational computations, particle updates, and collision
 * detection.
 */
class Parsim {
public:
  //--------------------------------------------------------------------------
  // Constructors and Destructor
  //--------------------------------------------------------------------------

  /**
   * @brief Constructs a Parsim simulation instance.
   *
   * Initializes simulation parameters such as random seed, simulation domain
   * side length, grid dimensions, number of particles, timesteps, and MPI
   * communicator.
   *
   * @param inputSeed Random seed for initialization.
   * @param side Side length of the simulation domain.
   * @param ncside Number of cells per side of the grid.
   * @param nPart Total number of particles.
   * @param tsteps Number of timesteps to simulate.
   * @param worldComm MPI communicator for the simulation.
   */
  Parsim(int32_t inputSeed, double side, int32_t ncside, uint64_t nPart,
         int32_t tsteps, MPI_Comm worldComm);

  /**
   * @brief Destructor for Parsim.
   */
  ~Parsim();

  //--------------------------------------------------------------------------
  // Public Simulation Control Methods
  //--------------------------------------------------------------------------

  /**
   * @brief Executes the simulation over a series of timesteps.
   *
   * Iterates through timesteps and, for each, performs:
   * - Center-of-mass computation
   * - Gravitational force calculations (intra- and extracell)
   * - Particle position and velocity updates
   * - Collision detection and resolution
   * Synchronization among threads is managed via OpenMP barriers.
   *
   * @param visited Set to track visited particles during collision grouping.
   * @param q Queue used for BFS grouping of colliding particles.
   * @return Total number of collisions detected over all timesteps.
   */
  uint64_t simulate(std::unordered_set<Particle *> &visited,
                    std::queue<Particle *> &q);

  /**
   * @brief Allocates memory for the simulation grid.
   *
   * Initializes the grid structure based on the current simulation parameters.
   */
  void allocateGrid();

  /**
   * @brief Populates the simulation grid with particles.
   *
   * Distributes particles across the grid cells according to initial
   * conditions.
   */
  void populateGrid();

  /**
   * @brief Gathers information about neighboring domains.
   *
   * Determines the neighboring MPI processes or grid cells required for ghost
   * exchange.
   */
  void getAllNeighbors();

  /**
   * @brief Initializes MPI communication parameters.
   *
   * Sets up MPI communicators and other communication-related parameters.
   */
  void initializeCommunication();

  /**
   * @brief Creates custom MPI data types.
   *
   * Defines MPI data types for sending and receiving center-of-mass and
   * particle data.
   */
  void createMPIType();

  /**
   * @brief Retrieves the first (global) particle.
   *
   * @return Particle representing the first particle in the global context.
   */
  inline Particle getParticle0() { return globalParticle0; }

  /**
   * @brief Returns the MPI rank of the current process.
   *
   * @return Rank of the current process.
   */
  inline int32_t getWorldRank() { return worldRank; }

  /**
   * @brief Returns the total number of MPI processes.
   *
   * @return Total MPI process count.
   */
  inline int32_t getWorldSize() { return worldSize; }

  //--------------------------------------------------------------------------
  // Private Simulation Core Functions
  //--------------------------------------------------------------------------

private:
  //--------------------------------------------------------------------------
  // Simulation Parameters and MPI Variables
  //--------------------------------------------------------------------------

  mutable int32_t inputSeed;
  const double side;    ///< Domain side length.
  const int32_t ncside; ///< Number of cells per side.
  const uint64_t nPart; ///< Total number of particles.
  const int32_t tsteps; ///< Number of timesteps.
  Particle localParticle0, globalParticle0;
  bool has_p0 = false;

  // MPI communicator information.
  const int32_t periods[2] = {1, 1}; // 2D toroidal wrapping.
  int32_t worldRank, worldSize;
  int32_t cartRank, cartSize;
  int32_t dims[2] = {0, 0}, coords[2];
  int32_t outdegree, indegree;
  int32_t weighted;
  MPI_Comm worldComm, cartComm, graphComm;

  // MPI graph communication.
  std::vector<int32_t> weights, ranks;

  // MPI communication buffers for center-of-mass exchange.
  std::vector<CenterOfMass> sendComBuf;
  std::vector<int32_t> sendCounts;
  std::vector<int32_t> sComDispls;
  std::vector<CenterOfMass> recvComBuf;
  std::vector<int32_t> recvCounts;
  std::vector<int32_t> rComDispls;

  // Particle grouping by neighbor.
  std::unordered_map<int32_t, std::vector<Particle>> particleNeighborMap;
  std::unordered_map<int32_t, std::vector<GhostDirection>> neighborDirectionMap;

  // Communication buffers for particles.
  std::vector<Particle> parSendBuffer;
  std::vector<int32_t> parSendCounts;
  std::vector<int32_t> parSDispls;
  std::vector<Particle> parRecvBuffer;
  std::vector<int32_t> parRecvCounts;
  std::vector<int32_t> parRDispls;

  // Particles and grid.
  std::vector<Particle> particles;
  Simulation::FlatGrid2D<Cell> grid;

  // Other simulation flags.
  bool isCheckeredTwoProc = false;

  // Precomputed MPI datatypes.
  MPI_Datatype cmType;
  MPI_Datatype particleType;

  // MPI adjacent cells: mapping from ghost directions.
  std::unordered_map<GhostDirection, int32_t> neighbors;
  const std::unordered_map<GhostDirection, std::pair<int32_t, int32_t>>
      directions = {
          {GhostDirection::N, {-1, 0}}, {GhostDirection::NE, {-1, +1}},
          {GhostDirection::E, {0, +1}}, {GhostDirection::SE, {+1, +1}},
          {GhostDirection::S, {+1, 0}}, {GhostDirection::SW, {+1, -1}},
          {GhostDirection::W, {0, -1}}, {GhostDirection::NW, {-1, -1}}};

  subdomain_t sub; ///< Subdomain information.

  // Random seed for particle initialization.
  uint32_t seed;

  // MPI requests.
  MPI_Request requestCOM = MPI_REQUEST_NULL, requestPar = MPI_REQUEST_NULL;

  // Invert cardinals.
  const std::unordered_map<GhostDirection, GhostDirection> complementMap = {
      {N, S}, {S, N}, {E, W}, {W, E}, {NE, SW}, {NW, SE}, {SE, NW}, {SW, NE}};

  //--------------------------------------------------------------------------
  // Simulation Core Computations (Grid, Forces, and Collisions)
  //--------------------------------------------------------------------------

  /**
   * @brief Computes the center of mass for each cell.
   *
   * Iterates over a contiguous range of cells and calculates the cumulative
   * center of mass for particles in each cell.
   */
  void computeCenterOfMass();

  /**
   * @brief Computes the intra-cell gravitational interactions.
   *
   * Calculates pairwise gravitational forces between particles within the same
   * cell.
   *
   * @param particles Reference to the vector of particle pointers in the cell.
   */
  void computeIncellPull(auto &particles);

  /**
   * @brief Computes the gravitational pull from a neighboring cell.
   *
   * Applies the extracell gravitational force, based on the neighbor cell's
   * center-of-mass, to each particle in the provided vector.
   *
   * @param particles Reference to the vector of particle pointers in the source
   * cell.
   * @param cm Center-of-mass of the neighboring cell.
   * @param adjRow Row offset of the adjacent cell.
   * @param adjCol Column offset of the adjacent cell.
   * @param side Domain side length.
   */
  void computeExtracellPull(std::vector<Particle *> &particles,
                            const CenterOfMass &cm, int32_t adjRow,
                            int32_t adjCol, double side);

  /**
   * @brief Computes the gravitational forces between cells.
   *
   * For each cell in a linear range, this function computes both intra-cell and
   * extracell gravitational interactions with neighboring cells.
   */
  void computeGravitationalPull();

  /**
   * @brief Updates particle positions and velocities.
   *
   * Processes each cell in the specified range, updating particle positions and
   * velocities. Particles moving out of their current cell are queued for
   * insertion into the destination cell.
   */
  void computePositionVelocity();

  /**
   * @brief Detects and resolves particle collisions.
   *
   * Scans cells in the given range to identify collisions using a
   * sweep-and-prune approach. When collisions are found, particles are grouped
   * via BFS and the collision count is incremented.
   *
   * @param visited Set to track visited particles during collision grouping.
   * @param q Queue used for BFS grouping of colliding particles.
   * @return Total number of collision groups processed.
   */
  uint64_t computeCollisions(std::unordered_set<Particle *> &visited,
                             std::queue<Particle *> &q);

  /**
   * @brief Recomputes the grid state after particle movements.
   *
   * Processes incoming particles for each cell, resets forces, and clears
   * collision information in preparation for the next timestep.
   */
  void recomputeGrid();

  //--------------------------------------------------------------------------
  // Particle and Grid Initialization
  //--------------------------------------------------------------------------

  /**
   * @brief Computes the simulation subdomain for the current MPI process.
   */
  void computeSubdomain();

  /**
   * @brief Initializes the random seed.
   *
   * Sets the global seed variable by offsetting the provided input seed with a
   * constant.
   *
   * @param inputSeed The input seed value.
   */
  inline void initR4Uni(int32_t inputSeed) { seed = inputSeed + 987654321; }

  /**
   * @brief Generates a uniformly distributed random number in the range [0,1).
   *
   * @return A random double in the range [0,1).
   */
  double rndUniform01();

  /**
   * @brief Generates a normally distributed random number in the range [0,1).
   *
   * Uses the Box-Muller transform to produce a normally distributed value.
   * Repeats until the value is within the desired range.
   *
   * @return A random double in the range [0,1).
   */
  double rndNormal01();

  /**
   * @brief Initializes particle properties (position, velocity, mass, etc.).
   *
   * Depending on the sign of the seed parameter, particles are initialized
   * using either a uniform or a normal distribution. Positions are scaled
   * by the domain side length, velocities are adjusted by cell size, and masses
   * are computed from simulation parameters.
   */
  void initParticles();

  //--------------------------------------------------------------------------
  // MPI Communication and Data Exchange
  //--------------------------------------------------------------------------

  /**
   * @brief Exchanges ghost cell data between neighboring MPI processes.
   */
  void exchange_ghost_cells();

  /**
   * @brief Exchanges particles that have moved between grid cells.
   */
  void exchange_moving_particles();

  /**
   * @brief Exchanges center-of-mass data between MPI processes.
   */
  void exchangeCOM();

  /**
   * @brief Exchanges particle data among MPI processes.
   */
  void exchangePar();

  /**
   * @brief Gathers the first particle's data (particle0) from all processes.
   */
  void gatherParticle0();

  //--------------------------------------------------------------------------
  // Auxiliary and Debug/Utility Functions
  //--------------------------------------------------------------------------

  /**
   * @brief Prints the send buffer information for debugging.
   */
  void debug_print_send() const;

  /**
   * @brief Prints the receive buffer information for debugging.
   */
  void debug_print_recv() const;

  /**
   * @brief Prints the MPI topology for debugging purposes.
   */
  void debug_print_topology() const;

  /**
   * @brief Prints the received data along with neighbor information.
   */
  void debug_print_recv_with_neighbors() const;

  /**
   * @brief Extracts and updates ghost regions from received data.
   */
  void extractAndUpdateGhostRegions();

  /**
   * @brief Updates the ghost region data in a given direction.
   *
   * @param dir The ghost direction.
   * @param data Pointer to the center-of-mass data.
   * @param count Number of data elements.
   */
  void updateGhostRegion(GhostDirection dir, const CenterOfMass *data,
                         int32_t count);

  /**
   * @brief Reduces the local collision count across MPI processes.
   *
   * @param local_collision The collision count from the local process.
   * @return The reduced (global) collision count.
   */
  int32_t reduceCollisionCount(int32_t local_collision);

  /**
   * @brief Performs final clean-up and finalization of the simulation.
   */
  void finalize();

  /**
   * @brief Determines the owner MPI rank for a cell at given grid coordinates.
   *
   * @param newRow Row index of the cell.
   * @param newCol Column index of the cell.
   * @return MPI rank of the owner process.
   */
  int32_t getOwnerRank(int32_t newRow, int32_t newCol);
};

#include "flatGrid2D.tpp"

} // namespace Simulation
#endif // PARSIM_HPP
