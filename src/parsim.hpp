#ifndef PARSIM_HPP
#define PARSIM_HPP

#include <cstdlib>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>

#include "parsimUtils.hpp"
#include "pointCalculations.hpp"

namespace Simulation {

class Parsim {
public:
  // Constructor & Destructor
  Parsim(int32_t seed, double side, int32_t ncside, uint64_t nPart,
         int32_t tsteps);
  ~Parsim();

  /**
   * @brief Runs the simulation for a specified number of timesteps.
   *
   * Iterates over the simulation timesteps, calling the necessary functions
   * (center-of-mass computation, gravitational pull, particle update, grid
   * recomputation, and collision detection) in sequence. OpenMP barriers
   * synchronize threads between phases.
   *
   * @param grid Flat grid array.
   * @param side Domain side length.
   * @param ncside Grid side length.
   * @param tsteps Number of timesteps to simulate.
   * @param visited Set to track visited particles for collision grouping.
   * @param q Queue used for BFS grouping.
   * @param idxStart Starting linear index for this thread's work.
   * @param idxEnd Ending linear index for this thread's work.
   * @return Total collisions detected over all timesteps.
   */
  uint64_t simulate(std::unordered_set<Particle *> &visited,
                    std::queue<Particle *> &q, const int32_t &idxStart,
                    const int32_t &idxEnd);

  /**
   * @brief Allocates the particle array and a contiguous (flat) grid.
   *
   * Allocates an array of particles and a single block of (ncside * ncside)
   * cells, constructing each cell in-place and initializing its lock.
   *
   * @param particles Reference to pointer to the Particle array.
   * @param nPart Number of particles.
   * @param grid Reference to pointer to the flat Cell array.
   * @param ncside Grid side length.
   */
  void allocateGrid();
  /**
   * @brief Populates the flat grid with particles.
   *
   * For each particle, computes its target cell based on its (x,y) position,
   * converts the 2D (row, col) into a linear index, and inserts the particle
   * into that cell.
   *
   * @param par Particle array.
   * @param nPart Number of particles.
   * @param grid Flat grid array.
   * @param side Simulation domain side length.
   * @param ncside Grid side length.
   */
  void populateGrid();

  /**
   * @brief Recomputes the grid state after particle movements.
   *
   * Processes incoming particles for each cell in the specified linear index
   * range, then resets forces and clears collision information for all
   * particles in the cell.
   *
   * @param grid Flat grid array.
   * @param idxStart Starting linear index (inclusive).
   * @param idxEnd Ending linear index (exclusive).
   */
  void recomputeGrid(const int32_t &idxStart, const int32_t &idxEnd);

  Particle get_particle(uint64_t index) { return particles[index]; }

private:
  const int32_t seed;
  const double side;
  const int32_t ncside;
  const uint64_t nPart;
  const int32_t tsteps;

  std::vector<Particle> particles;
  std::vector<Cell> grid;
  /**
   * @brief Computes the center of mass for each cell.
   *
   * Iterates over the specified linear index range and, for each cell, computes
   * the cumulative center of mass from all particles residing in the cell.
   *
   * @param grid Flat grid array.
   * @param idxStart Starting linear index.
   * @param idxEnd Ending linear index.
   */
  void computeCenterOfMass(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Computes gravitational interactions within a single cell.
   *
   * For a given cell's particle list, computes the intra-cell gravitational
   * pull using a pairwise interaction approach.
   *
   * @param side Domain side length.
   * @param particles Reference to the vector of particle pointers.
   */
  void computeIncellPull(auto &particles);

  /**
   * @brief Computes the gravitational pull from an adjacent cell.
   *
   * Using the center-of-mass of a neighboring cell, this function computes
   * and applies the extracell gravitational force to all particles in the given
   * vector.
   *
   * @param side Domain side length.
   * @param particles Reference to the vector of particle pointers from the
   * source cell.
   * @param cm Center-of-mass of the neighboring cell.
   * @param wrapX Flag indicating if x-direction wrapping was applied.
   * @param wrapY Flag indicating if y-direction wrapping was applied.
   */
  void computeExtracellPull(auto &particles, CenterOfMass cm, bool wrapX,
                            bool wrapY);

  /**
   * @brief Computes gravitational forces between cells.
   *
   * Iterates over a contiguous linear range of cells, computing both the
   * intra-cell interactions and the extracell pull from each of the eight
   * neighbors. The 2D coordinates are derived from the linear index.
   *
   * @param grid Flat grid array.
   * @param side Domain side length.
   * @param ncside Grid side length.
   * @param idxStart Starting linear index.
   * @param idxEnd Ending linear index.
   */
  void computeGravitationalPull(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Updates particle positions and velocities.
   *
   * Processes each cell in the given linear index range, updating particle
   * positions and velocities. If a particle moves to a new cell, it is removed
   * from the current cell and added to the incoming queue of the destination
   * cell (with proper locking).
   *
   * @param grid Flat grid array.
   * @param side Domain side length.
   * @param ncside Grid side length.
   * @param idxStart Starting linear index.
   * @param idxEnd Ending linear index.
   */
  void computePositionVelocity(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Detects and resolves collisions within cells.
   *
   * Scans the particles in each cell (within the given linear index range) to
   * detect collisions using a sweep-and-prune method. When collisions are
   * detected, groups of colliding particles are formed via BFS, and the
   * collision count is incremented.
   *
   * @param grid Flat grid array.
   * @param visited Set to track visited particles during grouping.
   * @param q Queue used for BFS grouping of colliding particles.
   * @param idxStart Starting linear index.
   * @param idxEnd Ending linear index.
   * @return Total number of collision groups processed.
   */
  uint64_t computeCollisions(std::unordered_set<Particle *> &visited,
                             std::queue<Particle *> &q, const int32_t idxStart,
                             const int32_t idxEnd);
};

} // namespace Simulation
#endif // PARSIM_HPP
