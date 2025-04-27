#ifndef PARSIM_HPP
#define PARSIM_HPP

#include <cstdlib>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "parsimUtils.hpp"
#include "pointCalculations.hpp"

namespace Simulation {

/**
 * @brief Represents the simulation engine for the Parsim project.
 *
 * This class manages particles on a grid by handling their initialization,
 * gravitational interactions, movement updates, and collision detection.
 * It uses parallel processing (with OpenMP barriers) to synchronize operations
 * across simulation timesteps.
 */
class Parsim {
public:
  /**
   * @brief Constructs a Parsim simulation instance.
   *
   * @param seed     Random seed for generating initial conditions.
   * @param side     Length of one side of the simulation domain.
   * @param ncside   Number of grid cells along one dimension.
   * @param nPart   Total number of particles to simulate.
   * @param tsteps   Total number of timesteps for the simulation.
   */
  Parsim(int32_t seed, double side, int32_t ncside, uint64_t nPart,
         int32_t tsteps);

  /**
   * @brief Destroys the Parsim simulation instance and frees resources.
   */
  ~Parsim();

  /**
   * @brief Runs the simulation over the specified number of timesteps.
   *
   * For each timestep, the simulation performs these phases:
   *   - Compute the center of mass for each grid cell.
   *   - Calculate intra-cell and inter-cell gravitational forces.
   *   - Update particle positions and velocities.
   *   - Recompute the grid after particles move.
   *   - Detect and group particle collisions using a BFS-based approach.
   *
   * OpenMP barriers ensure that threads are synchronized between these phases.
   *
   * @param visited  A set to track particles already processed during collision
   * grouping.
   * @param q        A queue used for grouping colliding particles via BFS.
   * @param idxStart The starting index for this thread's work.
   * @param idxEnd   The ending (exclusive) index for this thread's work.
   * @return The total number of collision events detected during the
   * simulation.
   */
  uint64_t simulate(std::unordered_set<Particle *> &visited,
                    std::queue<Particle *> &q, const int32_t &idxStart,
                    const int32_t &idxEnd);

  /**
   * @brief Allocates memory for particles and grid cells.
   *
   * This function creates an array of particles and a contiguous block of
   * cells. Each cell is constructed in place and its lock is initialized to
   * ensure safe concurrent access.
   */
  void allocateGrid();

  /**
   * @brief Assigns particles to grid cells based on their positions.
   *
   * For every particle, this function:
   *   - Computes the target cell based on its (x, y) position.
   *   - Converts the 2D cell coordinates into a linear index.
   *   - Inserts the particle into the corresponding cell.
   */
  void populateGrid();

  /**
   * @brief Refreshes grid data after particle movements.
   *
   * Processes each grid cell in the specified range by:
   *   - Integrating incoming particles.
   *   - Resetting forces acting on the particles.
   *   - Clearing previous collision information.
   *
   * @param idxStart The starting cell index (inclusive).
   * @param idxEnd   The ending cell index (exclusive).
   */
  void recomputeGrid(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Retrieves a copy of a particle using its index.
   *
   * @param index The linear index of the particle.
   * @return A copy of the particle at the specified index.
   */
  Particle get_particle(uint64_t index) { return particles[index]; }

private:
  const int32_t seed;   ///< Seed for random number generation.
  const double side;    ///< Length of one side of the simulation domain.
  const int32_t ncside; ///< Number of grid cells per row/column.
  const uint64_t nPart; ///< Total number of particles.
  const int32_t tsteps; ///< Total timesteps for the simulation.

  std::vector<Particle>
      particles;          ///< Container holding all simulation particles.
  std::vector<Cell> grid; ///< Flat container for grid cells.

  /**
   * @brief Computes the center of mass for a range of grid cells.
   *
   * Aggregates mass and position information from particles within each cell to
   * calculate the cell's center of mass.
   *
   * @param idxStart The starting cell index.
   * @param idxEnd   The ending cell index.
   */
  void computeCenterOfMass(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Calculates gravitational forces within a single cell.
   *
   * Evaluates pairwise gravitational interactions among particles residing in
   * the same cell.
   *
   * @param particles Collection of particle pointers within the cell.
   */
  void computeIncellPull(auto &particles);

  /**
   * @brief Applies gravitational influence from a neighboring cell.
   *
   * Uses the center of mass from an adjacent cell to compute and apply
   * gravitational pull on particles in the source cell. This function supports
   * periodic boundary conditions with wrapping in the x and/or y directions.
   *
   * @param particles Collection of particle pointers in the source cell.
   * @param cm        Center of mass of the neighboring cell.
   * @param wrapX    True if x-direction wrapping is applied.
   * @param wrapY    True if y-direction wrapping is applied.
   */
  void computeExtracellPull(auto &particles, CenterOfMass cm, bool wrapX,
                            bool wrapY);

  /**
   * @brief Computes gravitational forces between adjacent cells.
   *
   * For each cell within the specified range, this function calculates:
   *   - Intra-cell gravitational interactions.
   *   - Gravitational pulls from each of the eight neighboring cells.
   *
   * @param idxStart The starting cell index.
   * @param idxEnd   The ending cell index.
   */
  void computeGravitationalPull(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Updates particle positions and velocities.
   *
   * Processes each particle in the given cell range, updating its position and
   * velocity based on the computed forces. If a particle crosses cell
   * boundaries, it is safely re-assigned to the appropriate cell.
   *
   * @param idxStart The starting cell index.
   * @param idxEnd   The ending cell index.
   */
  void computePositionVelocity(const int32_t &idxStart, const int32_t &idxEnd);

  /**
   * @brief Detects and groups colliding particles.
   *
   * Uses a sweep-and-prune method to identify collisions within cells. When
   * collisions occur, a breadth-first search (BFS) groups the colliding
   * particles, and the total collision count is updated.
   *
   * @param visited Set to record particles already processed during collision
   * grouping.
   * @param q       Queue used for BFS to group colliding particles.
   * @param idxStart The starting cell index.
   * @param idxEnd   The ending cell index.
   * @return The number of distinct collision groups detected.
   */
  uint64_t computeCollisions(std::unordered_set<Particle *> &visited,
                             std::queue<Particle *> &q, const int32_t idxStart,
                             const int32_t idxEnd);
};

} // namespace Simulation

#endif // PARSIM_HPP
