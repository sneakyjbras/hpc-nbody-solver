#ifndef DS_HPP
#define DS_HPP

#include <algorithm>
#include <omp.h>
#include <queue>
#include <unordered_set>
#include <vector>

#include "macros.hpp"

namespace Simulation {
/**
 * @brief Represents a particle in the simulation.
 *
 * Contains properties such as position, velocity, mass, and force, along with a
 * list of colliding particles. Provides inline methods to reset its state,
 * update its motion, and manage collisions.
 */
struct Particle {
  double x, y;   ///< Position of the particle.
  double vx, vy; ///< Velocity components.
  double m;      ///< Mass of the particle.
  double fx, fy; ///< Accumulated force components.
  std::vector<Particle *>
      collidingParticles; ///< List of pointers to particles that are colliding
                          ///< with this one.
  /**
   * @brief Resets the particle's position to the origin.
   */
  inline void resetPosition() {
    x = 0.0;
    y = 0.0;
  }

  /**
   * @brief Resets the particle's mass to zero.
   */
  inline void resetMass() { m = 0.0; }

  /**
   * @brief Resets the particle's velocity to zero.
   */
  inline void resetVelocity() {
    vx = 0.0;
    vy = 0.0;
  }

  /**
   * @brief Resets the force accumulators for a new timestep.
   */
  inline void resetForce() {
    fx = 0.0;
    fy = 0.0;
  }

  /**
   * @brief Adds a force contribution to the particle.
   *
   * @param dfx Additional force in the x-direction.
   * @param dfy Additional force in the y-direction.
   */
  inline void addForce(double dfx, double dfy) {
    fx += dfx;
    fy += dfy;
  }

  /**
   * @brief Updates the particle's position based on its velocity and
   * acceleration.
   *
   * The position is updated using a basic kinematic equation that accounts for
   * the current velocity, applied force, and timestep DELTAT. The position is
   * then wrapped around the simulation space.
   *
   * @param side The side length of the simulation area.
   */
  inline void computePos(double side) {
    x = x + vx * DELTAT + 0.5 * (fx / m) * DELTAT * DELTAT;
    y = y + vy * DELTAT + 0.5 * (fy / m) * DELTAT * DELTAT;
    wrapPosition(side);
  }

  /**
   * @brief Wraps the particle's position around the simulation boundaries.
   *
   * If the particle moves out of the bounds of the simulation space, its
   * position is wrapped around to the opposite side, simulating periodic
   * boundary conditions.
   *
   * @param side The side length of the simulation area.
   */
  inline void wrapPosition(double side) {
    if (x < 0.0)
      x += side;
    else if (x >= side)
      x -= side;

    if (y < 0.0)
      y += side;
    else if (y >= side)
      y -= side;
  }

  /**
   * @brief Updates the particle's velocity based on the applied force.
   *
   * Uses the accumulated force and the particle's mass to compute the new
   * velocity, accounting for the timestep DELTAT.
   */
  inline void computeVelocity() {
    vx = vx + (fx / m) * DELTAT;
    vy = vy + (fy / m) * DELTAT;
  }

  /**
   * @brief Adds a particle to the colliding particles list.
   *
   * @param par Pointer to the colliding particle.
   */
  inline void addColliding(Particle *par) { collidingParticles.push_back(par); }

  /**
   * @brief Clears the list of colliding particles.
   */
  inline void clearColliding() { collidingParticles.clear(); }
};

/**
 * @brief Represents the center of mass for a group of particles.
 *
 * Contains the cumulative mass and weighted position (x, y) of the particles,
 * along with an index for debugging purposes. Provides inline methods to update
 * and reset the center of mass.
 */
struct CenterOfMass {
  double m;    ///< Total mass of the particles.
  double x, y; ///< Position of the center of mass.

  /**
   * @brief Accumulates mass into the center of mass.
   *
   * @param dm Additional mass to be added.
   */
  inline void computeMass(double dm) { m += dm; }

  /**
   * @brief Accumulates weighted position values into the center of mass.
   *
   * @param dx Weighted x-coordinate contribution.
   * @param dy Weighted y-coordinate contribution.
   */
  inline void computePos(double dx, double dy) {
    x += dx;
    y += dy;
  }

  /**
   * @brief Normalizes the center of mass position.
   *
   * Divides the accumulated position by the total mass, provided the mass is
   * non-zero. If the mass is zero, the position is reset to the origin.
   */
  inline void normalizePos() {
    if (m > 0) {
      x /= m;
      y /= m;
    } else {
      x = 0.0;
      y = 0.0;
    }
  }

  /**
   * @brief Resets the center of mass to its initial state.
   */
  inline void reset() {
    m = 0.0;
    x = 0.0;
    y = 0.0;
  }
};

/**
 * @brief Represents a cell in the simulation grid.
 *
 * A cell holds pointers to particles currently within its bounds as well as a
 * queue for incoming particles that have moved into the cell. It also maintains
 * the center of mass for the particles within the cell.
 */
struct Cell {
  std::vector<Particle *>
      particles; ///< Pointers to particles currently in the cell.
  std::vector<Particle *>
      queue;              ///< Pointers to particles queued for addition.
                          ///< //							  /
  struct CenterOfMass cm; ///< Center of mass of the particles in the cell.

  Cell() {}

  /**
   * @brief Constructs a cell with reserved capacity for particles.
   *
   * Reserves memory for both the particles vector and the incoming queue based
   * on the expected number of particles.
   *
   * @param numParticles Anticipated number of particles in the cell.
   */
  Cell(uint64_t numParticles) {
    particles.reserve(numParticles);
    queue.reserve(numParticles);
  }

  /**
   * @brief Adds an incoming particle to the cell's queue.
   *
   * Particles that move into the cell are first added to a queue and then later
   * processed.
   *
   * @param p Pointer to the particle to be added.
   */
  inline void addIncoming(Particle *p) { queue.push_back(p); }

  /**
   * @brief Processes the incoming queue.
   *
   * Transfers all particles from the incoming queue to the main particles
   * vector and then clears the queue.
   */
  void processIncoming() {
    for (uint64_t i = 0; i < queue.size(); i++) {
      particles.push_back(queue[i]);
    }
    queue.clear();
  }

  /**
   * @brief Adds a particle directly to the cell.
   *
   * Directly inserts a particle pointer into the particles vector, bypassing
   * the incoming queue.
   *
   * @param particle Pointer to the particle to be added.
   */
  inline void addParticle(Particle *particle) { particles.push_back(particle); }

  /**
   * @brief Removes a particle from the cell.
   *
   * Efficiently removes the particle at the given index by swapping it with the
   * last particle and then removing the last element.
   *
   * @param index Index of the particle to remove.
   */
  void removeParticle(uint64_t index) {
    if (index < particles.size()) {
      std::swap(particles[index], particles.back());
      particles.pop_back();
    }
  }
};
} // namespace Simulation
#endif // DS_HPP
