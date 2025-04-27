#ifndef DS_HPP
#define DS_HPP

#include "macros.hpp"
#include <algorithm>
#include <omp.h>
#include <queue>
#include <unordered_set>
#include <vector>

namespace Simulation {

/**
 * @brief Represents a particle within the simulation.
 *
 * This structure encapsulates a particle's state, including its spatial
 * coordinates, velocity, mass, and accumulated force. It also maintains a list
 * of colliding particles. Inline methods are provided for state resetting and
 * motion updates.
 */
struct Particle {
  double x, y;   ///< X and Y coordinates of the particle.
  double vx, vy; ///< Velocity components along the x and y axes.
  double m;      ///< Mass of the particle.
  double fx, fy; ///< Accumulated force components in the x and y directions.

  /// List of pointers to particles colliding with this particle.
  std::vector<Particle *> collidingParticles;

  /**
   * @brief Resets the particle's position to the origin (0,0).
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
   * @brief Resets the particle's velocity components to zero.
   */
  inline void resetVelocity() {
    vx = 0.0;
    vy = 0.0;
  }

  /**
   * @brief Resets the accumulated force components to zero.
   */
  inline void resetForce() {
    fx = 0.0;
    fy = 0.0;
  }

  /**
   * @brief Adds an incremental force to the particle's current force.
   *
   * @param dfx Additional force in the x-direction.
   * @param dfy Additional force in the y-direction.
   */
  inline void addForce(double dfx, double dfy) {
    fx += dfx;
    fy += dfy;
  }

  /**
   * @brief Computes the new position of the particle using basic kinematics.
   *
   * The updated position is calculated using the current velocity, force, mass,
   * and a fixed timestep (DELTAT). The position is then wrapped around the
   * simulation boundaries.
   *
   * @param side The side length of the simulation domain.
   */
  inline void computePos(double side) {
    x = x + vx * DELTAT + 0.5 * (fx / m) * DELTAT * DELTAT;
    y = y + vy * DELTAT + 0.5 * (fy / m) * DELTAT * DELTAT;
    wrapPosition(side);
  }

  /**
   * @brief Applies periodic boundary conditions to the particle's position.
   *
   * If the particle exceeds the simulation boundaries, it is wrapped to the
   * opposite side.
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
   * @brief Updates the particle's velocity based on the accumulated force.
   *
   * The velocity update uses the formula: new velocity = current velocity +
   * (force / mass) * DELTAT.
   */
  inline void computeVelocity() {
    vx = vx + (fx / m) * DELTAT;
    vy = vy + (fy / m) * DELTAT;
  }

  /**
   * @brief Adds a pointer to another particle that is colliding with this one.
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
 * This structure accumulates the mass and weighted positions of particles,
 * allowing the center of mass to be computed and normalized.
 */
struct CenterOfMass {
  double m;    ///< Total mass accumulated.
  double x, y; ///< Weighted sum of x and y coordinates.

  /**
   * @brief Adds an incremental mass to the current total.
   *
   * @param dm Additional mass to add.
   */
  inline void computeMass(double dm) { m += dm; }

  /**
   * @brief Accumulates weighted position contributions.
   *
   * @param dx Weighted x-coordinate contribution.
   * @param dy Weighted y-coordinate contribution.
   */
  inline void computePos(double dx, double dy) {
    x += dx;
    y += dy;
  }

  /**
   * @brief Normalizes the position by dividing the accumulated position by the
   * mass.
   *
   * If the mass is zero, the position is reset to the origin.
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
 * A cell is responsible for storing pointers to particles that reside within
 * its boundaries, managing an incoming queue for particles that move in, and
 * computing its own center of mass. It uses an OpenMP lock for safe concurrent
 * access.
 */
struct Cell {
  std::vector<Particle *>
      particles; ///< Pointers to particles currently in the cell.
  std::vector<Particle *> queue; ///< Queue for incoming particles.
  CenterOfMass cm;               ///< Center of mass for the cell.
  omp_lock_t lock;               ///< OpenMP lock to protect cell data.

  /**
   * @brief Default constructor.
   *
   * Initializes the OpenMP lock for the cell.
   */
  Cell() { omp_init_lock(&lock); }

  /**
   * @brief Destructor.
   *
   * Destroys the OpenMP lock.
   */
  ~Cell() { omp_destroy_lock(&lock); }

  /**
   * @brief Constructs a cell with pre-reserved capacity.
   *
   * Pre-reserves memory for both the particles vector and the incoming queue
   * based on the expected number of particles.
   *
   * @param numParticles Estimated number of particles for the cell.
   */
  Cell(uint64_t numParticles) {
    particles.reserve(numParticles);
    queue.reserve(numParticles);
  }

  /**
   * @brief Adds a particle pointer to the incoming queue.
   *
   * Particles that move into the cell are initially stored in this queue.
   *
   * @param p Pointer to the incoming particle.
   */
  inline void addIncoming(Particle *p) { queue.push_back(p); }

  /**
   * @brief Processes the incoming particle queue.
   *
   * Transfers all particles from the queue to the main particles vector and
   * then clears the queue.
   */
  void processIncoming() {
    for (uint64_t i = 0; i < queue.size(); i++) {
      particles.push_back(queue[i]);
    }
    queue.clear();
  }

  /**
   * @brief Directly adds a particle to the cell.
   *
   * Inserts a particle pointer into the particles vector.
   *
   * @param particle Pointer to the particle to be added.
   */
  inline void addParticle(Particle *particle) { particles.push_back(particle); }

  /**
   * @brief Removes a particle from the cell.
   *
   * Removes the particle at the specified index by swapping it with the last
   * element and then popping the last element.
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
