#ifndef DS_HPP
#define DS_HPP

#include <algorithm>
#include <omp.h>
#include <queue>
#include <unordered_set>
#include <vector>

#include "macros.hpp"

//------------------------------------------------------------------------------
// GhostDirection Enumeration
//------------------------------------------------------------------------------
// Defines the eight possible directions (N, NE, E, SE, S, SW, W, NW)
// which can be used to identify neighboring cells or ghost cells in the grid.
enum GhostDirection { N = 0, NE, E, SE, S, SW, W, NW };

//------------------------------------------------------------------------------
// Subdomain Structure
//------------------------------------------------------------------------------
// Represents a rectangular subdomain within the simulation grid.
// Stores the starting and ending indices for rows and columns, as well as
// the number of local rows and columns in the subdomain.
struct subdomain_t {
  int32_t rowStart, rowEnd;
  int32_t colStart, colEnd;
  int32_t localRows, localCols;
};

//------------------------------------------------------------------------------
// Particle Structure
//------------------------------------------------------------------------------
// Represents an individual particle within the simulation. This structure
// holds physical properties like position, velocity, mass, and force, as well
// as a list of pointers to other particles that are colliding with it.
struct Particle {
  uint64_t idx;  ///< Unique identifier for the particle.
  double x, y;   ///< Position coordinates.
  double vx, vy; ///< Velocity components.
  double m;      ///< Mass of the particle.
  double fx, fy; ///< Force accumulators for the x and y directions.
  std::vector<Particle *> collidingParticles; ///< List of particles currently
                                              ///< colliding with this one.

  //--------------------------------------------------------------------------
  // Constructors
  //--------------------------------------------------------------------------
  // Default constructor initializes the particle with default values.
  Particle()
      : idx(1488), x(0.0), y(0.0), vx(0.0), vy(0.0), m(0.0), fx(0.0), fy(0.0) {}

  // Parameterized constructor allows initializing all properties.
  Particle(uint64_t idx, double x, double y, double vx, double vy, double m,
           double fx, double fy)
      : idx(idx), x(x), y(y), vx(vx), vy(vy), m(m), fx(fx), fy(fy) {
    collidingParticles.reserve(
        3); // Reserve initial space for collision tracking.
  }

  //--------------------------------------------------------------------------
  // Particle Reset Methods
  //--------------------------------------------------------------------------
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
   * @brief Resets the force accumulators to zero.
   *
   * Should be called at the beginning of each simulation timestep.
   */
  inline void resetForce() {
    fx = 0.0;
    fy = 0.0;
  }

  //--------------------------------------------------------------------------
  // Particle Update Methods
  //--------------------------------------------------------------------------
  /**
   * @brief Adds incremental force contributions to the particle.
   *
   * @param dfx Increment in the x-direction.
   * @param dfy Increment in the y-direction.
   */
  inline void addForce(double dfx, double dfy) {
    fx += dfx;
    fy += dfy;
  }

  /**
   * @brief Updates the particle's position using a basic kinematic equation.
   *
   * The new position is calculated based on current velocity, force-derived
   * acceleration, and a timestep (DELTAT). After updating, the position is
   * wrapped to simulate periodic boundary conditions.
   *
   * @param side The length of one side of the square simulation domain.
   */
  inline void computePos(double side) {
    x = x + vx * DELTAT + 0.5 * (fx / m) * DELTAT * DELTAT;
    y = y + vy * DELTAT + 0.5 * (fy / m) * DELTAT * DELTAT;
    wrapPosition(side);
  }

  /**
   * @brief Wraps the particle's position around the simulation boundaries.
   *
   * Ensures that if a particle moves beyond the simulation domain, it appears
   * on the opposite side (periodic boundary conditions).
   *
   * @param side The length of one side of the square simulation domain.
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
   * Uses the relation: velocity += (force / mass) * DELTAT.
   */
  inline void computeVelocity() {
    vx = vx + (fx / m) * DELTAT;
    vy = vy + (fy / m) * DELTAT;
  }

  //--------------------------------------------------------------------------
  // Collision Management Methods
  //--------------------------------------------------------------------------
  /**
   * @brief Adds a colliding particle to the list.
   *
   * @param par Pointer to the particle that is colliding.
   */
  inline void addColliding(Particle *par) { collidingParticles.push_back(par); }

  /**
   * @brief Clears the list of colliding particles.
   */
  inline void clearColliding() { collidingParticles.clear(); }
};

//------------------------------------------------------------------------------
// CenterOfMass Structure
//------------------------------------------------------------------------------
// Represents the cumulative center of mass for a collection of particles.
// This structure accumulates the total mass and weighted positions (x, y) and
// can be normalized to compute the center of mass.
struct CenterOfMass {
  uint64_t idx;     ///< Identifier for debugging or reference.
  double m;         ///< Cumulative mass.
  double x, y;      ///< Accumulated weighted position coordinates.
  GhostDirection d; ///< Direction (e.g., for ghost cell identification).

  // Default constructor initializes values to zero.
  CenterOfMass() : idx(0), m(0), x(0.0), y(0.0), d(GhostDirection::N) {}

  // Parameterized constructor.
  CenterOfMass(uint64_t idx, double m, double x, double y, GhostDirection d)
      : idx(idx), m(m), x(x), y(y), d(d) {}

  /**
   * @brief Adds additional mass to the cumulative total.
   *
   * @param dm Mass to add.
   */
  inline void computeMass(double dm) { m += dm; }

  /**
   * @brief Adds weighted position contributions.
   *
   * @param dx Contribution in the x-direction.
   * @param dy Contribution in the y-direction.
   */
  inline void computePos(double dx, double dy) {
    x += dx;
    y += dy;
  }

  /**
   * @brief Normalizes the accumulated position by dividing by the total mass.
   *
   * If no mass has been accumulated, resets the position to the origin.
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
   * @brief Resets the cumulative mass and position to initial values.
   */
  inline void reset() {
    m = 0.0;
    x = 0.0;
    y = 0.0;
  }
};

//------------------------------------------------------------------------------
// Cell Structure
//------------------------------------------------------------------------------
// Represents a single cell in the simulation grid. Each cell maintains a
// list of pointers to particles currently within it, an incoming queue for
// particles that have moved into the cell, and its center of mass. An OpenMP
// lock is also included to ensure thread-safe updates.
struct Cell {
  std::vector<Particle *> particles; ///< Pointers to particles in the cell.
  std::vector<Particle *> queue; ///< Queue for particles arriving in the cell.
  CenterOfMass cm; ///< Center of mass for the particles in the cell.
  omp_lock_t lock; ///< Lock for synchronizing cell updates.

  /**
   * @brief Default constructor.
   *
   * Initializes the OpenMP lock for thread safety.
   */
  Cell() { omp_init_lock(&lock); }

  /**
   * @brief Destructor.
   *
   * Destroys the OpenMP lock.
   */
  ~Cell() { omp_destroy_lock(&lock); }

  /**
   * @brief Constructor with reserved capacity.
   *
   * Preallocates memory for the particles and incoming queue vectors to improve
   * performance when the expected number of particles is known.
   *
   * @param num_particles Expected number of particles in the cell.
   */
  Cell(uint64_t num_particles) {
    particles.reserve(num_particles);
    queue.reserve(num_particles);
  }

  /**
   * @brief Adds an incoming particle to the cell's queue.
   *
   * This function is used to queue particles that have moved into the cell
   * and will be processed in a later phase.
   *
   * @param p Pointer to the incoming particle.
   */
  inline void addIncoming(Particle *p) { queue.push_back(p); }

  /**
   * @brief Processes the incoming queue.
   *
   * Transfers all particles from the incoming queue into the main particles
   * vector and clears the queue afterwards.
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
   * Bypasses the incoming queue and directly inserts the particle pointer into
   * the particles vector.
   *
   * @param particle Pointer to the particle to be added.
   */
  inline void addParticle(Particle *particle) { particles.push_back(particle); }

  /**
   * @brief Removes a particle from the cell.
   *
   * Efficiently removes the particle at the specified index by swapping it
   * with the last element and then removing the last element.
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

#endif // DS_HPP
