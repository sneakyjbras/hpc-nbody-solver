#ifndef PARSIM_UTILS_HPP
#define PARSIM_UTILS_HPP

#include <iomanip>
#include <iostream>

#include "ds.hpp"
#include "flatGrid2D.hpp"

namespace Simulation {

/**
 * @brief Prints the details of a single particle for debugging purposes.
 *
 * This function outputs the particle's index, mass, position (x, y), and
 * velocity (vx, vy) using fixed-point notation with a preset precision defined
 * by DECIMAL_CASES.
 *
 * @param p     The particle whose details are to be printed.
 * @param i The index of the particle (for additional reference in the
 * output).
 */
void debugParticle(const Particle p);

/**
 * @brief Iterates over an array of particles and prints their details.
 *
 * Calls debugParticle() for each particle in the provided array to output
 * their properties for debugging.
 *
 * @param particles Pointer to an array of Particle objects.
 * @param nPart     The total number of particles in the array.
 */
void debugParticles(Particle *particles, uint64_t nPart);

/**
 * @brief Prints the center of mass details of each cell in a 2D grid.
 *
 * For every cell in the provided grid, this function outputs its index along
 * with the corresponding center of mass components (x, y) and mass. The output
 * is formatted using fixed-point notation with a preset precision defined by
 * DECIMAL_CASES.
 *
 * @param grid A FlatGrid2D object representing the grid of cells.
 */
void debugCells(FlatGrid2D<Cell> grid);

/**
 * @brief Prints debugging information for the force interaction between a
 * particle and a cell.
 *
 * Outputs a formatted message showing the force magnitude and its x and y
 * components for an interaction between a particle (identified by label1) and a
 * cell (identified by label2).
 *
 * @param label1 Label or identifier for the particle.
 * @param label2 Label or identifier for the cell.
 * @param dx     Force component in the x-direction.
 * @param dy     Force component in the y-direction.
 * @param mag    Magnitude of the force.
 */
void debugCellForce(const uint64_t &label1, const uint64_t &label2, double dx,
                    double dy, double mag);

/**
 * @brief Prints debugging information for the force interaction between two
 * particles.
 *
 * Outputs a formatted message with the labels of the interacting particles
 * along with the force magnitude and its x and y components.
 *
 * @param label1 Label or identifier for the first particle.
 * @param label2 Label or identifier for the second particle.
 * @param dx     Force component in the x-direction.
 * @param dy     Force component in the y-direction.
 * @param mag    Magnitude of the force.
 */
void debugParForce(const uint64_t &label1, const uint64_t &label2, double dx,
                   double dy, double mag);

/**
 * @brief Prints debugging information for a collision event between two
 * particles.
 *
 * Outputs a formatted message including the labels of the colliding particles
 * and the collision distance.
 *
 * @param label1 Label or identifier for the first particle.
 * @param label2 Label or identifier for the second particle.
 * @param d      The collision distance (passed by reference).
 */
void debugCollision(const uint64_t &label1, const uint64_t &label2, double &d);

} // namespace Simulation

#endif // PARSIM_UTILS_HPP
