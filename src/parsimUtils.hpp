#ifndef PARSIMUTILS_HPP
#define PARSIMUTILS_HPP

#include <iomanip>
#include <iostream>

#include "ds.hpp"

namespace Simulation {
/**
 * @brief Prints the details of a single particle for debugging.
 *
 * Outputs the particle's index, mass, position (x, y), and velocity (vx, vy)
 * using fixed-point notation with a preset precision.
 *
 * @param p     The particle to be debugged.
 * @param i The index of the particle.
 */
void debugParticle(const Particle p, uint64_t i);

/**
 * @brief Iterates over an array of particles and prints their details.
 *
 * Calls debugParticle() for each particle in the provided array.
 *
 * @param particles Array of particles to be debugged.
 * @param nPart    Total number of particles.
 */
void debugParticles(Particle *particles, uint64_t nPart);

/**
 * @brief Prints the center of mass details of each cell in a 2D grid.
 *
 * For every cell, the function outputs its index and the corresponding center
 * of mass components (x, y) and mass using fixed formatting with a preset
 * precision.
 *
 * @param cells   2D grid of cells.
 * @param ncside  Number of cells per side of the grid.
 */
void debugCells(Cell **cells, uint32_t ncside);

/**
 * @brief Prints debugging information for force interactions.
 *
 * Outputs a formatted message containing the labels of interacting entities,
 * the force magnitude, and its x and y components.
 *
 * @param label1 First label for the force interaction.
 * @param label2 Second label for the force interaction.
 * @param dx     Force component in the x-direction.
 * @param dy     Force component in the y-direction.
 * @param mag    Magnitude of the force.
 */
void debugForce(const std::string &label1, const std::string &label2, double dx,
                double dy, double mag);

/**
 * @brief Prints debugging information for a collision event.
 *
 * Outputs a formatted message with the provided labels and the collision
 * distance.
 *
 * @param label1 First label for the collision event.
 * @param label2 Second label for the collision event.
 * @param d      The collision distance.
 */
void debugCollision(const std::string &label1, const std::string &label2,
                    double &d);
} // namespace Simulation
#endif // PARSIMUTILS_HPP
