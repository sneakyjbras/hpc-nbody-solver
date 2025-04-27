#ifndef INITPARTICLES_HPP
#define INITPARTICLES_HPP

#include "ds.hpp"

namespace Simulation {
/**
 * @brief Initializes the random seed.
 *
 * Sets the global seed variable by offsetting the given input seed with a
 * constant.
 *
 * @param inputSeed The input seed value.
 */
void initR4Uni(int32_t inputSeed);

/**
 * @brief Generates a uniformly distributed random number in the range [0,1).
 *
 * Uses bitwise operations on a global seed to produce a pseudo-random value,
 * which is then normalized.
 *
 * @return double A uniformly distributed random number in the range [0,1).
 */
double rndUniform01();

/**
 * @brief Generates a normally distributed random number within the range [0,1).
 *
 * Uses the Box-Muller transform to produce a normally distributed value, then
 * shifts and scales it. If the resulting value is outside the [0,1) range, the
 * process is repeated until a valid value is obtained.
 *
 * @return double A normally distributed random number in the range [0,1).
 */
double rndNormal01();

/**
 * @brief Initializes particle properties such as position, velocity, and mass.
 *
 * Depending on the sign of the seed parameter, particles are initialized using
 * either a uniform or a normal distribution. Positions are scaled by the
 * simulation area's side length, velocities are adjusted based on the grid cell
 * size, and masses are computed from a formula involving the total number of
 * particles, grid dimensions, and physical constants.
 *
 * @param seed      Random seed for initialization. A negative seed value
 * selects the normal distribution.
 * @param side      The side length of the simulation area.
 * @param ncside    Number of cells per side of the grid.
 * @param nPart    Total number of particles.
 * @param par       Array of particles to be initialized.
 */
void initParticles(long seed, double side, long ncside, long long nPart,
                   Particle *par);
} // namespace Simulation
#endif // INITPARTICLES_HPP
