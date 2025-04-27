/**
 * @file      macros.hpp
 * @brief     Centralised compile‑time constants for the N‑body simulation.
 *
 * This header collects all numerical parameters that must remain identical
 * across every translation unit.  All constants are expressed in SI units
 * unless explicitly stated.  Keeping them in one place helps avoid the subtle
 * bugs that arise when different compilation units use slightly different
 * values.
 *
 * @note  No symbols from this header produce linker artefacts; everything is
 *        defined as a pre‑processor macro so the compiler can fully fold the
 *        computations at compile‑time.
 */

#define MACROS_HPP
#ifdef MACROS_HPP

// ---------------------------------------------------------------------------
// Standard library includes
// ---------------------------------------------------------------------------
// `_USE_MATH_DEFINES` must be defined *before* <cmath> on MSVC in order to
// expose mathematical constants such as `M_PI`.  Other tool‑chains simply
// ignore it, so defining it unconditionally is harmless.
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
/// @brief Gravitational constant **G**
///  \f$[\text{m}^3\;\text{kg}^{-1}\;\text{s}^{-2}]\f$.
///
/// Used in Newton’s law of universal gravitation:
/// \f$F = G\,m_1 m_2 / r^2\f$.
#define G 6.67408e-11

// ---------------------------------------------------------------------------
// Softening parameter
// ---------------------------------------------------------------------------
/// @brief Softening length \f$\varepsilon\f$ used to avoid singularities
///        when two particles approach each other too closely.
#define EPSILON 0.005

/// @brief The square of the softening length, pre‑computed for efficiency.
#define EPSILON2 (EPSILON * EPSILON)

// ---------------------------------------------------------------------------
// Integration parameters
// ---------------------------------------------------------------------------
/// @brief Fixed time‑step size \f$\Delta t\f$ (seconds) for the integrator.
#define DELTAT 0.1

// ---------------------------------------------------------------------------
// I/O formatting
// ---------------------------------------------------------------------------
/// @brief Controls the number of digits printed after the decimal separator
///        when streaming floating‑point values.
#define DECIMAL_CASES 6

// ---------------------------------------------------------------------------
// Miscellaneous tags
// ---------------------------------------------------------------------------
/// @brief Compile‑time tag identifying the first partition of the simulation
///        domain (useful in MPI/parallel contexts).
#define TAG_PART0 1488

#endif // MACROS_HPP
