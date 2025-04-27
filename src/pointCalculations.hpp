#ifndef POINTCALCULATIONS_HPP
#define POINTCALCULATIONS_HPP

#include <algorithm>
#include <cmath>

/**
 * @brief Wraps the index around the grid boundaries.
 *
 * If the index is outside the range [0, ncside), this function adjusts it using
 * modular arithmetic so that it falls within the valid range.
 *
 * @param idx    The index to be wrapped. It is updated in-place.
 * @param ncside The number of cells per side of the grid.
 */
inline void computeWrap(int32_t &idx, int32_t ncside) {
  if (idx < 0 || idx >= ncside) {
    idx = ((idx % ncside) + ncside) % ncside;
  }
}

/**
 * @brief Computes the grid cell index for a given position.
 *
 * Given a position in one dimension, this function calculates the corresponding
 * cell index by dividing the position by the cell side length. The result is
 * then wrapped to ensure it falls within the grid boundaries.
 *
 * @param pos     The position (either x or y) in the simulation space.
 * @param side    The total length of the simulation space along that dimension.
 * @param ncside  The number of cells per side of the grid.
 * @return int32_t The index of the cell corresponding to the position.
 */
inline int32_t getIndex(double pos, double side, int32_t ncside) {
  double cellSide = side / ncside;
  int32_t idx = static_cast<int32_t>(pos / cellSide);
  computeWrap(idx, ncside);
  return idx;
}

#endif // POINTCALCULATIONS_HPP
