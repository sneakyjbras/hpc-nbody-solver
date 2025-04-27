#ifndef POINTCALCULATIONS_HPP
#define POINTCALCULATIONS_HPP

#include <algorithm>
#include <cmath>

/**
 * @brief Wraps an index so that it lies within the valid grid boundaries.
 *
 * This function adjusts the given index using modular arithmetic. If the index
 * is less than 0 or greater than or equal to the number of cells per side
 * (ncside), it wraps the index into the range [0, ncside). The index is updated
 * in-place.
 *
 * @param idx    The index to be wrapped. This parameter is modified to contain
 *               the wrapped index.
 * @param ncside The number of cells along one side of the grid.
 */
inline void computeWrap(int32_t &idx, int32_t ncside) {
  if (idx < 0 || idx >= ncside) {
    idx = ((idx % ncside) + ncside) % ncside;
  }
}

/**
 * @brief Calculates the grid cell index for a given position.
 *
 * This function computes the cell index corresponding to a position along one
 * dimension. It divides the position by the length of a single cell to
 * determine the index. The resulting index is then wrapped into the valid range
 * using computeWrap().
 *
 * @param pos    The coordinate (either x or y) within the simulation space.
 * @param side   The total length of the simulation space along the given
 * dimension.
 * @param ncside The number of cells per side of the grid.
 * @return int32_t The index of the grid cell that contains the position.
 */
inline int32_t getIndex(double pos, double side, int32_t ncside) {
  double cellSide = side / ncside;
  int32_t idx = static_cast<int32_t>(pos / cellSide);
  computeWrap(idx, ncside);
  return idx;
}

#endif // POINTCALCULATIONS_HPP
