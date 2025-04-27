#ifndef FLATGRID2D_HPP
#define FLATGRID2D_HPP

#include <cassert>
#include <functional>
#include <iomanip> // For output formatting
#include <iostream>
#include <vector>

namespace Simulation {

/**
 * @brief A flat (1D vector) representation of a 2D grid.
 *
 * This template class provides a two-dimensional grid abstraction backed by a
 * one-dimensional std::vector. It supports element access, grid manipulation,
 * and debugging routines. The grid is defined by local (subdomain) and total
 * dimensions.
 *
 * @tparam T The type of the grid elements.
 */
template <typename T> class FlatGrid2D {
public:
  //--------------------------------------------------------------------------
  // Constructors and Destructor
  //--------------------------------------------------------------------------

  /**
   * @brief Default constructor.
   *
   * Initializes an empty grid.
   */
  FlatGrid2D();

  /**
   * @brief Constructs a grid with specified local dimensions.
   *
   * The total dimensions are computed based on the local dimensions.
   *
   * @param localRows Number of rows in the local grid.
   * @param localCols Number of columns in the local grid.
   */
  FlatGrid2D(int32_t localRows, int32_t localCols);

  /**
   * @brief Destructor.
   *
   * Cleans up resources used by the grid.
   */
  ~FlatGrid2D();

  //--------------------------------------------------------------------------
  // Element Access Methods
  //--------------------------------------------------------------------------

  /**
   * @brief Accesses an element at a given row and column.
   *
   * Performs bounds checking and prints an error message then aborts
   * if the indices are out of range.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the element at (i, j).
   */
  inline T &at(int32_t i, int32_t j) {
    if (!(i >= 0 && i < totalRows)) {
      std::cerr << "Row index out of bounds: i = " << i
                << ", expected in range [0, " << totalRows - 1 << "]\n";
      std::abort();
    }
    if (!(j >= 0 && j < totalCols)) {
      std::cerr << "Column index out of bounds: j = " << j
                << ", expected in range [0, " << totalCols - 1 << "]\n";
      std::abort();
    }
    return data[i * totalCols + j];
  }

  /**
   * @brief Accesses an element at a given row and column (const version).
   *
   * Uses assert to ensure that the indices are within the valid range.
   *
   * @param i Row index.
   * @param j Column index.
   * @return Const reference to the element at (i, j).
   */
  inline const T &at(int32_t i, int32_t j) const {
    assert(i >= 0 && i < totalRows);
    assert(j >= 0 && j < totalCols);
    return data[i * totalCols + j];
  }

  /**
   * @brief Computes the 1D index corresponding to the 2D coordinates.
   *
   * @param i Row index.
   * @param j Column index.
   * @return The 1D index in the underlying vector.
   */
  inline int32_t index(int32_t i, int32_t j) const { return i * totalCols + j; }

  /**
   * @brief Provides direct access to the underlying data array.
   *
   * @return Pointer to the beginning of the data array.
   */
  inline T *rawData() { return data.data(); }

  /**
   * @brief Provides const access to the underlying data array.
   *
   * @return Const pointer to the beginning of the data array.
   */
  inline const T *rawData() const { return data.data(); }

  //--------------------------------------------------------------------------
  // Grid Dimension Accessors
  //--------------------------------------------------------------------------

  /**
   * @brief Retrieves the total number of rows in the grid.
   *
   * @return Total row count.
   */
  inline int32_t getTotalRows() const { return totalRows; }

  /**
   * @brief Retrieves the total number of columns in the grid.
   *
   * @return Total column count.
   */
  inline int32_t getTotalCols() const { return totalCols; }

  /**
   * @brief Retrieves the number of local rows.
   *
   * Local rows might represent the subdomain for parallel computations.
   *
   * @return Number of local rows.
   */
  inline int32_t getLocalRows() const { return localRows; }

  /**
   * @brief Retrieves the number of local columns.
   *
   * Local columns might represent the subdomain for parallel computations.
   *
   * @return Number of local columns.
   */
  inline int32_t getLocalCols() const { return localCols; }

  //--------------------------------------------------------------------------
  // Grid Manipulation Methods
  //--------------------------------------------------------------------------

  /**
   * @brief Fills the grid with a specified value.
   *
   * @param value The value to fill each element of the grid.
   */
  void fill(const T &value);

  //--------------------------------------------------------------------------
  // Debugging Methods
  //--------------------------------------------------------------------------

  /**
   * @brief Prints the grid for debugging purposes.
   *
   * A custom formatter can be provided to convert elements to strings.
   * If no formatter is provided, the default output operator is used.
   *
   * @param formatter Optional function to format grid elements.
   */
  void
  printGrid(std::function<std::string(const T &)> formatter = nullptr) const;

  /**
   * @brief Debug method to print the halo (boundary) of the grid.
   *
   * Useful in distributed simulations to verify halo exchanges.
   *
   * @param worldRank The rank of the current process (for distributed
   * contexts).
   */
  void debugPrintHalo(uint32_t worldRank) const;

  /**
   * @brief Debug method to print the ghost halo of the grid.
   *
   * Helps in debugging the exchange of ghost cells between neighboring
   * processes.
   *
   * @param worldRank The rank of the current process (for distributed
   * contexts).
   */
  void debugGhostHalo(uint32_t worldRank) const;

private:
  int32_t localRows, localCols; ///< Local grid dimensions (subdomain size).
  int32_t totalRows, totalCols; ///< Total grid dimensions.
  std::vector<T> data; ///< Underlying data storage in a flat 1D vector.
};

} // namespace Simulation

#include "flatGrid2D.tpp"

#endif // FLATGRID2D_HPP
