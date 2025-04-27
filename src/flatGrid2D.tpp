#ifndef FLATGRID2D_TPP
#define FLATGRID2D_TPP

namespace Simulation {


template <typename T>
FlatGrid2D<T>::FlatGrid2D()
    : FlatGrid2D(0, 0) {
}

template <typename T>
FlatGrid2D<T>::FlatGrid2D(int32_t localRows, int32_t localCols)
    : localRows(localRows), localCols(localCols), totalRows(localRows + 2),
      totalCols(localCols + 2), data(totalRows * totalCols) {}
template <typename T> FlatGrid2D<T>::~FlatGrid2D(){};

template <typename T> void FlatGrid2D<T>::fill(const T &value) {
  std::fill(data.begin(), data.end(), value);
}

template <typename T>
void FlatGrid2D<T>::printGrid(
    std::function<std::string(const T &)> formatter) const {
  std::cout << "\nFlatGrid2D Debug Print:\n";
  for (int32_t i = 0; i < totalRows; i++) {
    for (int32_t j = 0; j < totalCols; j++) {
      if (formatter) {
        std::cout << std::setw(6) << formatter(at(i, j)) << " ";
      } else {
        std::cout << "[x] ";
      }
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

template <typename T>
void FlatGrid2D<T>::debugPrintHalo(uint32_t worldRank) const {
  if (worldRank != 0)
    return;

  std::cout << "Ghost halo for rank 0:\n";

  // Top row
  std::cout << "Top row: ";
  for (int32_t j = 1; j <= localCols; ++j) {
    std::cout << at(0, j).cm.m << " ";
  }
  std::cout << "\n";

  // Left and right columns (excluding corners already printed)
  for (int32_t i = 1; i < localRows - 1; ++i) {
    std::cout << "Row " << i << " sides: "
              << "Left=" << at(i, 0).cm.m
              << ", Right=" << at(i, localCols - 1).cm.m << "\n";
  }

  // Bottom row
  std::cout << "Bottom row: ";
  for (int32_t j = 1; j <= localCols; ++j) {
    std::cout << at(localRows - 1, j).cm.m << " ";
  }
  std::cout << "\n";
}

template <typename T>
void FlatGrid2D<T>::debugGhostHalo(uint32_t worldRank) const {
  if (worldRank != 2)
    return;

  std::cout << "Ghost halo for rank: " << worldRank << std::endl;

  // Top row
  std::cout << "Top row: ";
  for (int32_t j = 0; j < totalCols; ++j) {
    std::cout << at(0, j).cm.x << " ";
  }
  std::cout << "\n";

  // Left and right columns (excluding corners already printed)
  for (int32_t i = 1; i < totalRows - 1; ++i) {
    std::cout << "Row " << i << " sides: "
              << "Left=" << at(i, 0).cm.x
              << ", Right=" << at(i, totalCols - 1).cm.x << "\n";
  }

  // Bottom row
  std::cout << "Bottom row: ";
  for (int32_t j = 0; j < totalCols; ++j) {
    std::cout << at(totalRows - 1, j).cm.x << " ";
  }
  std::cout << "\n";
}


} // namespace Simulation

#endif // FLATGRID2D_TPP
