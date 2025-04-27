#include "ds.hpp"
using namespace std;

#include "parsim_utils.hpp"

namespace Simulation {

void debugParticle(const Particle p) {
  std::cout << std::fixed
            << std::setprecision(
                   DECIMAL_CASES); // Set fixed-point and precision once
  std::cout << "Particle " << p.idx << ": "
            << "mass=" << p.m << " x=" << p.x << " y=" << p.y << " vx=" << p.vx
            << " vy=" << p.vy << std::endl;
}

void debugParticles(Particle *particles, uint64_t nPart) {
  for (uint64_t i = 0; i < nPart; i++) {
    Particle par = particles[i];
    debugParticle(par);
  }
}

void debugCells(FlatGrid2D<Cell> grid) {
  uint32_t localRows = grid.getLocalRows();
  uint32_t localCols = grid.getLocalCols();
  for (uint32_t row = 1; row <= localRows; row++) {
    for (uint32_t col = 1; col <= localCols; col++) {
      const CenterOfMass cm = grid.at(row, col).cm;
      // Print using fixed formatting with three decimal places
      std::cout << "Cell " << (row - 1) * localCols + (col - 1)
                << " x: " << std::fixed << std::setprecision(DECIMAL_CASES)
                << cm.x << " y: " << std::fixed
                << std::setprecision(DECIMAL_CASES) << cm.y
                << " m: " << std::fixed << std::setprecision(DECIMAL_CASES)
                << cm.m << std::endl;
    }
  }
}

void debugCellForce(const uint64_t &label1, const uint64_t &label2, double dx,
                    double dy, double mag) {
  std::cout << "P" << label1 << "/" << "C" << label2 << " mag: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << mag << " fx: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dx << " fy: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dy << std::endl;
}

void debugParForce(const uint64_t &label1, const uint64_t &label2, double dx,
                   double dy, double mag) {
  std::cout << "P" << label1 << "/" << "P" << label2 << " mag: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << mag << " fx: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dx << " fy: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dy << std::endl;
}

void debugCollision(const uint64_t &label1, const uint64_t &label2, double &d) {
  std::cout << "[Collision] " << "P" << label1 << "/" << "P" << label2
            << " Distance: " << std::setprecision(DECIMAL_CASES) << d << endl;
}

} // namespace Simulation
