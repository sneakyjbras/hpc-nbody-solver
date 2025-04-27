using namespace std;

#include "parsimUtils.hpp"

namespace Simulation {

void debugParticle(const Particle p, uint64_t i) {
  std::cout << std::fixed
            << std::setprecision(
                   DECIMAL_CASES); // Set fixed-point and precision once
  std::cout << "Particle " << i << ": "
            << "mass=" << p.m << " x=" << p.x << " y=" << p.y << " vx=" << p.vx
            << " vy=" << p.vy << std::endl;
}

void debugParticles(Particle *particles, uint64_t nPart) {
  for (uint64_t i = 0; i < nPart; i++) {
    debugParticle(particles[i], i);
  }
}

void debugCells(Cell **cells, uint32_t ncside) {
  for (uint32_t row = 0; row < ncside; row++) {
    for (uint32_t col = 0; col < ncside; col++) {
      const CenterOfMass cm = cells[row][col].cm;
      // Print using fixed formatting with three decimal places
      std::cout << "Cell " << row * ncside + col << " x: " << std::fixed
                << std::setprecision(DECIMAL_CASES) << cm.x
                << " y: " << std::fixed << std::setprecision(DECIMAL_CASES)
                << cm.y << " m: " << std::fixed
                << std::setprecision(DECIMAL_CASES) << cm.m << std::endl;
    }
  }
}

void debugForce(const std::string &label1, const std::string &label2, double dx,
                double dy, double mag) {
  std::cout << label1 << "/" << label2 << " mag: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << mag << " fx: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dx << " fy: " << std::fixed
            << std::setprecision(DECIMAL_CASES) << dy << std::endl;
}

void debugCollision(const std::string &label1, const std::string &label2,
                    double &d) {
  std::cout << "[Collision] " << label1 << "/" << label2
            << " Distance: " << std::setprecision(DECIMAL_CASES) << d << endl;
}
} // namespace Simulation
