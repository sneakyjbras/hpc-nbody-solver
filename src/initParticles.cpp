#include "initParticles.hpp"

namespace Simulation {
unsigned int seed;
void initR4Uni(int inputSeed) { seed = inputSeed + 987654321; }
double rndUniform01() {
  int seed_in = seed;
  seed ^= (seed << 13);
  seed ^= (seed >> 17);
  seed ^= (seed << 5);
  return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}
double rndNormal01() {
  double u1, u2, z, result;
  do {
    u1 = rndUniform01();
    u2 = rndUniform01();
    z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
  } while (result < 0 || result >= 1);
  return result;
}

void initParticles(long userseed, double side, long ncside, long long nPart,
                   Particle *par) {
  double (*rnd01)() = rndUniform01;
  long long i;

  if (userseed < 0) {
    rnd01 = rndNormal01;
    userseed = -userseed;
  }

  initR4Uni(userseed);

  for (i = 0; i < nPart; i++) {
    par[i].x = rnd01() * side;
    par[i].y = rnd01() * side;
    par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
    par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;

    par[i].m = rnd01() * 0.01 * (ncside * ncside) / nPart / G * EPSILON2;
  }
}
} // namespace Simulation
