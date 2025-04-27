#include <omp.h>
#include <unordered_map>

using namespace std;

#include "initParticles.hpp"
#include "parsim.hpp"

namespace Simulation {

Parsim::Parsim(int32_t seed, double side, int32_t ncside, uint64_t nPart,
               int32_t tsteps)
    : seed(seed), side(side), ncside(ncside), nPart(nPart), tsteps(tsteps) {
#pragma omp single
  {
    particles = std::vector<Particle>(nPart);
  }
  initParticles(seed, side, ncside, nPart, particles.data());
}
Parsim::~Parsim() {}

uint64_t Parsim::simulate(std::unordered_set<Particle *> &visited,
                          std::queue<Particle *> &q, const int32_t &idxStart,
                          const int32_t &idxEnd) {
  uint64_t collisions = 0;
  for (int32_t t = 0; t < tsteps; ++t) {
    computeCenterOfMass(idxStart, idxEnd);
#pragma omp barrier
    computeGravitationalPull(idxStart, idxEnd);
#pragma omp barrier
    computePositionVelocity(idxStart, idxEnd);
#pragma omp barrier
    recomputeGrid(idxStart, idxEnd);
#pragma omp barrier
    collisions += computeCollisions(visited, q, idxStart, idxEnd);
#pragma omp barrier
  }
  return collisions;
}

void Parsim::allocateGrid() {

  uint32_t totalCells = 0;
#pragma omp single
  {
    totalCells = static_cast<uint64_t>(ncside) *
                 ncside; // doesn't need to be initialized by N threads
    grid = std::vector<Cell>(
        totalCells); // doesn't need to be initialized by N threads
  }

#pragma omp for schedule(static)
  for (uint64_t i = 0; i < totalCells; ++i) {
    grid[i].particles.reserve((seed >= 0 ? 4 : 8) * nPart / totalCells);
  }
}

void Parsim::populateGrid() {
#pragma omp for schedule(static) // contention of having N threads adding
                                 // particles to the grid or having just one
  for (uint64_t i = 0; i < nPart; i++) {
    // Reset forces for the particle.
    particles[i].fx = 0;
    particles[i].fy = 0;

    // Compute the column and row indices based on the particle's position.
    int32_t col = getIndex(particles[i].x, side, ncside);
    int32_t row = getIndex(particles[i].y, side, ncside);

    // Convert 2D (row, col) into a linear index.
    int32_t idx = row * ncside + col;

    // Insert the particle into the appropriate cell.
    omp_set_lock(&grid[idx].lock);
    grid[idx].addParticle(&particles[i]);
    omp_unset_lock(&grid[idx].lock);
  }
}

void Parsim::recomputeGrid(const int32_t &idxStart, const int32_t &idxEnd) {
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    cell.processIncoming(); // merge incoming particles into this cell
    // Reset forces and clear collision lists for all particles in the cell
    for (Particle *p : cell.particles) {
      p->resetForce();
      p->clearColliding();
    }
  }
}

void Parsim::computeCenterOfMass(const int32_t &idxStart,
                                 const int32_t &idxEnd) {
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    const auto &particles = cell.particles;
    cell.cm.reset();
    // Sum masses and weighted positions of all particles in this cell
    for (Particle *p : particles) {
      double mass = p->m;
      cell.cm.computeMass(mass);
      cell.cm.computePos(mass * p->x, mass * p->y);
    }
    cell.cm.normalizePos();
  }
}

void Parsim::computeIncellPull(auto &particles) {
  const double EPS = 1e-18;
  for (size_t i = 0; i < particles.size(); ++i) {
    Particle *pA = particles[i];
    double ax = pA->x, ay = pA->y;
    for (size_t j = i + 1; j < particles.size(); ++j) {
      Particle *pB = particles[j];
      // Compute wrapped differences preserving sign.
      double dx = pB->x - ax;
      if (fabs(dx) > side - fabs(dx))
        dx = (dx < 0) ? dx + side : dx - side;

      double dy = pB->y - ay;
      if (fabs(dy) > side - fabs(dy))
        dy = (dy < 0) ? dy + side : dy - side;

      double d2 = dx * dx + dy * dy;
      if (d2 < EPS)
        continue; // Avoid division by zero

      // Compute inverse distance using one sqrt call.
      double invD = 1.0 / sqrt(d2); // or use a fast inverse sqrt approximation
      double invD3 = invD * invD * invD;
      double fx = dx * G * (pA->m * pB->m) * invD3;
      double fy = dy * G * (pA->m * pB->m) * invD3;

      pA->addForce(fx, fy);
      pB->addForce(-fx, -fy);
    }
  }
}

void Parsim::computeExtracellPull(auto &particles, CenterOfMass cm, bool wrapX,
                                  bool wrapY) {
  const double EPS = 1e-18;
  for (size_t i = 0; i < particles.size(); ++i) {
    Particle *par = particles[i];
    double deltaX = cm.x - par->x;
    double deltaY = cm.y - par->y;

    if (wrapX) {
      if (fabs(deltaX) > side - fabs(deltaX))
        deltaX = (deltaX < 0) ? deltaX + side : deltaX - side;
    }
    if (wrapY) {
      if (fabs(deltaY) > side - fabs(deltaY))
        deltaY = (deltaY < 0) ? deltaY + side : deltaY - side;
    }

    double d2 = deltaX * deltaX + deltaY * deltaY;
    if (d2 < EPS)
      continue; // Avoid division by zero

    double invD = 1.0 / sqrt(d2); // Only one sqrt here.
    double invD3 = invD * invD * invD;
    double fx = deltaX * G * (par->m * cm.m) * invD3;
    double fy = deltaY * G * (par->m * cm.m) * invD3;

    par->addForce(fx, fy);
  }
}

void Parsim::computeGravitationalPull(const int32_t &idxStart,
                                      const int32_t &idxEnd) {
  // 8 neighboring cell offsets (row_offset, col_offset) for adjacent cells
  const std::pair<int32_t, int32_t> offsets[8] = {
      {-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    if (cell.particles.empty())
      continue;
    // Intra-cell interactions (same cell)
    computeIncellPull(cell.particles);
    // Compute neighbor interactions using center-of-mass of adjacent cells
    int32_t row = idx / ncside;
    int32_t col = idx % ncside;
    for (auto offset : offsets) {
      int32_t adjRow = row + offset.first;
      int32_t adjCol = col + offset.second;
      bool wrapX = false, wrapY = false;
      // Wrap around boundaries if needed (toroidal grid)
      if (adjRow < 0 || adjRow >= ncside) {
        wrapY = true;
        computeWrap(adjRow, ncside); // wrap row index within [0, ncside)
      }
      if (adjCol < 0 || adjCol >= ncside) {
        wrapX = true;
        computeWrap(adjCol, ncside); // wrap col index within [0, ncside)
      }
      // Compute linear index of the neighbor and apply extra-cell pull
      int32_t neighborIdx = adjRow * ncside + adjCol;
      Cell &neighbor = grid[neighborIdx];
      computeExtracellPull(cell.particles, neighbor.cm, wrapX, wrapY);
    }
  }
}

// Example function that updates particles and then merges incoming particles
// into the global grid, using a per-cell lock for safe updates.
void Parsim::computePositionVelocity(const int32_t &idxStart,
                                     const int32_t &idxEnd) {
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    if (cell.particles.empty())
      continue;
    // Iterate backwards so we can remove elements safely by index
    for (int64_t i = cell.particles.size() - 1; i >= 0; --i) {
      Particle *par = cell.particles[i];
      par->computePos(side);
      par->computeVelocity();
      // Determine the particle's new cell coordinates
      int32_t new_row = getIndex(par->y, side, ncside);
      int32_t new_col = getIndex(par->x, side, ncside);
      int32_t newIdx = new_row * ncside + new_col;
      if (newIdx != idx) {
        // Remove particle from current cell
        cell.removeParticle(i);
        // Lock the destination cell and add particle to its incoming queue
        omp_set_lock(&grid[newIdx].lock);
        grid[newIdx].addIncoming(par);
        omp_unset_lock(&grid[newIdx].lock);
      }
    }
  }
}

uint64_t Parsim::computeCollisions(std::unordered_set<Particle *> &visited,
                                   std::queue<Particle *> &q,
                                   const int32_t idxStart,
                                   const int32_t idxEnd) {
  uint64_t collisionCount = 0;
  // 1. Detect collisions and mark colliding particles
  const double COLLISION_RADIUS = 5e-3;
  const double thresh2 = COLLISION_RADIUS * COLLISION_RADIUS;
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    uint64_t numParticles = cell.particles.size();
    if (numParticles < 2)
      continue;
    // Sort particles by x-coordinate for sweep-and-prune
    std::sort(cell.particles.begin(), cell.particles.end(),
              [](const Particle *a, const Particle *b) { return a->x < b->x; });
    // Check each pair for collision
    for (uint64_t i = 0; i < numParticles; ++i) {
      Particle *pA = cell.particles[i];
      double ax = pA->x, ay = pA->y;
      for (uint64_t j = i + 1; j < numParticles; ++j) {
        Particle *pB = cell.particles[j];
        double dx = pB->x - ax;
        if (dx > COLLISION_RADIUS)
          break; // no more possible collisions in this cell
        double dy = fabs(ay - pB->y);
        if (dy > COLLISION_RADIUS)
          continue;
        double dist2 = dx * dx + dy * dy;
        if (dist2 <= thresh2) {
          pA->addColliding(pB);
          pB->addColliding(pA);
        }
      }
    }
  }
  // 2. Resolve collisions (group particles and remove from cells)
  for (int32_t idx = idxStart; idx < idxEnd; ++idx) {
    Cell &cell = grid[idx];
    int64_t numParticles = cell.particles.size();
    if (numParticles < 2)
      continue;
    // Iterate backwards when removing particles from the vector
    for (int64_t i = numParticles - 1; i >= 0; --i) {
      Particle *par = cell.particles[i];
      if (par->collidingParticles.empty())
        continue;
      cell.removeParticle(i);
      par->resetMass();
      if (visited.count(par))
        continue;
      // BFS to group all connected colliding particles
      q.push(par);
      visited.insert(par);
      while (!q.empty()) {
        Particle *current = q.front();
        q.pop();
        // (Optional: gather current in a group vector if needed)
        for (Particle *nbr : current->collidingParticles) {
          if (!visited.count(nbr)) {
            visited.insert(nbr);
            q.push(nbr);
          }
        }
      }
      collisionCount++;
    }
  }
  // Clean up for next timestep
  visited.clear();
  while (!q.empty())
    q.pop();
  return collisionCount;
}
} // namespace Simulation
