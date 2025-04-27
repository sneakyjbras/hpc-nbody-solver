#include "ds.hpp"
#include "parsimUtils.hpp"
#include <cstdint>
using namespace std;

#include "parsim.hpp"

namespace Simulation {

// Define thread_local variables.
thread_local int32_t Parsim::rowStart;
thread_local int32_t Parsim::rowEnd;
thread_local int32_t Parsim::colStart;
thread_local int32_t Parsim::colEnd;
thread_local uint16_t Parsim::tid;

Parsim::Parsim(int32_t inputSeed, double side, int32_t ncside, uint64_t nPart,
               int32_t tsteps, uint16_t nThreads, MPI_Comm worldComm)
    : inputSeed(inputSeed), side(side), ncside(ncside), nPart(nPart),
      tsteps(tsteps), nThreads(nThreads), worldComm(worldComm) {
  int32_t err;

  // Create global communicator rank and size.
  err = MPI_Comm_rank(worldComm, &worldRank);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Comm_rank failed in global communicator."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Comm_size(worldComm, &worldSize);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Comm_size failed in global communicator."
              << std::endl;
    MPI_Abort(worldComm, err);
  }

  // Set flag for special 2-process checkered decomposition.
  isCheckeredTwoProc = (worldSize == 2);

  // Create cartesian dimensions.
  err = MPI_Dims_create(worldSize, 2, dims);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Dims_create failed." << std::endl;
    MPI_Abort(worldComm, err);
  }

  // Create cartesian communicator.
  err = MPI_Cart_create(worldComm, 2, dims, periods, 0, &cartComm);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Cart_create failed." << std::endl;
    MPI_Abort(worldComm, err);
  }

  // Create cartesian rank and size.
  err = MPI_Comm_rank(cartComm, &cartRank);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Comm_rank failed in cartesian communicator."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Comm_size(cartComm, &cartSize);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Comm_size failed in cartesian communicator."
              << std::endl;
    MPI_Abort(worldComm, err);
  }

  // Get cartesian coordinates.
  err = MPI_Cart_coords(cartComm, cartRank, 2, coords);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Cart_coords failed." << std::endl;
    MPI_Abort(worldComm, err);
  }
}

Parsim::~Parsim() {}

// Public

uint64_t Parsim::simulate(std::unordered_set<Particle *> &visited,
                          std::queue<Particle *> &q) {
  uint64_t collisions = 0;
  for (int32_t t = 0; t < tsteps; ++t) {
    computeCenterOfMass();
#pragma omp barrier
#pragma omp master
    {
      exchangeCOM();
    }
#pragma omp barrier
    computeGravitationalPull();
#pragma omp barrier
    computePositionVelocity();
#pragma omp barrier
#pragma omp master
    {
      exchangePar();
    }
#pragma omp barrier
    recomputeGrid();
#pragma omp barrier
    collisions += computeCollisions(visited, q);
#pragma omp barrier
  }

  // debugParticles(particles.data(), particles.size());
  return collisions;
}

// Assuming FlatGrid2D, Cell, and other types are already defined.

void Parsim::allocateGrid() {
  uint64_t reserveCount = 0;
  uint32_t neighCount = ranks.size();
#pragma omp single
  {

    for (uint32_t i = 0; i < neighCount; ++i)
      omp_init_lock(&rankLocks[i]);

    // Initialize the grid with the local subdomain dimensions.
    grid = FlatGrid2D<Cell>(sub.localRows, sub.localCols);

    // Estimate the average number of particles per cell.
    // nPart/worldSize is the approximate number of particles per process.
    std::uint64_t avgParticlesPerCell =
        (nPart / worldSize) /
        (static_cast<std::uint64_t>(sub.localRows) * sub.localCols);
    // Use a factor of 4 for uniform or 8 for normal distributions (based on
    // inputSeed).
    reserveCount = (inputSeed >= 0 ? 4 : 8) * (avgParticlesPerCell + 1);
  }
  // For each cell in the local grid, reserve space for its particle vector.
  // Also, initialize each cell’s COM index.
#pragma omp for schedule(static) collapse(2)
  for (int32_t i = 1; i <= sub.localRows; ++i) {
    for (int32_t j = 1; j <= sub.localCols; ++j) {
      grid.at(i, j).particles.reserve(reserveCount);
      grid.at(i, j).cm.idx =
          static_cast<uint64_t>((i - 1) * sub.localCols + (j - 1));
    }
  }

#pragma omp master
  {
    // Reserve capacity for COM (Center-Of-Mass) exchange buffers.
    // Reserve enough space to cover two rows (top and bottom), two columns
    // (left and right), plus 4 for the corners.
    sendComBuf.reserve(static_cast<std::uint64_t>(sub.localRows) * 2 +
                       static_cast<std::uint64_t>(sub.localCols) * 2 + 4);
    recvComBuf.reserve(static_cast<std::uint64_t>(sub.localRows) * 2 +
                       static_cast<std::uint64_t>(sub.localCols) * 2 + 4);

    // Reserve space for counts and displacement arrays for COM communication.
    sendCounts.reserve(neighCount);
    recvCounts.reserve(neighCount);
    sComDispls.reserve(neighCount);
    rComDispls.reserve(neighCount);

    // Reserve neighbor exchange vectors for particle migration.
    parSendCounts.reserve(neighCount);
    parSDispls.reserve(neighCount);
    parRecvCounts.reserve(neighCount);
    parRDispls.reserve(neighCount);

    // Estimate maximum number of local particles (with an extra margin).
    std::uint64_t maxLocalParts =
        (nPart / worldSize) * (inputSeed >= 0 ? 4 : 8);
    // Reserve buffers for sending and receiving particles during migration.
    parSendBuffer.reserve(maxLocalParts);
    parRecvBuffer.reserve(maxLocalParts);
  }
#pragma omp single
  {
    // Prepare the neighbor map (particleNeighborMap) for tracking particles
    // that cross domain boundaries.
    particleNeighborMap.clear();
    for (int32_t nbr : ranks) {
      particleNeighborMap[nbr].reserve(reserveCount / neighCount * 4);
    }
  }
}

void Parsim::populateGrid() {
#pragma omp for schedule(static)
  for (uint64_t i = 0; i < particles.size(); i++) {
    // Get global grid indices for the particle
    int32_t col = getIndex(particles[i].x, side, ncside);
    int32_t row = getIndex(particles[i].y, side, ncside);

    // Convert global to local coordinates
    int32_t localCol = col - sub.colStart;
    int32_t localRow = row - sub.rowStart;

    // Add particle to appropriate cell
    omp_set_lock(&grid.at(localRow + 1, localCol + 1).lock);
    grid.at(localRow + 1, localCol + 1).addParticle(&particles[i]);
    omp_unset_lock(&grid.at(localRow + 1, localCol + 1).lock);
  }
}
void Parsim::getAllNeighbors() {
  int32_t err;
  // --- Compute neighbor rank for each ghost direction ---
  for (const auto &[dir, offset] : directions) {
    int32_t newRow = (coords[0] + offset.first + dims[0]) % dims[0];
    int32_t newCol = (coords[1] + offset.second + dims[1]) % dims[1];

    int32_t neighborCoords[2] = {newRow, newCol};
    int32_t neighborRank;
    err = MPI_Cart_rank(cartComm, neighborCoords, &neighborRank);
    if (err != MPI_SUCCESS) {
      std::cerr << "Error: MPI_Cart_rank failed for neighbor direction "
                << static_cast<int32_t>(dir) << std::endl;
      MPI_Abort(worldComm, err);
    }
    neighbors[dir] = neighborRank;
  }

  // --- Build neighborDirectionMap, ranks, and weights ---
  for (const auto &neighbor : neighbors) {
    // If this neighbor hasn't been added yet, add its rank and set weight=1.
    if (!neighborDirectionMap.contains(neighbor.second)) {
      ranks.push_back(neighbor.second);
      weights.push_back(1);
    }
    // Associate this direction with the neighbor.
    neighborDirectionMap[neighbor.second].push_back(neighbor.first);
  }
  std::sort(ranks.begin(), ranks.end());
  outdegree = static_cast<int32_t>(weights.size());
  indegree = outdegree;

  // --- Create the distributed graph communicator ---
  err = MPI_Dist_graph_create(cartComm,
                              1, // Number of source nodes (only one here)
                              &worldRank,     // Source node array
                              &outdegree,     // Degree of source node
                              ranks.data(),   // Destination ranks
                              weights.data(), // Weights for each edge
                              MPI_INFO_NULL,  // No special info
                              0,              // No reordering
                              &graphComm);    // Output communicator
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Dist_graph_create failed." << std::endl;
    MPI_Abort(worldComm, err);
  }
}

void Parsim::initializeCommunication() {

  int32_t err;

  // Process each neighbor in the neighborDirectionMap.
  for (const auto &rank : ranks) {
    int32_t count = 0;
    // Create a local mapping from ghost direction to a lambda returning the
    // count.
    auto countFunc =
        std::unordered_map<GhostDirection, std::function<int32_t()>>{
            {GhostDirection::N, [this]() { return sub.localCols; }},
            {GhostDirection::S, [this]() { return sub.localCols; }},
            {GhostDirection::E, [this]() { return sub.localRows; }},
            {GhostDirection::W, [this]() { return sub.localRows; }},
            {GhostDirection::NE, []() { return 1; }},
            {GhostDirection::SE, []() { return 1; }},
            {GhostDirection::SW, []() { return 1; }},
            {GhostDirection::NW, []() { return 1; }}};

    for (const auto &direction : neighborDirectionMap[rank]) {
      // Exclude diagonal ghost regions for 2-process checkered mode.
      if (isCheckeredTwoProc &&
          (direction == GhostDirection::NE || direction == GhostDirection::NW ||
           direction == GhostDirection::SE ||
           direction == GhostDirection::SW)) {
        continue;
      }
      count += countFunc.at(direction)();
    }
    sendCounts.push_back(count);
  }

  // Compute displacements based on the sendCounts.
  sComDispls.resize(sendCounts.size());
  if (!sendCounts.empty()) {
    sComDispls[0] = 0;
    for (uint32_t i = 1; i < sendCounts.size(); i++) {
      sComDispls[i] = sComDispls[i - 1] + sendCounts[i - 1];
    }
  }

  int32_t indegree, outdegree, weighted;
  // Get the number of in-neighbors and out-neighbors.
  err = MPI_Dist_graph_neighbors_count(graphComm, &indegree, &outdegree,
                                       &weighted);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Dist_graph_neighbors_count failed." << std::endl;
    MPI_Abort(worldComm, err);
  }

  std::vector<int32_t> sources(indegree);
  std::vector<int32_t> sourceWeights(indegree);
  std::vector<int32_t> destinations(outdegree);
  std::vector<int32_t> destWeights(outdegree);

  // Retrieve the neighbor lists.
  err = MPI_Dist_graph_neighbors(graphComm, indegree, sources.data(),
                                 sourceWeights.data(), outdegree,
                                 destinations.data(), destWeights.data());
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Dist_graph_neighbors failed." << std::endl;
    MPI_Abort(worldComm, err);
  }

  std::vector<int32_t> computedRecvCounts(sendCounts.size(), 0);
  err = MPI_Neighbor_alltoall(sendCounts.data(), 1, MPI_INT,
                              computedRecvCounts.data(), 1, MPI_INT, graphComm);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Neighbor_alltoall failed." << std::endl;
    MPI_Abort(worldComm, err);
  }

  recvCounts = computedRecvCounts;

  // Compute receive displacements from recvCounts.
  rComDispls.resize(recvCounts.size());
  if (!recvCounts.empty()) {
    rComDispls[0] = 0;
    for (uint32_t i = 1; i < recvCounts.size(); ++i) {
      rComDispls[i] = rComDispls[i - 1] + recvCounts[i - 1];
    }
  }

  // Resize the receive buffer to hold all incoming data.
  uint32_t totalRecv = 0;
  for (uint32_t cnt : recvCounts)
    totalRecv += cnt;
  recvComBuf.resize(totalRecv);
}

void Parsim::createMPIType() {
  int32_t err;

  // --- Create MPI datatype for CenterOfMass ---
  CenterOfMass dummyCom;
  const int32_t comFieldCount = 5;
  int32_t comBlockLengths[comFieldCount] = {1, 1, 1, 1, 1};
  MPI_Aint comDisplacements[comFieldCount] = {};
  MPI_Datatype comTypes[comFieldCount] = {MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE,
                                          MPI_DOUBLE, MPI_INT};

  MPI_Aint comBaseAddress;
  err = MPI_Get_address(&dummyCom, &comBaseAddress);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for CenterOfMass base."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyCom.idx, &comDisplacements[0]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for dummyCom.idx." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyCom.m, &comDisplacements[1]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for dummyCom.m." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyCom.x, &comDisplacements[2]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for dummyCom.x." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyCom.y, &comDisplacements[3]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for dummyCom.y." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyCom.d, &comDisplacements[4]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for dummyCom.d." << std::endl;
    MPI_Abort(worldComm, err);
  }
  for (int32_t i = 0; i < comFieldCount; ++i) {
    comDisplacements[i] -= comBaseAddress;
  }
  err = MPI_Type_create_struct(comFieldCount, comBlockLengths, comDisplacements,
                               comTypes, &cmType);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Type_create_struct failed for CenterOfMass."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Type_commit(&cmType);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Type_commit failed for CenterOfMass." << std::endl;
    MPI_Abort(worldComm, err);
  }

  // --- Create MPI datatype for Particle ---
  const int32_t pFieldCount = 8;
  int32_t pBlockLengths[pFieldCount] = {1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype pTypes[pFieldCount] = {MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE,   MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE,   MPI_DOUBLE};
  MPI_Aint pDisplacements[pFieldCount] = {};

  Particle dummyParticle;
  MPI_Aint pBase;
  err = MPI_Get_address(&dummyParticle, &pBase);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle base."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.idx, &pDisplacements[0]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.idx." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.x, &pDisplacements[1]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.x." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.y, &pDisplacements[2]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.y." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.vx, &pDisplacements[3]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.vx." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.vy, &pDisplacements[4]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.vy." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.m, &pDisplacements[5]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.m." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.fx, &pDisplacements[6]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.fx." << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Get_address(&dummyParticle.fy, &pDisplacements[7]);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Get_address failed for Particle.fy." << std::endl;
    MPI_Abort(worldComm, err);
  }
  for (int32_t i = 0; i < pFieldCount; ++i) {
    pDisplacements[i] -= pBase;
  }

  MPI_Datatype rawParticleType;
  err = MPI_Type_create_struct(pFieldCount, pBlockLengths, pDisplacements,
                               pTypes, &rawParticleType);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Type_create_struct failed for Particle."
              << std::endl;
    MPI_Abort(worldComm, err);
  }

  MPI_Aint pExtent = sizeof(Particle);
  MPI_Datatype resizedType;
  err = MPI_Type_create_resized(rawParticleType, 0, pExtent, &resizedType);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Type_create_resized failed for Particle."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  err = MPI_Type_commit(&resizedType);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Type_commit failed for Particle." << std::endl;
    MPI_Abort(worldComm, err);
  }

  MPI_Type_free(&rawParticleType);
  particleType = resizedType;
}

// Private

// 1: Simulation Steps

void Parsim::computeCenterOfMass() {
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
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
}

void Parsim::computeIncellPull(auto &particles) {
  static constexpr double EPS = 1e-18;
  for (uint64_t i = 0; i < particles.size(); ++i) {
    Particle *pA = particles[i];
    double ax = pA->x, ay = pA->y;
    for (uint64_t j = i + 1; j < particles.size(); ++j) {
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

      // if (worldRank == 0)
      // if (!pA->idx)
      // debugParForce(pA->idx, pB->idx, fx, fy, invD);
    }
  }
}

void Parsim::computeExtracellPull(std::vector<Particle *> &particles,
                                  const CenterOfMass &cm, int32_t adjRow,
                                  int32_t adjCol, double side) {
  static constexpr double EPS = 1e-18;

  for (uint64_t i = 0; i < particles.size(); ++i) {
    Particle *par = particles[i];
    double deltaX = cm.x - par->x;
    double deltaY = cm.y - par->y;

    bool wrapX = false;
    if (adjCol < 1) {
      if (cm.x > par->x) {
        wrapX = true;
      }
    } else if (adjCol > sub.localCols) {
      if (cm.x < par->x) {
        wrapX = true;
      }
    }

    bool wrapY = false;
    if (adjRow < 1) {
      if (cm.y > par->y) {
        wrapY = true;
      }
    } else if (adjRow > sub.localRows) {
      if (cm.y < par->y) {
        wrapY = true;
      }
    }
    if (wrapX) {
      deltaX = (deltaX < 0) ? deltaX + side : deltaX - side;
    }
    if (wrapY) {
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

    // if (worldRank == 0)
    //   if (!par->idx)
    //     debugCellForce(par->idx, cm.idx, fx, fy, invD);
  }
}

void Parsim::computeGravitationalPull() {
  // 8 neighbor offsets: top-left, top, top-right, left, right, bottom-left,
  // bottom, bottom-right
  const std::pair<int32_t, int32_t> offsets[8] = {
      {-1, -1}, {0, -1}, {1, -1}, {-1, 0}, {1, 0}, {-1, 1}, {0, 1}, {1, 1}};

  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
      if (cell.particles.empty()) {
        continue;
      }
      // Intra-cell interactions
      computeIncellPull(cell.particles);
    }
  }
#pragma omp barrier
#pragma omp master
  {
    extractAndUpdateGhostRegions();
  }
#pragma omp barrier
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
      if (cell.particles.empty()) {
        continue;
      }

      // For each neighbor offset
      for (auto offset : offsets) {
        int32_t adjRow = i + offset.first;
        int32_t adjCol = j + offset.second;

        // Retrieve the neighbor cell (could be out-of-bounds => ghost)
        Cell &neighbor = grid.at(adjRow, adjCol);

        // Pass adjRow and adjCol to extracell
        computeExtracellPull(cell.particles, neighbor.cm, adjRow, adjCol, side);
      }
    }
  }
}

void Parsim::recomputeGrid() {
  // First, process local cells: merge any internal incoming particles
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
      cell.processIncoming(); // Merge internal incoming particles.
      // Reset forces and clear collision lists.
      for (Particle *p : cell.particles) {
        p->resetForce();
        p->clearColliding();
      }
    }
  }
#pragma omp master
  {
    MPI_Wait(&requestPar, MPI_STATUS_IGNORE);
  }
#pragma omp barrier

  // Now process the receive buffer.
#pragma omp for schedule(static)
  for (uint32_t i = 0; i < parRecvBuffer.size(); ++i) {
    Particle *newParticle = new Particle(parRecvBuffer[i]);

    // Compute global grid indices for the particle.
    int32_t col = getIndex(newParticle->x, side, ncside);
    int32_t row = getIndex(newParticle->y, side, ncside);

    // Convert global indices to local indices.
    int32_t localCol = col - sub.colStart;
    int32_t localRow = row - sub.rowStart;

    // Add the new particle to the appropriate cell.
    omp_set_lock(&grid.at(localRow + 1, localCol + 1).lock);
    grid.at(localRow + 1, localCol + 1).addParticle(newParticle);
    omp_unset_lock(&grid.at(localRow + 1, localCol + 1).lock);
  }

#pragma omp master
  {
    // Clear com vectors
    sendComBuf.clear();
    recvComBuf.clear();
    sendCounts.clear();
    recvCounts.clear();
    sComDispls.clear();
    rComDispls.clear();

    // Clear par vectors
    parSendBuffer.clear();
    parSendCounts.clear();
    parSDispls.clear();

    parRecvBuffer.clear();
    parRecvCounts.clear();
    parRDispls.clear();

    // Clear neighbor map
    for (auto &entry : particleNeighborMap) {
      entry.second.clear();
    }
    particleNeighborMap.clear();
  }
}

// Example function that updates particles and then merges incoming particles
// into the global grid, using a per-cell lock for safe updates.
void Parsim::computePositionVelocity() {
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
      if (cell.particles.empty())
        continue;
      // Iterate backwards so we can remove elements safely by index
      for (int64_t i = cell.particles.size() - 1; i >= 0; --i) {
        Particle *par = cell.particles[i];
        par->computePos(side);
        par->computeVelocity();

        // Determine the particle's new cell coordinates
        int32_t newRow = getIndex(par->y, side, ncside);
        int32_t newCol = getIndex(par->x, side, ncside);

        if ((newRow >= sub.rowStart && newRow < sub.rowEnd) &&
            (newCol >= sub.colStart && newCol < sub.colEnd)) {
          // Convert global to local coordinates
          int32_t localCol = newCol - sub.colStart;
          int32_t localRow = newRow - sub.rowStart;
          if ((localCol + 1 != j) || (localRow + 1 != i)) {
            cell.removeParticle(i);
            omp_set_lock(&grid.at(localRow + 1, localCol + 1).lock);
            grid.at(localRow + 1, localCol + 1).addIncoming(par);
            omp_unset_lock(&grid.at(localRow + 1, localCol + 1).lock);
          }
        } else {

          // Compute the final destination rank directly.
          int32_t destRank = getOwnerRank(newRow, newCol);
          // Remove the particle from the current cell.
          cell.removeParticle(i);
          // Directly add it to a new container keyed by destRank.
          omp_set_lock(&rankLocks[destRank]);
          particleNeighborMap[destRank].push_back(*par);
          omp_unset_lock(&rankLocks[destRank]);
        }
      }
    }
  }
}

uint64_t Parsim::computeCollisions(std::unordered_set<Particle *> &visited,
                                   std::queue<Particle *> &q) {
  uint64_t collisionCount = 0;
  // 1. Detect collisions and mark colliding particles
  static constexpr double COLLISION_RADIUS = 5e-3;
  static constexpr double thresh2 = COLLISION_RADIUS * COLLISION_RADIUS;
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
      uint64_t numParticles = cell.particles.size();
      if (numParticles < 2)
        continue;
      // Sort particles by x-coordinate for sweep-and-prune
      std::sort(
          cell.particles.begin(), cell.particles.end(),
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

          // debugCollision(pA->idx, pB->idx, dist2);
          if (dist2 <= thresh2) {
            pA->addColliding(pB);
            pB->addColliding(pA);
          }
        }
      }
    }
  }
  // 2. Resolve collisions (group particles and remove from cells)
  for (int32_t i = rowStart; i <= rowEnd; i++) {
    for (int32_t j = colStart; j <= colEnd; j++) {
      Cell &cell = grid.at(i, j);
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
  }
  // Clean up for next timestep
  visited.clear();
  while (!q.empty())
    q.pop();
  return collisionCount;
}

// 2: Particle and Grid

void Parsim::computeSubdomain() {
  // MPI Subdomain
  int32_t rowsPerProc = ncside / dims[0];
  int32_t extraRows = ncside % dims[0];

  int32_t colsPerProc = ncside / dims[1];
  int32_t extraCols = ncside % dims[1];

  sub.rowStart = coords[0] * rowsPerProc + std::min(coords[0], extraRows);
  sub.rowEnd = sub.rowStart + rowsPerProc + (coords[0] < extraRows ? 1 : 0);

  sub.colStart = coords[1] * colsPerProc + std::min(coords[1], extraCols);
  sub.colEnd = sub.colStart + colsPerProc + (coords[1] < extraCols ? 1 : 0);

  sub.localRows = sub.rowEnd - sub.rowStart;
  sub.localCols = sub.colEnd - sub.colStart;
}

void Parsim::computeThreadIndexes(const uint16_t &threadId) {
  // Total threads in this MPI process.
  uint16_t T = nThreads;

  // Factor T into two numbers tRow and tCol such that tRow * tCol == T.
  // For power-of-two T, one common strategy is to try to get them as equal as
  // possible.
  uint16_t tRow = (uint16_t)std::sqrt(T);
  while (T % tRow != 0) { // Adjust tRow until it divides T evenly.
    tRow--;
  }
  uint16_t tCol = T / tRow;

  // Determine this thread's 2D index.
  uint16_t threadRow = threadId / tCol;
  uint16_t threadCol = threadId % tCol;

  // Partition the local subdomain grid (sub.localRows x sub.localCols) among
  // tRow and tCol.
  uint32_t baseRows = sub.localRows / tRow;
  uint32_t extraRows = sub.localRows % tRow;
  uint32_t baseCols = sub.localCols / tCol;
  uint32_t extraCols = sub.localCols % tCol;

  // Calculate the starting row for this thread.
  uint32_t startRow =
      threadRow * baseRows + std::min((uint32_t)threadRow, extraRows);
  uint32_t rowsForThread = baseRows + (threadRow < extraRows ? 1 : 0);
  uint32_t endRow = startRow + rowsForThread - 1;

  // Calculate the starting column for this thread.
  uint32_t startCol =
      threadCol * baseCols + std::min((uint32_t)threadCol, extraCols);
  uint32_t colsForThread = baseCols + (threadCol < extraCols ? 1 : 0);
  uint32_t endCol = startCol + colsForThread - 1;

  // Convert 0-indexed coordinates to 1-indexed, assuming the rest of your code
  // uses 1-indexing.
  rowStart = startRow + 1;
  rowEnd = endRow + 1;
  colStart = startCol + 1;
  colEnd = endCol + 1;
}

double Parsim::rndUniform01() {
  int32_t seedIn = seed;
  seed ^= (seed << 13);
  seed ^= (seed >> 17);
  seed ^= (seed << 5);
  return 0.5 + 0.2328306e-09 * (seedIn + (int32_t)seed);
}

double Parsim::rndNormal01() {
  double u1, u2, z, result;
  do {
    u1 = rndUniform01();
    u2 = rndUniform01();
    z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
  } while (result < 0 || result >= 1);
  return result;
}

void Parsim::initParticles() {
  std::function<double()> rnd01 = [this]() { return rndUniform01(); };
  uint64_t i;

  if (inputSeed < 0) {
    rnd01 = [this]() { return rndNormal01(); };
    inputSeed = -inputSeed;
  }

  initR4Uni(inputSeed);

  for (i = 0; i < nPart; i++) {
    double x = rnd01() * side;
    double y = rnd01() * side;

    // Compute the column and row indices based on the particle's position.
    int32_t col = getIndex(x, side, ncside);
    int32_t row = getIndex(y, side, ncside);
    double vx = (rnd01() - 0.5) * side / ncside / 5.0;
    double vy = (rnd01() - 0.5) * side / ncside / 5.0;
    double m = rnd01() * 0.01 * (ncside * ncside) / nPart / G * EPSILON2;
    if ((row >= sub.rowStart && row < sub.rowEnd) &&
        (col >= sub.colStart && col < sub.colEnd)) {
      Particle particle(i, x, y, vx, vy, m, 0, 0);
      particles.push_back(particle);
    }
  }
}

// 3: MPI Communication

void Parsim::exchangeCOM() {
  // Loop over each neighbor.
  for (const auto &rank : ranks) {
    // For each ghost direction associated with this neighbor.
    for (const auto &direction : neighborDirectionMap[rank]) {
      // Skip diagonal directions in 2-process checkered mode.
      if (isCheckeredTwoProc &&
          (direction == GhostDirection::NE || direction == GhostDirection::NW ||
           direction == GhostDirection::SE ||
           direction == GhostDirection::SW)) {
        continue;
      }
      // Table-driven mapping: each ghost direction is associated with a lambda
      // that pushes the appropriate COM messages into sendComBuf.
      static const std::unordered_map<
          GhostDirection, std::function<void(Parsim *, GhostDirection)>>
          sendFunc = {
              {GhostDirection::N,
               [](Parsim *self, GhostDirection d) {
                 // Send the entire top row.
                 for (int32_t j = 1; j <= self->sub.localCols; ++j) {
                   CenterOfMass msg = self->grid.at(1, j).cm;
                   msg.d = d;
                   self->sendComBuf.push_back(msg);
                 }
               }},
              {GhostDirection::NE,
               [](Parsim *self, GhostDirection d) {
                 // Send top-right corner.
                 CenterOfMass msg = self->grid.at(1, self->sub.localCols).cm;
                 msg.d = d;
                 self->sendComBuf.push_back(msg);
               }},
              {GhostDirection::E,
               [](Parsim *self, GhostDirection d) {
                 // Send the entire right column.
                 for (int32_t i = 1; i <= self->sub.localRows; ++i) {
                   CenterOfMass msg = self->grid.at(i, self->sub.localCols).cm;
                   msg.d = d;
                   self->sendComBuf.push_back(msg);
                 }
               }},
              {GhostDirection::SE,
               [](Parsim *self, GhostDirection d) {
                 // Send bottom-right corner.
                 CenterOfMass msg =
                     self->grid.at(self->sub.localRows, self->sub.localCols).cm;
                 msg.d = d;
                 self->sendComBuf.push_back(msg);
               }},
              {GhostDirection::S,
               [](Parsim *self, GhostDirection d) {
                 // Send the entire bottom row.
                 for (int32_t j = 1; j <= self->sub.localCols; ++j) {
                   CenterOfMass msg = self->grid.at(self->sub.localRows, j).cm;
                   msg.d = d;
                   self->sendComBuf.push_back(msg);
                 }
               }},
              {GhostDirection::SW,
               [](Parsim *self, GhostDirection d) {
                 // Send bottom-left corner.
                 CenterOfMass msg = self->grid.at(self->sub.localRows, 1).cm;
                 msg.d = d;
                 self->sendComBuf.push_back(msg);
               }},
              {GhostDirection::W,
               [](Parsim *self, GhostDirection d) {
                 // Send the entire left column.
                 for (int32_t i = 1; i <= self->sub.localRows; ++i) {
                   CenterOfMass msg = self->grid.at(i, 1).cm;
                   msg.d = d;
                   self->sendComBuf.push_back(msg);
                 }
               }},
              {GhostDirection::NW, [](Parsim *self, GhostDirection d) {
                 // Send top-left corner.
                 CenterOfMass msg = self->grid.at(1, 1).cm;
                 msg.d = d;
                 self->sendComBuf.push_back(msg);
               }}};

      sendFunc.at(direction)(this, direction);
    }
  }

  int32_t err = MPI_Ineighbor_alltoallv(
      sendComBuf.data(), sendCounts.data(), sComDispls.data(), cmType,
      recvComBuf.data(), recvCounts.data(), rComDispls.data(), cmType,
      graphComm, &requestCOM);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Ineighbor_alltoallv failed." << std::endl;
    MPI_Abort(worldComm, err);
  }
}

void Parsim::exchangePar() {
  //  Fill parSendCounts and parSendBuffer based on your neighbor map...
  for (const auto &nbr : ranks) {
    std::vector<Particle> particles = particleNeighborMap[nbr];
    parSendCounts.push_back(static_cast<int32_t>(particles.size()));
    for (auto par : particles) {
      parSendBuffer.push_back(par); // Pack particle by value.
    }
  }

  // Compute send displacements.
  parSDispls.resize(parSendCounts.size());
  if (!parSendCounts.empty()) {
    parSDispls[0] = 0;
    for (uint32_t i = 1; i < parSendCounts.size(); i++) {
      parSDispls[i] = parSDispls[i - 1] + parSendCounts[i - 1];
    }
  }

  // Exchange counts with neighbors.
  int32_t indegree, outdegree, weighted;
  MPI_Dist_graph_neighbors_count(graphComm, &indegree, &outdegree, &weighted);
  int32_t degree = indegree;

  parRecvCounts.resize(degree);
  parRDispls.resize(degree);
  MPI_Neighbor_alltoall(parSendCounts.data(), 1, MPI_INT, parRecvCounts.data(),
                        1, MPI_INT, graphComm);

  uint32_t totalRecv = 0;
  for (uint32_t count : parRecvCounts) {
    totalRecv += count;
  }
  if (degree > 0) {
    parRDispls[0] = 0;
    for (int32_t i = 1; i < degree; i++) {
      parRDispls[i] = parRDispls[i - 1] + parRecvCounts[i - 1];
    }
  }

  // Resize receive buffer (if needed, but if totalRecv is always below a
  // threshold, you could also pre-reserve).
  parRecvBuffer.resize(totalRecv);

  // Perform the MPI exchange.
  int32_t err = MPI_Ineighbor_alltoallv(
      parSendBuffer.data(), parSendCounts.data(), parSDispls.data(),
      particleType, parRecvBuffer.data(), parRecvCounts.data(),
      parRDispls.data(), particleType, graphComm, &requestPar);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Ineighbor_alltoallv failed." << std::endl;
    MPI_Abort(worldComm, err);
  }
}

void Parsim::gatherParticle0() {
  // 1) Local search for particle 0
  Particle localParticle0; // idx == -1 by default
  for (int32_t i = 1; i <= sub.localRows; ++i) {
    for (int32_t j = 1; j <= sub.localCols; ++j) {
      Cell &cell = grid.at(i, j);
      for (Particle *p : cell.particles) {
        if (p->idx == 0) {
          localParticle0 = *p;
          break;
        }
      }
      if (localParticle0.idx == 0)
        break;
    }
    if (localParticle0.idx == 0)
      break;
  }

  // 3) Point-to-point send/receive
  if (worldRank != 0) {
    // Non-root ranks send only if they found it
    if (localParticle0.idx == 0) {
      int32_t err = MPI_Send(&localParticle0, // buffer
                             1,               // count
                             particleType,    // MPI_Datatype
                             0,               // dest = rank 0
                             TAG_PART0,       // tag
                             worldComm);
      if (err != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Send failed sending particle 0 from rank "
                  << worldRank << std::endl;
        MPI_Abort(worldComm, err);
      }
    }
  } else {
    // On rank 0, either keep local copy or receive from whoever has it
    if (localParticle0.idx == 0) {
      // rank 0 itself found it
      globalParticle0 = localParticle0;
    } else {
      MPI_Status status;
      int32_t err =
          MPI_Recv(&globalParticle0, // receive into globalParticle0
                   1,                // count
                   particleType,     // MPI_Datatype
                   MPI_ANY_SOURCE,   // wildcard: accept from any rank
                   MPI_ANY_TAG, // wildcard: accept any tag (or use TAG_PART0)
                   worldComm, &status);
      if (err != MPI_SUCCESS) {
        std::cerr << "Error: MPI_Recv failed on rank 0 waiting for particle 0"
                  << std::endl;
        MPI_Abort(worldComm, err);
      }
    }
  }
}

int32_t Parsim::reduceCollisionCount(int32_t localCollision) {
  uint32_t globalCollision = 0;
  int32_t err = MPI_Reduce(&localCollision, &globalCollision, 1, MPI_INT,
                           MPI_SUM, 0, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Reduce failed in reduceCollisionCount."
              << std::endl;
    MPI_Abort(worldComm, err);
  }
  return globalCollision;
}

// 4: Util functions

void Parsim::extractAndUpdateGhostRegions() {
  MPI_Wait(&requestCOM, MPI_STATUS_IGNORE);
  int32_t n = 0, e = 0, s = 0, w = 0;

  // Create a local mapping of ghost directions to assignment lambdas.
  const std::unordered_map<GhostDirection,
                           std::function<void(const CenterOfMass &)>>
      assignCell = {{GhostDirection::N,
                     [this, &n](const CenterOfMass &cm) {
                       grid.at(0, n + 1).cm = cm;
                       n++;
                     }},
                    {GhostDirection::S,
                     [this, &s](const CenterOfMass &cm) {
                       grid.at(sub.localRows + 1, s + 1).cm = cm;
                       s++;
                     }},
                    {GhostDirection::E,
                     [this, &e](const CenterOfMass &cm) {
                       grid.at(e + 1, sub.localCols + 1).cm = cm;
                       e++;
                     }},
                    {GhostDirection::W,
                     [this, &w](const CenterOfMass &cm) {
                       grid.at(w + 1, 0).cm = cm;
                       w++;
                     }},
                    {GhostDirection::NE,
                     [this](const CenterOfMass &cm) {
                       grid.at(0, sub.localCols + 1).cm = cm;
                     }},
                    {GhostDirection::SE,
                     [this](const CenterOfMass &cm) {
                       grid.at(sub.localRows + 1, sub.localCols + 1).cm = cm;
                     }},
                    {GhostDirection::SW,
                     [this](const CenterOfMass &cm) {
                       grid.at(sub.localRows + 1, 0).cm = cm;
                     }},
                    {GhostDirection::NW, [this](const CenterOfMass &cm) {
                       grid.at(0, 0).cm = cm;
                     }}};

  for (const auto &recv : recvComBuf) {
    // Use your complement lookup (which returns a GhostDirection)
    GhostDirection myGhostRegion = complementMap.at(recv.d);

    // Skip diagonal cases if in 2-process checkered mode.
    if (isCheckeredTwoProc && (myGhostRegion == GhostDirection::NE ||
                               myGhostRegion == GhostDirection::NW ||
                               myGhostRegion == GhostDirection::SE ||
                               myGhostRegion == GhostDirection::SW)) {
      continue;
    }

    // Look up and invoke the assignment lambda.
    assignCell.at(myGhostRegion)(recv);
  }

  // Handle the 2-process special case for missing corners.
  if (isCheckeredTwoProc) {
    grid.at(0, 0).cm = grid.at(0, sub.localCols).cm;     // NW corner
    grid.at(0, sub.localCols + 1).cm = grid.at(0, 1).cm; // NE corner
    grid.at(sub.localRows + 1, 0).cm =
        grid.at(sub.localRows + 1, sub.localCols).cm; // SW corner
    grid.at(sub.localRows + 1, sub.localCols + 1).cm =
        grid.at(sub.localRows + 1, 1).cm; // SE corner
  }
}

int32_t Parsim::getOwnerRank(int32_t newRow, int32_t newCol) {
  // Determine the number of rows and columns per process.
  int32_t rowsPerProc = ncside / dims[0];
  int32_t extraRows = ncside % dims[0];
  int32_t colsPerProc = ncside / dims[1];
  int32_t extraCols = ncside % dims[1];

  // Determine the row coordinate (procRow) of the owning process.
  int32_t procRow = 0, sum = 0;
  for (int32_t i = 0; i < dims[0]; i++) {
    int32_t blockSize = rowsPerProc + (i < extraRows ? 1 : 0);
    if (newRow < sum + blockSize) {
      procRow = i;
      break;
    }
    sum += blockSize;
  }

  // Determine the column coordinate (procCol) of the owning process.
  int32_t procCol = 0;
  sum = 0;
  for (int32_t j = 0; j < dims[1]; j++) {
    int32_t blockSize = colsPerProc + (j < extraCols ? 1 : 0);
    if (newCol < sum + blockSize) {
      procCol = j;
      break;
    }
    sum += blockSize;
  }

  // Convert Cartesian coordinates to rank.
  int32_t coords[2] = {procRow, procCol};
  int32_t ownerRank;
  int32_t err = MPI_Cart_rank(cartComm, coords, &ownerRank);
  if (err != MPI_SUCCESS) {
    std::cerr << "Error: MPI_Cart_rank failed in getOwnerRank." << std::endl;
    MPI_Abort(worldComm, err);
  }
  return ownerRank;
}

} // namespace Simulation
