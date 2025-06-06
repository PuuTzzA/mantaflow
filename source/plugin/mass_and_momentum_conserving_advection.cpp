#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <array>

namespace Manta
{
    template <typename T>
    using Sparse2DMap = std::unordered_map<IndexInt, std::unordered_map<IndexInt, T>>;
    using Reverse2dMap = std::unordered_map<IndexInt, std::unordered_set<IndexInt>>;

    //! Semi-Lagrange interpolation kernel
    KERNEL()
    void fillHelper(Grid<Real> &dst)
    {
        dst(i, j, k) = 1;
    }

    KERNEL()
    void fillHelperMAC(MACGrid &dst)
    {
        dst(i, j, k) = Vec3(1., 1., 1.);
    }

    PYTHON()
    void fillWithOnes(GridBase *grid)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fillHelper(*((Grid<Real> *)grid)).run();
            std::cout << "filled Grid<Real> with ones" << std::endl;
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fillHelperMAC(*((MACGrid *)grid)).run();
            std::cout << "filled MACGrid with ones" << std::endl;
        }
    }

    bool isValid(int i, int j, int k, const FlagGrid &flags, Vec3i &gs)
    {
        return (flags.isFluid(i, j, k) || flags.isOutflow(i, j, k) || flags.isInflow(i, j, k)) && i >= 0 && i <= gs[0] - 1 && j >= 0 && j <= gs[1] - 1;
    }

    Vec3 customTrace(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(pos))
        {
            throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        Vec3 nextPos = pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);

        if (!flags.isObstacle(nextPos))
        {
            return nextPos;
        }

        Vec3 segmentStart = pos;
        Vec3 segmentEnd = nextPos;
        Vec3 lastKnownFluidPos = pos;

        Vec3 direction = segmentEnd - segmentStart;
        Real totalDistance = norm(direction);

        if (totalDistance < 1e-9)
        {
            return pos;
        }

        int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.25)));
        for (int i = 1; i <= numSearchSteps; ++i)
        {
            Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
            Vec3 currentTestPoint = segmentStart + t * direction;

            if (!flags.isObstacle(currentTestPoint))
            {
                lastKnownFluidPos = currentTestPoint;
            }
        }
        return lastKnownFluidPos;
    }

    Vec3 customTraceWater(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(pos))
        {
            throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        Vec3 nextPos = pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);

        if (!flags.isObstacle(nextPos))
        {
            return nextPos;
        }

        Vec3 segmentStart = pos;
        Vec3 segmentEnd = nextPos;
        Vec3 lastKnownFluidPos = pos;

        Vec3 direction = segmentEnd - segmentStart;
        Real totalDistance = norm(direction);

        if (totalDistance < 1e-9)
        {
            return pos;
        }

        int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.25)));
        for (int i = 1; i <= numSearchSteps; ++i)
        {
            Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
            Vec3 currentTestPoint = segmentStart + t * direction;

            if (!flags.isObstacle(currentTestPoint))
            {
                lastKnownFluidPos = currentTestPoint;
            }
        }
        return lastKnownFluidPos;
    }

    std::vector<std::tuple<Vec3i, Real>> getInterpolationstencilAndWeights(const FlagGrid &flags, Vec3 x, Vec3i &gs, Vec3 &offset)
    {
        x -= offset;

        int i = std::floor(x[0]);
        int j = std::floor(x[1]);

        Real fx = x[0] - i;
        Real fy = x[1] - j;

        Real w00 = 0., w10 = 0., w01 = 0., w11 = 0.;

        if (isValid(i, j, 0, flags, gs))
        {
            w00 = (1 - fx) * (1 - fy);
        }
        if (isValid(i + 1, j, 0, flags, gs))
        {
            w10 = fx * (1 - fy);
        }
        if (isValid(i, j + 1, 0, flags, gs))
        {
            w01 = (1 - fx) * fy;
        }
        if (isValid(i + 1, j + 1, 0, flags, gs))
        {
            w11 = fx * fy;
        }

        Real tot = w00 + w01 + w10 + w11;

        if (tot < 1e-5)
        {
            return {};
        }

        w00 /= tot;
        w10 /= tot;
        w01 /= tot;
        w11 /= tot;

        std::vector<std::tuple<Vec3i, Real>> result{};
        result.reserve(4);

        const Real EPSILON = 1e-5;

        if (w00 > EPSILON)
        {
            result.push_back({Vec3i{i, j, 0}, w00});
        }
        if (w10 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j, 0}, w10});
        }
        if (w01 > EPSILON)
        {
            result.push_back({Vec3i{i, j + 1, 0}, w01});
        }
        if (w11 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j + 1, 0}, w11});
        }

        return result;
    }

    KERNEL(bnd = 0)
    template <class T>
    void advectGammaCum(const MACGrid &vel, Grid<T> &grid, Grid<T> &newGrid, float dt, Vec3i gridSize, Vec3 &offset, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(i, j, k))
        {
            // newGrid(i, j, k) = 1;
        }
        else
        {
            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTrace(newPos, vel, -dt, flags, gs);

            auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
            for (const auto &[n, w] : neighboursAndWeights)
            {
                newGrid(i, j, k) += w * grid(n);
            }
        }
    }

    KERNEL(bnd = 0)
    template <class T>
    void advectGammaCumWater(const MACGrid &vel, Grid<T> &grid, Grid<T> &newGrid, float dt, Vec3i gridSize, Vec3 &offset, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(i, j, k))
        {
            // newGrid(i, j, k) = 1;
        }
        else
        {
            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTraceWater(newPos, vel, -dt, flags, gs);

            auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
            for (const auto &[n, w] : neighboursAndWeights)
            {
                newGrid(i, j, k) += w * grid(n);
            }
        }
    }

    KERNEL()
    template <class T>
    void setNewGammaCum(Grid<Real> &gammaCum, std::vector<Real> gamma, Vec3i gridSize)
    {
        gammaCum(i, j, k) = gamma[i * gridSize[1] + j];
    }

    void recalculateBeta(std::vector<Real> &beta, const Sparse2DMap<Real> &weights)
    {
        std::fill(beta.begin(), beta.end(), 0);

        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                beta[cellI] += value;
            }
        }
    }

    void recalculateGamma(std::vector<Real> &gamma, const Sparse2DMap<Real> &weights)
    {
        std::fill(gamma.begin(), gamma.end(), 0.0);
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                gamma[cellJ] += value;
            }
        }
    }

    template <class GridType>
    void fnMassMomentumConservingAdvect(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        // For testing of the "normal" advection step that is used in this function
        /* Grid<T> testGrid(parent);
        Vec3i testGS = parent->getGridSize();
        advectGammaCum<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags, testGS);
        grid.swap(testGrid);
        return; */

        // Advect the cummulative Gamma the same way as later the rest
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Grid<Real> newGammaCum(parent);
        advectGammaCum<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags, gridSize);
        gammaCumulative.swap(newGammaCum);

        // main advection part
        long unsigned numCells = gridSize[0] * gridSize[1] * gridSize[2];

        GridType newGrid(parent);
        // weights[k][p] = weight from cell k to cell p (cell indeces k/p = i * gridSize[0] + j)
        Sparse2DMap<Real> weights;
        Reverse2dMap reverseWeights;
        std::vector<Real> beta(numCells, 0.);
        std::vector<Real> gamma(numCells, 0.);

        int bnd = 0;

        // Step 1: backwards step
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;

                if (flags.isObstacle(i, j, k))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                newPos = customTrace(newPos, vel, -dt, flags, gridSize);

                auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
                for (const auto &[n, w] : neighboursAndWeights)
                {
                    IndexInt cellI = n[0] * gridSize[1] + n[1];
                    weights[cellI][cellJ] = w;
                    reverseWeights[cellJ].insert(cellI);
                    beta[cellI] += w;
                }
            }
        }

        // Step 2: forward step for all beta < 1
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;
                if (flags.isObstacle(i, j, k))
                {
                    continue;
                }

                IndexInt cellI = i * gridSize[1] + j;

                if (beta[cellI] < 1 - EPSILON)
                {
                    Vec3 posForward = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                    posForward = customTrace(posForward, vel, dt, flags, gridSize);

                    Real amountToDistribute = 1. - beta[cellI];

                    auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, posForward, gridSize, offset);
                    for (const auto &[n, w] : neighboursAndWeights)
                    {
                        IndexInt cellJ = n[0] * gridSize[1] + n[1];

                        weights[cellI][cellJ] += w * amountToDistribute;
                        reverseWeights[cellJ].insert(cellI);
                    }
                }
            }
        }

        // Step 3: Clamp gamma to the cumulative gamma
        recalculateGamma(gamma, weights);
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;

                if (flags.isObstacle(i, j, k))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                if (gamma[cellJ] < EPSILON)
                    continue; // avoid division by 0

                Real factor = gammaCumulative(i, j, k) / gamma[cellJ];

                for (IndexInt cellI : reverseWeights[cellJ])
                {
                    weights[cellI][cellJ] *= factor;
                }
            }
        }

        // Step 4: Clamp beta to 1 for conservation
        recalculateBeta(beta, weights);
        for (IndexInt cellI = 0; cellI < numCells; cellI++)
        {
            if (flags.isObstacle(cellI / gridSize[1], cellI % gridSize[1], 0))
            {
                continue;
            }

            if (beta[cellI] < EPSILON)
                continue; // avoid division by 0

            Real factor = 1 / beta[cellI];

            for (auto &[_, value] : weights[cellI])
            {
                value *= factor;
            }
        }

        // Step 5 calculate the an intermediate result
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, weight] : innerMap)
            {
                int k = 0;

                IndexInt cellI_i = cellI / gridSize[1];
                IndexInt cellI_j = cellI % gridSize[1];

                IndexInt cellJ_i = cellJ / gridSize[1];
                IndexInt cellJ_j = cellJ % gridSize[1];

                newGrid(cellJ_i, cellJ_j, k) += weight * grid(cellI_i, cellI_j, k);
            }
        }

        // Step 6: Diffuse gamma using Gaus seidel Sweep
        recalculateGamma(gamma, weights);
        for (int _ = 0; _ < 5; _++)
        {
            // X-Dimension
            for (IndexInt y = bnd; y < gridSize[1] - bnd; y++)
            {
                for (IndexInt x = bnd; x < gridSize[0] - bnd - 1; x++)
                {
                    int k = 0;
                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = (x + 1) * gridSize[1] + y;

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x + 1, y, k))
                    {
                        continue;
                    }

                    Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;

                    gamma[cellI] += fluxGamma;
                    gamma[cellI_1] -= fluxGamma;

                    T gammaToMove = newGrid(x + 1, y, k) * (fluxGamma / gamma[cellI_1]);
                    newGrid(x, y, k) += gammaToMove;
                    newGrid(x + 1, y, k) -= gammaToMove;
                }
            }

            // Y-Dimension
            for (IndexInt x = bnd; x < gridSize[0] - bnd; x++)
            {
                for (IndexInt y = bnd; y < gridSize[1] - bnd - 1; y++)
                {
                    int k = 0;
                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = cellI + 1;

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x, y + 1, k))
                    {
                        continue;
                    }

                    Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;

                    gamma[cellI] += fluxGamma;
                    gamma[cellI_1] -= fluxGamma;

                    T gammaToMove = newGrid(x, y + 1, k) * (fluxGamma / gamma[cellI_1]);
                    newGrid(x, y, k) += gammaToMove;
                    newGrid(x, y + 1, k) -= gammaToMove;
                }
            }
        }

        setNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWater(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        std::cout << "Water version" << std::endl;

        // For testing of the "normal" advection step that is used in this function
        /* Grid<T> testGrid(parent);
        Vec3i testGS = parent->getGridSize();
        advectGammaCumWater<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags, testGS);
        grid.swap(testGrid);
        return; */

        // Advect the cummulative Gamma the same way as later the rest
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Grid<Real> newGammaCum(parent);
        advectGammaCumWater<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags, gridSize);
        gammaCumulative.swap(newGammaCum);

        // main advection part
        long unsigned numCells = gridSize[0] * gridSize[1] * gridSize[2];

        GridType newGrid(parent);
        // weights[k][p] = weight from cell k to cell p (cell indeces k/p = i * gridSize[0] + j)
        Sparse2DMap<Real> weights;
        Reverse2dMap reverseWeights;
        std::vector<Real> beta(numCells, 0.);
        std::vector<Real> gamma(numCells, 0.);

        int bnd = 0;

        // Step 1: backwards step
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;

                if (!flags.isFluid(i, j, k))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                newPos = customTraceWater(newPos, vel, -dt, flags, gridSize);

                auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
                for (const auto &[n, w] : neighboursAndWeights)
                {
                    IndexInt cellI = n[0] * gridSize[1] + n[1];
                    weights[cellI][cellJ] = w;
                    reverseWeights[cellJ].insert(cellI);
                    beta[cellI] += w;
                }
            }
        }

        // Step 2: forward step for all beta < 1
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;
                if (!flags.isFluid(i, j, k))
                {
                    continue;
                }

                IndexInt cellI = i * gridSize[1] + j;

                if (beta[cellI] < 1 - EPSILON)
                {
                    Vec3 posForward = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                    posForward = customTraceWater(posForward, vel, dt, flags, gridSize);

                    Real amountToDistribute = 1. - beta[cellI];

                    auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, posForward, gridSize, offset);
                    for (const auto &[n, w] : neighboursAndWeights)
                    {
                        IndexInt cellJ = n[0] * gridSize[1] + n[1];

                        weights[cellI][cellJ] += w * amountToDistribute;
                        reverseWeights[cellJ].insert(cellI);
                    }
                }
            }
        }

        // Step 3: Clamp gamma to the cumulative gamma
        recalculateGamma(gamma, weights);
        for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
        {
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                int k = 0;

                if (!flags.isFluid(i, j, k))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                if (gamma[cellJ] < EPSILON)
                    continue; // avoid division by 0

                Real factor = gammaCumulative(i, j, k) / gamma[cellJ];

                for (IndexInt cellI : reverseWeights[cellJ])
                {
                    weights[cellI][cellJ] *= factor;
                }
            }
        }

        // Step 4: Clamp beta to 1 for conservation
        recalculateBeta(beta, weights);
        for (IndexInt cellI = 0; cellI < numCells; cellI++)
        {
            if (!flags.isFluid(cellI / gridSize[1], cellI % gridSize[1], 0))
            {
                continue;
            }

            if (beta[cellI] < EPSILON)
                continue; // avoid division by 0

            Real factor = 1 / beta[cellI];

            for (auto &[_, value] : weights[cellI])
            {
                value *= factor;
            }
        }

        // Step 5 calculate the an intermediate result
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, weight] : innerMap)
            {
                int k = 0;

                IndexInt cellI_i = cellI / gridSize[1];
                IndexInt cellI_j = cellI % gridSize[1];

                IndexInt cellJ_i = cellJ / gridSize[1];
                IndexInt cellJ_j = cellJ % gridSize[1];

                newGrid(cellJ_i, cellJ_j, k) += weight * grid(cellI_i, cellI_j, k);
            }
        }

        // Step 6: Diffuse gamma using Gaus seidel Sweep
        recalculateGamma(gamma, weights);
        for (int _ = 0; _ < 5; _++)
        {
            // X-Dimension
            for (IndexInt y = bnd; y < gridSize[1] - bnd; y++)
            {
                for (IndexInt x = bnd; x < gridSize[0] - bnd - 1; x++)
                {
                    int k = 0;
                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = (x + 1) * gridSize[1] + y;

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x + 1, y, k))
                    {
                        continue;
                    }

                    Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;

                    gamma[cellI] += fluxGamma;
                    gamma[cellI_1] -= fluxGamma;

                    T gammaToMove = newGrid(x + 1, y, k) * (fluxGamma / gamma[cellI_1]);
                    newGrid(x, y, k) += gammaToMove;
                    newGrid(x + 1, y, k) -= gammaToMove;
                }
            }

            // Y-Dimension
            for (IndexInt x = bnd; x < gridSize[0] - bnd; x++)
            {
                for (IndexInt y = bnd; y < gridSize[1] - bnd - 1; y++)
                {
                    int k = 0;
                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = cellI + 1;

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x, y + 1, k))
                    {
                        continue;
                    }

                    Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;

                    gamma[cellI] += fluxGamma;
                    gamma[cellI_1] -= fluxGamma;

                    T gammaToMove = newGrid(x, y + 1, k) * (fluxGamma / gamma[cellI_1]);
                    newGrid(x, y, k) += gammaToMove;
                    newGrid(x, y + 1, k) -= gammaToMove;
                }
            }
        }

        setNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }

    KERNEL()
    void MAC2Grids(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
    {
        Vec3 data = vel(i, j, k);
        velX(i, j, k) = data.x;
        velY(i, j, k) = data.y;
        velZ(i, j, k) = data.z;
    }

    KERNEL()
    void Grids2MAC(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ, const FlagGrid &flags)
    {
        vel(i, j, k) = Vec3(velX(i, j, k), velY(i, j, k), velZ(i, j, k));
    }

    void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative, bool water)
    {
        Grid<Real> velX(parent);
        Grid<Real> velY(parent);
        Grid<Real> velZ(parent);

        Grid<Real> gammaX(parent);
        Grid<Real> gammaY(parent);
        Grid<Real> gammaZ(parent);

        MAC2Grids(grid, velX, velY, velZ);
        MAC2Grids(gammaCumulative, gammaX, gammaY, gammaZ);

        /* Vec3 offsetX = Vec3(0.5, 0.0, 0.0);
        Vec3 offsetY = Vec3(0.0, 0.5, 0.0);
        Vec3 offsetZ = Vec3(0.0, 0.0, 0.5); */

        Vec3 offsetX = Vec3(0.0, 0.5, 0.5);
        Vec3 offsetY = Vec3(0.5, 0.0, 0.5);
        Vec3 offsetZ = Vec3(0.5, 0.5, 0.0);

        if (!water)
        {
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velX, gammaX, offsetX);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velY, gammaY, offsetY);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velZ, gammaZ, offsetZ);
        }
        else
        {
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, vel, velX, gammaX, offsetX);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, vel, velY, gammaY, offsetY);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, vel, velZ, gammaZ, offsetZ);
        }

        Grids2MAC(grid, velX, velY, velZ, flags);
        Grids2MAC(gammaCumulative, gammaX, gammaY, gammaZ, flags);
        return;
    }

    PYTHON()
    void massMomentumConservingAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvect<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), false);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");
    }

    PYTHON()
    void massMomentumConservingAdvectWater(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvectWater<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), true);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");
    }

    PYTHON()
    std::string calculateMass(const Grid<Real> *grid)
    {
        Vec3i gridSize = grid->getParent()->getGridSize();
        int bnd = 1;
        Real sum = 0;
        Real min = std::numeric_limits<Real>::infinity();
        Real max = -std::numeric_limits<Real>::infinity();
        for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
        {
            for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
            {
                Real val = grid->operator()(i, j, 0);
                sum += val;
                min = std::min(min, val);
                max = std::max(max, val);
            }
        }
        return std::to_string(sum) + "," + std::to_string(min) + "," + std::to_string(max);
    }
}