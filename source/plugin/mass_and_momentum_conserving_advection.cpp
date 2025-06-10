#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
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
        return (flags.isFluid(i, j, k) || flags.isOutflow(i, j, k) || flags.isInflow(i, j, k)) && i >= 0 && i <= gs[0] - 1 && j >= 0 && j <= gs[1] - 1 && k >= 0 && k <= gs[2] - 1;
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
            else
            {
                break;
            }
        }
        return lastKnownFluidPos;
    }

    inline Vec3 RK4(Vec3 pos, Real dt, const MACGrid &vel)
    {
        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        return pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    }

    Vec3 customTraceWaterBack(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags_n, Vec3i &gs, Vec3 &offset)
    {
        if (flags_n.isObstacle(pos))
        {
            throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 nextPos = RK4(pos, dt, vel);

        if (flags_n.isFluid(nextPos))
        {
            return nextPos;
        }

        if (offset[0] == 0.0) // MAC Grid X-component
        {
            Vec3 nextPos1 = RK4(pos - Vec3(0.5, 0, 0), dt, vel);
            bool firstInside = flags_n.isFluid(nextPos1);

            Vec3 nextPos2 = RK4(pos + Vec3(0.5, 0, 0), dt, vel);
            bool secondInside = flags_n.isFluid(nextPos2);

            if (firstInside && secondInside)
            {
                return (nextPos1 + nextPos2) * 0.5;
            }
            if (firstInside)
            {
                return nextPos1 + Vec3(0.5, 0, 0);
            }
            if (secondInside)
            {
                return nextPos2 - Vec3(0.5, 0, 0);
            }
        }

        if (offset[1] == 0.0) // MAC Grid Y-component
        {
            Vec3 nextPos1 = RK4(pos - Vec3(0, 0.5, 0), dt, vel);
            bool firstInside = flags_n.isFluid(nextPos1);

            Vec3 nextPos2 = RK4(pos + Vec3(0, 0.5, 0), dt, vel);
            bool secondInside = flags_n.isFluid(nextPos2);

            if (firstInside && secondInside)
            {
                return (nextPos1 + nextPos2) * 0.5;
            }
            if (firstInside)
            {
                return nextPos1 + Vec3(0, 0.5, 0);
            }
            if (secondInside)
            {
                return nextPos2 - Vec3(0, 0.5, 0);
            }
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

            if (flags_n.isFluid(currentTestPoint))
            {
                lastKnownFluidPos = currentTestPoint;
            }
        }
        return lastKnownFluidPos;
    }

    std::tuple<Vec3, bool> customTraceWaterForward(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(pos))
        {
            throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 nextPos = RK4(pos, dt, vel);

        if (flags.isFluid(nextPos))
        {
            return {nextPos, false};
        }

        Vec3 segmentStart = pos;
        Vec3 segmentEnd = nextPos;
        Vec3 lastKnownFluidPos = pos;

        Vec3 direction = segmentEnd - segmentStart;
        Real totalDistance = norm(direction);

        if (totalDistance < 1e-9)
        {
            return {nextPos, true};
        }

        int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.25)));
        bool foundFluid = false;
        for (int i = 1; i <= numSearchSteps; ++i)
        {
            Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
            Vec3 currentTestPoint = segmentStart + t * direction;

            if (flags.isFluid(currentTestPoint))
            {
                lastKnownFluidPos = currentTestPoint;
                foundFluid = true;
            }
            else if (foundFluid)
            {
                return {lastKnownFluidPos, true};
            }
        }
        return {lastKnownFluidPos, true};
    }

    std::vector<std::tuple<Vec3i, Real>> getInterpolationstencilAndWeights(const FlagGrid &flags, Vec3 x, Vec3i &gs, Vec3 &offset)
    {
        x -= offset;

        int i = std::floor(x[0]);
        int j = std::floor(x[1]);
        int k = std::floor(x[2]);

        Real fx = x[0] - i;
        Real fy = x[1] - j;
        Real fz = x[2] - k;

        Real w000, w100, w010, w110;
        Real w001, w101, w011, w111;

        w000 = w100 = w010 = w110 = 0.;
        w001 = w101 = w011 = w111 = 0.;

        if (isValid(i, j, k, flags, gs))
        {
            w000 = (1 - fx) * (1 - fy) * (1 - fz);
        }
        if (isValid(i + 1, j, k, flags, gs))
        {
            w100 = fx * (1 - fy) * (1 - fz);
        }
        if (isValid(i, j + 1, k, flags, gs))
        {
            w010 = (1 - fx) * fy * (1 - fz);
        }
        if (isValid(i + 1, j + 1, k, flags, gs))
        {
            w110 = fx * fy * (1 - fz);
        }
        if (gs.z > 1)
        {
            if (isValid(i, j, k + 1, flags, gs))
            {
                w001 = (1 - fx) * (1 - fy) * fz;
            }
            if (isValid(i + 1, j, k + 1, flags, gs))
            {
                w101 = fx * (1 - fy) * fz;
            }
            if (isValid(i, j + 1, k + 1, flags, gs))
            {
                w011 = (1 - fx) * fy * fz;
            }
            if (isValid(i + 1, j + 1, k + 1, flags, gs))
            {
                w111 = fx * fy * fz;
            }
        }

        Real tot = w000 + w010 + w100 + w110;
        if (gs.z > 1)
        {
            tot += w001 + w011 + w101 + w111;
        }

        if (tot < 1e-5)
        {
            return {};
        }

        w000 /= tot;
        w100 /= tot;
        w010 /= tot;
        w110 /= tot;
        if (gs.z > 1)
        {
            w001 /= tot;
            w101 /= tot;
            w011 /= tot;
            w111 /= tot;
        }

        std::vector<std::tuple<Vec3i, Real>> result{};
        result.reserve(gs.z > 1 ? 8 : 4);

        const Real EPSILON = 1e-5;

        if (w000 > EPSILON)
        {
            result.push_back({Vec3i{i, j, 0}, w000});
        }
        if (w100 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j, 0}, w100});
        }
        if (w010 > EPSILON)
        {
            result.push_back({Vec3i{i, j + 1, 0}, w010});
        }
        if (w110 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j + 1, 0}, w110});
        }
        if (gs.z > 1)
        {
            if (w001 > EPSILON)
                result.push_back({Vec3i{i, j, k + 1}, w001});
            if (w101 > EPSILON)
                result.push_back({Vec3i{i + 1, j, k + 1}, w101});
            if (w011 > EPSILON)
                result.push_back({Vec3i{i, j + 1, k + 1}, w011});
            if (w111 > EPSILON)
                result.push_back({Vec3i{i + 1, j + 1, k + 1}, w111});
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
        if (!flags.isFluid(i, j, k))
        {
            // newGrid(i, j, k) = 1;
        }
        else
        {
            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTraceWaterBack(newPos, vel, -dt, flags, gs, offset);

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

                if (flags.isObstacle(i, j, k) || flags.isOutflow(i, j, k))
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
                if (flags.isObstacle(i, j, k) || flags.isOutflow(i, j, k))
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

                if (flags.isObstacle(i, j, k) || flags.isOutflow(i, j, k))
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
            if (flags.isObstacle(cellI / gridSize[1], cellI % gridSize[1], 0) || flags.isOutflow(cellI / gridSize[1], cellI % gridSize[1], 0))
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

    inline bool isValidInOne(Vec3 firstPos, Vec3 secondPos, const FlagGrid &flags, const FlagGrid &flags_2, Vec3i gs)
    {
        if (firstPos[0] < 0 || firstPos[0] >= gs[0] || firstPos[1] < 0 || firstPos[1] >= gs[1] || firstPos[2] < 0 || firstPos[2] >= gs[2])
        {
            return false;
        }
        if (secondPos[0] < 0 || secondPos[0] >= gs[0] || secondPos[1] < 0 || secondPos[1] >= gs[1] || secondPos[2] < 0 || secondPos[2] >= gs[2])
        {
            return false;
        }
        return (flags_2.isFluid(firstPos) || flags_2.isFluid(firstPos)) && (flags_2.isFluid(secondPos) || flags_2.isFluid(secondPos));
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWater(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        std::cout << "Water version" << std::endl;

        // For testing of the "normal" advection step that is used in this function
        /*         Grid<T> testGrid(parent);
                Vec3i testGS = parent->getGridSize();
                advectGammaCumWater<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags_n, testGS);
                grid.swap(testGrid);
                return; */

        // Advect the cummulative Gamma the same way as later the rest
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Grid<Real> newGammaCum(parent);
        advectGammaCumWater<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags_n, gridSize);
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

                if (!(flags_n.isFluid(i, j, k)))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                newPos = customTraceWaterBack(newPos, vel, -dt, flags_n, gridSize, offset);

                auto neighboursAndWeights = getInterpolationstencilAndWeights(flags_n, newPos, gridSize, offset);
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
                if (!(flags_n_plus_one.isFluid(i, j, k)))
                {
                    continue;
                }

                IndexInt cellI = i * gridSize[1] + j;

                if (beta[cellI] < 1 - EPSILON)
                {
                    Vec3 posForward = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                    auto [newPosForward, hitObstacle] = customTraceWaterForward(posForward, vel, dt, flags_n_plus_one, gridSize);
                    posForward = newPosForward;

                    Real amountToDistribute = 1. - beta[cellI];

                    if (!hitObstacle)
                    {
                        auto neighboursAndWeights = getInterpolationstencilAndWeights(flags_n_plus_one, posForward, gridSize, offset);
                        for (const auto &[n, w] : neighboursAndWeights)
                        {
                            IndexInt cellJ = n[0] * gridSize[1] + n[1];

                            weights[cellI][cellJ] += w * amountToDistribute;
                            reverseWeights[cellJ].insert(cellI);
                        }
                    }
                    else
                    {
                        int x = std::floor(posForward[0]);
                        int y = std::floor(posForward[1]);
                        int k = 0;

                        std::vector<IndexInt> candidates;
                        if (flags_n_plus_one.isFluid(x, y, k))
                            candidates.push_back(x * gridSize[1] + y);
                        if (flags_n_plus_one.isFluid(x + 1, y, k))
                            candidates.push_back((x + 1) * gridSize[1] + y);
                        if (flags_n_plus_one.isFluid(x, y + 1, k))
                            candidates.push_back(x * gridSize[1] + (y + 1));
                        if (flags_n_plus_one.isFluid(x + 1, y + 1, k))
                            candidates.push_back((x + 1) * gridSize[1] + (y + 1));

                        if (!candidates.empty())
                        {
                            Real w = 1.0 / candidates.size();
                            for (auto cellJ : candidates)
                            {
                                weights[cellI][cellJ] += w * amountToDistribute;
                                reverseWeights[cellJ].insert(cellI);
                            }
                        }
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

                if (!(flags_n_plus_one.isFluid(i, j, k) || flags_n.isFluid(i, j, k)))
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
            if (!(flags_n.isFluid(cellI / gridSize[1], cellI % gridSize[1], 0) || flags_n_plus_one.isFluid(cellI / gridSize[1], cellI % gridSize[1], 0)))
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

                    if (!isValidInOne(Vec3(x, y, k), Vec3(x + 1, y, k), flags_n, flags_n_plus_one, gridSize))
                    {
                        continue;
                    }

                    if (gamma[cellI] < EPSILON || gamma[cellI_1] < EPSILON){
                        std::cout << "gamma is zero" << std::endl;
                        continue;
                    }

                    Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;

                    gamma[cellI] += fluxGamma;
                    gamma[cellI_1] -= fluxGamma;

                    T gammaToMove = newGrid(x + 1, y, k) * (fluxGamma / gamma[cellI_1]);

                    if (!std::isfinite(gammaToMove))
                    {
                        std::cerr << "gammaToMove is not finite! "
                                  << "fluxGamma=" << fluxGamma << ", "
                                  << "gamma[cellI_1]=" << gamma[cellI_1] << "\n";
                        std::abort();
                    }

                    newGrid(x, y, k) += gammaToMove;
                    newGrid(x + 1, y, k) -= gammaToMove;
                }
            }

            /*             // Y-Dimension
                        for (IndexInt x = bnd; x < gridSize[0] - bnd; x++)
                        {
                            for (IndexInt y = bnd; y < gridSize[1] - bnd - 1; y++)
                            {
                                int k = 0;
                                IndexInt cellI = x * gridSize[1] + y;
                                IndexInt cellI_1 = cellI + 1;

                                if (!isValidInOne(Vec3(x, y, k), Vec3(x, y + 1, k), flags_n, flags_n_plus_one, gridSize))
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
                        } */
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

    void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const FlagGrid &flags_n_plus_one, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative, bool water)
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
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velX, gammaX, offsetX);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velY, gammaY, offsetY);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velZ, gammaZ, offsetZ);
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
            fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), false);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, MAC)");
    }

    PYTHON()
    void massMomentumConservingAdvectWater(const FlagGrid *flags_n, const FlagGrid *flags_n_plus_one, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvectWater<Grid<Real>>(flags_n->getParent(), *flags_n, *flags_n_plus_one, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags_n->getParent(), *flags_n, *flags_n_plus_one, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), true);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, MAC)");
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

    KERNEL(points)
    void knAdvectParticlesForward(BasicParticleSystem &particles, const MACGrid &vel, Real dt)
    {
        Vec3 pos = particles.getPos(idx);
        particles.setPos(idx, pos + dt * vel.getInterpolatedHi(pos, 2));
    }

    PYTHON()
    void advectParticlesForward(BasicParticleSystem *particles, const MACGrid *vel)
    {
        knAdvectParticlesForward(*particles, *vel, vel->getParent()->getDt());
    }

    // TEST, für mehr performance, geht noch nicht
    KERNEL()
    void backwardsStep(Sparse2DMap<Real> &weights, Reverse2dMap &reverseWeights, Vec3 &offset, const FlagGrid &flags, Real dt, const MACGrid &vel, Vec3i &gridSize)
    {
        if (!flags.isObstacle(i, j, k))
        {
            IndexInt cellJ = i * gridSize[1] + j;

            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTrace(newPos, vel, -dt, flags, gridSize);

            auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
            for (const auto &[n, w] : neighboursAndWeights)
            {
                IndexInt cellI = n[0] * gridSize[1] + n[1];
                weights[cellI][cellJ] = w;
                reverseWeights[cellJ].insert(cellI);
            }
        }
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectFast(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
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
        backwardsStep(weights, reverseWeights, offset, flags, dt, vel, gridSize);
        recalculateBeta(beta, weights);

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
}