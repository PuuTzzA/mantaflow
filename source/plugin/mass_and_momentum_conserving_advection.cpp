#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
#include "mass_and_momentum_conserving_advection.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <numeric>
#include <algorithm>

namespace Manta
{
    template <typename T>
    using Sparse2DMap = std::unordered_map<IndexInt, std::unordered_map<IndexInt, T>>;
    using Reverse2dMap = std::unordered_map<IndexInt, std::unordered_set<IndexInt>>;

    inline IndexInt vecToIdx(Vec3i vec, Vec3i gridSize)
    {
        return vec.x * gridSize.y + vec.y;
    }

    inline Vec3i idxToVec(IndexInt idx, Vec3i gridSize)
    {
        return Vec3i(idx / gridSize.y, idx % gridSize.y, 0);
    }

    inline bool isValidFluid(IndexInt i, IndexInt j, IndexInt k, const FlagGrid &flags, MACGridComponent component)
    {
        switch (component)
        {
        case MAC_X:
            return (flags.isFluid(i, j, k) || flags.isFluid(i - 1, j, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i - 1, j, k));
        case MAC_Y:
            return (flags.isFluid(i, j, k) || flags.isFluid(i, j - 1, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j - 1, k));
        case MAC_Z:
            return (flags.isFluid(i, j, k) || flags.isFluid(i, j, k - 1)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j, k - 1));
        default:
            return flags.isFluid(i, j, k);
        }
    }

    //! Semi-Lagrange interpolation kernel
    KERNEL()
    void knFillHelper(Grid<Real> &dst)
    {
        dst(i, j, k) = 1;
    }

    KERNEL()
    void knFillHelperMAC(MACGrid &dst)
    {
        dst(i, j, k) = Vec3(1., 1., 1.);
    }

    PYTHON()
    void fillWithOnes(GridBase *grid)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            knFillHelper(*((Grid<Real> *)grid)).run();
            std::cout << "filled Grid<Real> with ones" << std::endl;
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            knFillHelperMAC(*((MACGrid *)grid)).run();
            std::cout << "filled MACGrid with ones" << std::endl;
        }
    }

    KERNEL()
    void knFluidFillHelper(Grid<Real> &grid, const FlagGrid &flags)
    {
        if (flags.isFluid(i, j, k))
        {
            grid(i, j, k) = 1.;
        }
        else
        {
            grid(i, j, k) = 0.;
        }
    }

    PYTHON()
    void fillFluidWithOnes(GridBase *grid, const FlagGrid *flags)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            knFluidFillHelper(*((Grid<Real> *)grid), *flags);
            std::cout << "filled fluid cells in Grid<Real> with ones" << std::endl;
        }
        else
        {
            throw std::runtime_error("fill fluid with ones not implemented for this grid type!");
        }
    }

    // MASS_MOMENTUM_CONSERVING_ADVECTION MASS_MOMENTUM_CONSERVING_ADVECTION
    // MASS_MOMENTUM_CONSERVING_ADVECTION MASS_MOMENTUM_CONSERVING_ADVECTION
    // MASS_MOMENTUM_CONSERVING_ADVECTION MASS_MOMENTUM_CONSERVING_ADVECTION
    // MASS_MOMENTUM_CONSERVING_ADVECTION MASS_MOMENTUM_CONSERVING_ADVECTION
    // MASS_MOMENTUM_CONSERVING_ADVECTION MASS_MOMENTUM_CONSERVING_ADVECTION

    bool isValid(int i, int j, int k, const FlagGrid &flags, Vec3i &gs)
    {
        return (!flags.isObstacle(i, j, k)) && i >= 0 && i <= gs[0] - 1 && j >= 0 && j <= gs[1] - 1 && k >= 0 && k <= gs[2] - 1;
    }

    bool isValid(Vec3 pos, const FlagGrid &flags, Vec3i &gs)
    {
        return isValid(std::floor(pos.x), std::floor(pos.y), std::floor(pos.z), flags, gs);
    }

    inline bool isFluid(IndexInt i, IndexInt j, IndexInt k, const FlagGrid &flags, MACGridComponent component)
    {
        switch (component)
        {
        case MAC_X:
            return (flags.isFluid(i, j, k) || flags.isFluid(i - 1, j, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i - 1, j, k));
        case MAC_Y:
            return (flags.isFluid(i, j, k) || flags.isFluid(i, j - 1, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j - 1, k));
        case MAC_Z:
            return (flags.isFluid(i, j, k) || flags.isFluid(i, j, k - 1)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j, k - 1));
        default:
            return flags.isFluid(i, j, k);
        }
    }

    bool isValidWater(int i, int j, int k, const FlagGrid &flags, Vec3i &gs, MACGridComponent component)
    {
        return isFluid(i, j, k, flags, component) && i >= 0 && i <= gs[0] - 1 && j >= 0 && j <= gs[1] - 1 && k >= 0 && k <= gs[2] - 1;
    }

    inline Vec3 RK4(Vec3 pos, Real dt, const MACGrid &vel)
    {
        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        return pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    }

    Vec3 customTrace(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs)
    {
        if (flags.isObstacle(pos))
        {
            // throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 nextPos = RK4(pos, dt, vel);

        if (isValid(nextPos, flags, gs))
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

            if (isValid(currentTestPoint, flags, gs))
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

    Real getDx(const Grid<Real> &grid)
    {
        return grid.getParent()->getDx();
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
    void knAdvectGammaCum(const MACGrid &vel, Grid<T> &grid, Grid<T> &newGrid, float dt, Vec3i gridSize, Vec3 &offset, const FlagGrid &flags)
    {
        if (flags.isObstacle(i, j, k))
        {
            // newGrid(i, j, k) = 1;
        }
        else
        {
            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTrace(newPos, vel, -dt, flags, gridSize);

            auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
            for (const auto &[n, w] : neighboursAndWeights)
            {
                newGrid(i, j, k) += w * grid(n);
            }
        }
    }

    KERNEL()
    template <class T>
    void knSetNewGammaCum(Grid<Real> &gammaCum, std::vector<Real> gamma, Vec3i gridSize)
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
    void fnMassMomentumConservingAdvect(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, MACGridComponent component)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        // For testing of the "normal" advection step that is used in this function
        /* Grid<T> testGrid(parent);
        knAdvectGammaCum<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags);
        grid.swap(testGrid);
        return; */

        // Advect the cummulative Gamma the same way as later the rest
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Grid<Real> newGammaCum(parent);
        knAdvectGammaCum<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags);
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
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags, component))
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

        // Step 2: forward step for all beta < 1
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags, component))
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

        // Step 3: Clamp gamma to the cumulative gamma
        recalculateGamma(gamma, weights);
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags, component))
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

        // Step 4: Clamp beta to 1 for conservation
        recalculateBeta(beta, weights);
        for (IndexInt cellI = 0; cellI < numCells; cellI++)
        {
            int i = cellI / gridSize[1];
            int j = cellI % gridSize[1];
            int k = 0;
            if (!isValidFluid(i, j, k, flags, component))
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
        for (int _ = 0; _ < 7; _++)
        {
            std::array<Vec3i, 4> dirs{{{1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}}};
            GridType tempGrid(parent);
            std::vector<Real> tempGamma(numCells, 0.);

            for (auto &d : dirs)
            {
                FOR_IJK(grid)
                {
                    if (!isValidFluid(i, j, k, flags, component) || !isValidFluid(i + d.x, j + d.y, k + d.z, flags, component))
                    {
                        continue;
                    }

                    Vec3i vecMoved = Vec3i(i + d.x, j + d.y, k + d.z);
                    IndexInt indexHere = vecToIdx(Vec3i(i, j, k), gridSize);
                    IndexInt indexMoved = vecToIdx(vecMoved, gridSize);

                    Real averageGamma = (gamma[indexHere] + gamma[indexMoved]) / 2.;
                    Real denominator = 2. * gamma[indexMoved];
                    if (abs(denominator) < EPSILON)
                    {
                        continue;
                    }
                    Real phiToMove = newGrid(vecMoved) * (gamma[indexMoved] - gamma[indexHere]) / denominator;

                    if (!std::isfinite(phiToMove) || std::isnan(phiToMove))
                    {
                        continue;
                    }

                    tempGamma[indexMoved] = tempGamma[indexHere] = averageGamma;
                    tempGrid(vecMoved) = newGrid(vecMoved) - phiToMove;
                    tempGrid(i, j, k) = newGrid(i, j, k) + phiToMove;
                }
                newGrid.copyFrom(tempGrid);
                gamma = tempGamma;
            }
            continue;

            // X-Dimension
            for (IndexInt y = bnd; y < gridSize[1] - bnd; y++)
            {
                for (IndexInt x = bnd; x < gridSize[0] - bnd - 1; x++)
                {
                    int k = 0;

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x + 1, y, k))
                    {
                        continue;
                    }

                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = (x + 1) * gridSize[1] + y;

                    Real gammaAvg = (gamma[cellI_1] - gamma[cellI]) / 2.;
                    T phiToMove = newGrid(x + 1, y, k) * (gammaAvg / gamma[cellI_1]);

                    tempGamma[cellI_1] = gamma[cellI_1] - gammaAvg;
                    tempGamma[cellI] = gamma[cellI] + gammaAvg;

                    tempGrid(x + 1, y, k) = newGrid(x + 1, y, k) - phiToMove;
                    tempGrid(x, y, k) = newGrid(x, y, k) + phiToMove;
                }
            }
            newGrid.copyFrom(tempGrid);
            gamma = tempGamma;

            // Y-Dimension
            for (IndexInt x = bnd; x < gridSize[0] - bnd; x++)
            {
                for (IndexInt y = bnd; y < gridSize[1] - bnd - 1; y++)
                {
                    int k = 0;
                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x, y + 1, k))
                    {
                        continue;
                    }

                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = x * gridSize[1] + y + 1;

                    Real gammaAvg = (gamma[cellI_1] - gamma[cellI]) / 2.;
                    T phiToMove = newGrid(x, y + 1, k) * (gammaAvg / gamma[cellI_1]);

                    tempGamma[cellI_1] = gamma[cellI_1] - gammaAvg;
                    tempGamma[cellI] = gamma[cellI] + gammaAvg;

                    tempGrid(x, y + 1, k) = newGrid(x, y + 1, k) - phiToMove;
                    tempGrid(x, y, k) = newGrid(x, y, k) + phiToMove;
                }
            }
            newGrid.copyFrom(tempGrid);
            gamma = tempGamma;
        }

        knSetNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }

    // WATER WATER WATER WATER WATER WATER WATER WATER WATER
    // WATER WATER WATER WATER WATER WATER WATER WATER WATER
    // WATER WATER WATER WATER WATER WATER WATER WATER WATER
    // WATER WATER WATER WATER WATER WATER WATER WATER WATER
    // WATER WATER WATER WATER WATER WATER WATER WATER WATER

    inline Vec3 rungeKutta4(Vec3 pos, Real dt, const MACGrid &vel)
    {
        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        return pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    }

    bool getInterpolationStencilWithWeights(std::vector<std::tuple<Vec3i, Real>> &result, Vec3 pos, const FlagGrid &flags, Vec3 offset, MACGridComponent component)
    {
        result.clear();

        pos -= offset;

        int i = std::floor(pos[0]);
        int j = std::floor(pos[1]);
        int k = std::floor(pos[2]);

        Real fx = pos[0] - i;
        Real fy = pos[1] - j;
        Real fz = pos[2] - k;

        Real w000, w100, w010, w110;
        Real w001, w101, w011, w111;

        w000 = w100 = w010 = w110 = 0.;
        w001 = w101 = w011 = w111 = 0.;

        if (isValidFluid(i, j, k, flags, component))
        {
            w000 = (1 - fx) * (1 - fy) * (1 - fz);
        }
        if (isValidFluid(i + 1, j, k, flags, component))
        {
            w100 = fx * (1 - fy) * (1 - fz);
        }
        if (isValidFluid(i, j + 1, k, flags, component))
        {
            w010 = (1 - fx) * fy * (1 - fz);
        }
        if (isValidFluid(i + 1, j + 1, k, flags, component))
        {
            w110 = fx * fy * (1 - fz);
        }
        if (flags.getParent()->getGridSize().z > 1)
        {
            if (isValidFluid(i, j, k + 1, flags, component))
            {
                w001 = (1 - fx) * (1 - fy) * fz;
            }
            if (isValidFluid(i + 1, j, k + 1, flags, component))
            {
                w101 = fx * (1 - fy) * fz;
            }
            if (isValidFluid(i, j + 1, k + 1, flags, component))
            {
                w011 = (1 - fx) * fy * fz;
            }
            if (isValidFluid(i + 1, j + 1, k + 1, flags, component))
            {
                w111 = fx * fy * fz;
            }
        }

        Real tot = w000 + w010 + w100 + w110;
        if (flags.getParent()->getGridSize().z > 1)
        {
            tot += w001 + w011 + w101 + w111;
        }

        if (tot < 1e-5)
        {
            return false;
        }

        w000 /= tot;
        w100 /= tot;
        w010 /= tot;
        w110 /= tot;
        if (flags.getParent()->getGridSize().z > 1)
        {
            w001 /= tot;
            w101 /= tot;
            w011 /= tot;
            w111 /= tot;
        }

        const Real EPSILON = 1e-5;

        if (w000 > EPSILON)
        {
            result.push_back({Vec3i{i, j, k}, w000});
        }
        if (w100 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j, k}, w100});
        }
        if (w010 > EPSILON)
        {
            result.push_back({Vec3i{i, j + 1, k}, w010});
        }
        if (w110 > EPSILON)
        {
            result.push_back({Vec3i{i + 1, j + 1, k}, w110});
        }
        if (flags.getParent()->getGridSize().z > 1)
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

        return true;
    }

    std::vector<std::tuple<Vec3i, Real>> traceBack(Vec3 pos, Real dt, const MACGrid &vel, const FlagGrid &flags, Vec3 offset, MACGridComponent component)
    {
        pos += offset;
        Vec3 newPos = rungeKutta4(pos, -dt, vel);

        std::vector<std::tuple<Vec3i, Real>> resultVec{};
        resultVec.reserve(4);

        if (getInterpolationStencilWithWeights(resultVec, newPos, flags, offset, component))
        {
            return resultVec;
        }

        // Special trace back for MAC grid componets
        Vec3 neighbourOffset;
        switch (component)
        {
        case MAC_X:
            neighbourOffset = Vec3(0.5, 0.0, 0.0);
            break;
        case MAC_Y:
            neighbourOffset = Vec3(0.0, 0.5, 0.0);
            break;
        case MAC_Z:
            neighbourOffset = Vec3(0.0, 0.0, 0.5);
            break;
        default:
            break;
        }

        if (component != NONE)
        {
            std::vector<std::tuple<Vec3i, Real>> resultVec1{};
            std::vector<std::tuple<Vec3i, Real>> resultVec2{};

            Vec3 start1 = pos - neighbourOffset;
            Vec3 start2 = pos + neighbourOffset;

            if (isValidFluid(std::floor(start1.x), std::floor(start1.y), std::floor(start1.z), flags, NONE))
            {
                newPos = rungeKutta4(start1, -dt, vel);
                getInterpolationStencilWithWeights(resultVec1, newPos + neighbourOffset, flags, offset, component);
            }

            if (isValidFluid(std::floor(start2.x), std::floor(start2.y), std::floor(start2.z), flags, NONE))
            {
                newPos = rungeKutta4(start2, -dt, vel);
                getInterpolationStencilWithWeights(resultVec2, newPos - neighbourOffset, flags, offset, component);
            }

            resultVec.clear();
            resultVec.reserve(resultVec1.size() + resultVec2.size());
            resultVec.insert(resultVec.end(), resultVec1.begin(), resultVec1.end());
            resultVec.insert(resultVec.end(), resultVec2.begin(), resultVec2.end());

            Real sum = (resultVec1.size() > 0 ? 1. : 0.) + (resultVec2.size() > 0 ? 1. : 0.);
            std::transform(resultVec.begin(), resultVec.end(), resultVec.begin(), [sum](const std::tuple<Vec3i, Real> t)
                           { return std::make_tuple(std::get<0>(t), std::get<1>(t) / sum); });

            return resultVec;
        }

        // Fallback, try finding the closest surface point
        std::vector<std::tuple<Vec3i, Real>> testResultVec{};
        testResultVec.reserve(4);

        Vec3 current = pos;
        Vec3 direction = newPos - pos;
        Real totalDistance = norm(direction);

        if (totalDistance < 1e-9)
        {
            return {};
        }

        int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.25)));
        for (int i = 0; i <= numSearchSteps; ++i)
        {
            Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
            current = pos + t * direction;

            if (getInterpolationStencilWithWeights(testResultVec, current, flags, offset, component))
            {
                resultVec = testResultVec;
            }
            else
            {
                return resultVec;
            }
        }
        return {};
    }

    std::vector<std::tuple<Vec3i, Real>> traceForward(Vec3 pos, Real dt, const MACGrid &vel, const FlagGrid &flags, Vec3 offset, MACGridComponent component)
    {
        pos += offset;
        Vec3 newPos = rungeKutta4(pos, dt, vel);

        std::vector<std::tuple<Vec3i, Real>> resultVec{};
        resultVec.reserve(4);

        if (getInterpolationStencilWithWeights(resultVec, newPos, flags, offset, component))
        {
            return resultVec;
        }

        // Fallback, try finding the closest surface point
        std::vector<std::tuple<Vec3i, Real>> testResultVec{};
        testResultVec.reserve(4);

        Vec3 current = pos;
        Vec3 direction = newPos - pos;
        Real totalDistance = norm(direction);

        if (totalDistance < 1e-9)
        {
            return {};
        }

        int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.25)));
        for (int i = 0; i <= numSearchSteps; ++i)
        {
            Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
            current = pos + t * direction;

            if (getInterpolationStencilWithWeights(testResultVec, current, flags, offset, component))
            {
                resultVec = testResultVec;
            }
            else
            {
                return resultVec;
            }
        }
        return {};
    }

    std::vector<std::tuple<Vec3i, Real>> getClosestSurfacePoint(Vec3 originalPos, const Grid<Real> &phi, Vec3 offset, const FlagGrid &flags, MACGridComponent component)
    {
        Vec3 startingPoint = originalPos + offset;
        Vec3 closestSurfacePoint = startingPoint;

        const int projectionSteps = 5;
        for (int step = 0; step < projectionSteps; ++step)
        {
            Real phiVal = phi.getInterpolatedHi(closestSurfacePoint, 2);
            Vec3 grad = getGradient(phi, closestSurfacePoint.x, closestSurfacePoint.y, closestSurfacePoint.z);

            if (normSquare(grad) < 1e-12)
            {
                break;
            }
            normalize(grad);

            closestSurfacePoint -= grad * phiVal;
        }
        std::vector<std::tuple<Vec3i, Real>> surfaceNeighboursAndWeights;
        getInterpolationStencilWithWeights(surfaceNeighboursAndWeights, closestSurfacePoint, flags, offset, component);

        return surfaceNeighboursAndWeights;
    }

    KERNEL()
    template <class T>
    void knAdvectTraditional(const FlagGrid &flags, const MACGrid &vel, const Grid<T> &oldGrid, Grid<T> &newGrid, Vec3 &offset, Real dt, MACGridComponent component)
    {
        if (isValidFluid(i, j, k, flags, component))
        {
            newGrid(i, j, k) = 0;
            auto neighboursAndWeights = traceBack(Vec3(i, j, k), dt, vel, flags, offset, component);

            for (const auto &[n, w] : neighboursAndWeights)
            {
                newGrid(i, j, k) += w * oldGrid(n);
            }
        }
        else
        {
            newGrid(i, j, k) = 1.;
        }
        newGrid(i, j, k) = 1.;
    }

    inline void insertIntoWeights(Sparse2DMap<Real> &map, Reverse2dMap &rmap, Vec3i cellI, Vec3i cellJ, Vec3i gridSize, Real value)
    {
        IndexInt indexCellI = vecToIdx(cellI, gridSize);
        IndexInt indexCellJ = vecToIdx(cellJ, gridSize);

        map[indexCellI][indexCellJ] = value;
        rmap[indexCellJ].insert(indexCellI);
    }

    inline void addToWeights(Sparse2DMap<Real> &map, Reverse2dMap &rmap, Vec3i cellI, Vec3i cellJ, Vec3i gridSize, Real value)
    {
        IndexInt indexCellI = vecToIdx(cellI, gridSize);
        IndexInt indexCellJ = vecToIdx(cellJ, gridSize);

        map[indexCellI][indexCellJ] += value;
        rmap[indexCellJ].insert(indexCellI);
    }

    inline void scaleWeightBy(Sparse2DMap<Real> &map, Vec3i cellI, Vec3i cellJ, Real factor, Vec3i gridSize)
    {
        IndexInt indexCellI = vecToIdx(cellI, gridSize);
        IndexInt indexCellJ = vecToIdx(cellJ, gridSize);

        map[indexCellI][indexCellJ] *= factor;
    }

    void recalculateGamma(Grid<Real> &gamma, const Sparse2DMap<Real> &weights, Vec3i gridSize)
    {
        Grid<Real> newGrid(gamma.getParent());
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                newGrid(idxToVec(cellJ, gridSize)) += value;
            }
        }
        gamma.swap(newGrid);
    }

    void recalculateBeta(Grid<Real> &beta, const Sparse2DMap<Real> &weights, Vec3i gridSize)
    {
        Grid<Real> newGrid(beta.getParent());
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                newGrid(idxToVec(cellI, gridSize)) += value;
            }
        }
        beta.swap(newGrid);
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWater(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, const Grid<Real> &phi, MACGridComponent component = NONE)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();

        // Step 0: Advect Gamma with the same tratitional sceme
        Grid<Real> newGammaCum(parent);
        knAdvectTraditional(flags_n_plus_one, vel, gammaCumulative, newGammaCum, offset, dt, component);
        // knAdvectGammaCum(vel, grid, newGammaCum, dt, gridSize, offset, flags_n_plus_one);
        gammaCumulative.swap(newGammaCum);

        Sparse2DMap<Real> weights;
        Reverse2dMap reverseWeights;
        Grid<Real> beta(parent);
        Grid<Real> gamma(parent);
        Grid<Real> newGrid(parent);

        // Step 1: Backwards step
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags_n_plus_one, component))
            {
                continue;
            }

            auto neighboursAndWeights = traceBack(Vec3(i, j, k), dt, vel, flags_n, offset, component);

            if (neighboursAndWeights.empty()) // Find the nearest surface point and dump the excess momentum there
            {
                /* auto surfaceNeighboursAndWeights = getClosestSurfacePoint(Vec3(i, j, k), phi, offset, flags_n, component);
                for (const auto &[n, w] : surfaceNeighboursAndWeights)
                {
                    insertIntoWeights(weights, reverseWeights, n, Vec3i(i, j, k), gridSize, w);
                    beta(i, j, k) += w;
                } */
            }
            else
            {
                for (const auto &[n, w] : neighboursAndWeights)
                {
                    insertIntoWeights(weights, reverseWeights, n, Vec3i(i, j, k), gridSize, w);
                    beta(i, j, k) += w;
                }
            }
        }

        // Step 2: Forwards Step
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags_n, component))
            {
                continue;
            }

            if (beta(i, j, k) < 1 - EPSILON)
            {
                Real amountToDistribute = 1 - beta(i, j, k);

                auto neighboursAndWeights = traceForward(Vec3(i, j, k), dt, vel, flags_n_plus_one, offset, component);

                if (neighboursAndWeights.empty()) // Find the nearest surface point and dump the excess momentum there
                {
                    auto surfaceNeighboursAndWeights = getClosestSurfacePoint(Vec3(i, j, k), phi, offset, flags_n_plus_one, component);
                    if (surfaceNeighboursAndWeights.empty()) // Now really just dump it into the same cell
                    {
                        addToWeights(weights, reverseWeights, Vec3i(i, j, k), Vec3i(i, j, k), gridSize, amountToDistribute);
                    }
                    else
                    {
                        for (const auto &[n, w] : surfaceNeighboursAndWeights)
                        {
                            addToWeights(weights, reverseWeights, Vec3i(i, j, k), n, gridSize, w * amountToDistribute);
                        }
                    }
                }
                else
                {
                    for (const auto &[n, w] : neighboursAndWeights)
                    {
                        addToWeights(weights, reverseWeights, Vec3i(i, j, k), n, gridSize, w * amountToDistribute);
                    }
                }
            }
        }

        // Step 3: clamp gamma
        recalculateGamma(gamma, weights, gridSize);
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags_n_plus_one, component))
            {
                continue;
            }
            if (gamma(i, j, k) < EPSILON)
            {
                continue; // avoid division by 0
            }
            Real factor = gammaCumulative(i, j, k) / gamma(i, j, k);
            // factor = Manta::clamp(factor, static_cast<Real>(0.1), static_cast<Real>(10.));

            for (IndexInt cellI : reverseWeights[vecToIdx(Vec3i(i, j, k), gridSize)])
            {
                scaleWeightBy(weights, idxToVec(cellI, gridSize), Vec3i(i, j, k), factor, gridSize);
            }
        }

        // Step 4: clamp beta
        recalculateBeta(beta, weights, gridSize);
        FOR_IJK(grid)
        {
            if (!isValidFluid(i, j, k, flags_n, component))
            {
                continue;
            }
            if (beta(i, j, k) < EPSILON)
            {
                continue;
            }
            Real factor = 1 / beta(i, j, k);
            for (auto &[_, value] : weights[vecToIdx(Vec3i(i, j, k), gridSize)])
            {
                value *= factor;
            }
        }

        // Step 5: intermediate Result
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, weight] : innerMap)
            {
                newGrid(idxToVec(cellJ, gridSize)) += weight * grid(idxToVec(cellI, gridSize));
            }
        }

        // Step 6: Diffuse Gamma with per-axis sweeps
        recalculateGamma(gamma, weights, gridSize);
        for (int _ = 0; _ < 5; _++)
        {
            std::array<Vec3i, 2> dirs{{{1, 0, 0}, {0, 1, 0}}};
            Grid<Real> tempGrid(parent);
            Grid<Real> tempGamma(parent);

            for (auto &d : dirs)
            {
                FOR_IJK(grid)
                {
                    if (!isValidFluid(i, j, k, flags_n_plus_one, component) || !isValidFluid(i + d.x, j + d.y, k + d.z, flags_n_plus_one, component))
                    {
                        continue;
                    }

                    Vec3i idx_moved = Vec3i(i + d.x, j + d.y, k + d.z);

                    Real averageGamma = (gamma(i, j, k) + gamma(idx_moved)) / 2.;
                    Real denominator = 2. * gamma(idx_moved);
                    if (abs(denominator) < EPSILON)
                    {
                        continue;
                    }
                    Real phiToMove = newGrid(idx_moved) * (gamma(idx_moved) - gamma(i, j, k)) / denominator;

                    if (!std::isfinite(phiToMove) || std::isnan(phiToMove))
                    {
                        continue;
                    }

                    tempGamma(i, j, k) = tempGamma(idx_moved) = averageGamma;
                    tempGrid(idx_moved) = newGrid(idx_moved) - phiToMove;
                    tempGrid(i, j, k) = newGrid(i, j, k) + phiToMove;
                }

                gamma.copyFrom(tempGamma);
                newGrid.copyFrom(tempGrid);
            }
        }

        grid.swap(newGrid);
        gammaCumulative.swap(gamma);
    }

    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON

    void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const FlagGrid &flags_n_plus_one, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative, bool water, const Grid<Real> &phi)
    {
        Grid<Real> velX(parent);
        Grid<Real> velY(parent);
        Grid<Real> velZ(parent);

        Grid<Real> gammaX(parent);
        Grid<Real> gammaY(parent);
        Grid<Real> gammaZ(parent);

        knMAC2Grids(grid, velX, velY, velZ);
        knMAC2Grids(gammaCumulative, gammaX, gammaY, gammaZ);

        /* Vec3 offsetX = Vec3(0.5, 0.0, 0.0);
        Vec3 offsetY = Vec3(0.0, 0.5, 0.0);
        Vec3 offsetZ = Vec3(0.0, 0.0, 0.5); */

        Vec3 offsetX = Vec3(0.0, 0.5, 0.5);
        Vec3 offsetY = Vec3(0.5, 0.0, 0.5);
        Vec3 offsetZ = Vec3(0.5, 0.5, 0.0);

        if (!water)
        {
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velX, gammaX, offsetX, MAC_X);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velY, gammaY, offsetY, MAC_Y);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velZ, gammaZ, offsetZ, MAC_Z);
        }
        else
        {
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velX, gammaX, offsetX, phi, MAC_X);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velY, gammaY, offsetY, phi, MAC_Y);
            fnMassMomentumConservingAdvectWater<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velZ, gammaZ, offsetZ, phi, MAC_Z);
        }

        knGrids2MAC(grid, velX, velY, velZ, flags);
        knGrids2MAC(gammaCumulative, gammaX, gammaY, gammaZ, flags);
    }

    PYTHON()
    void massMomentumConservingAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvect<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5), NONE);
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), false, *((Grid<Real> *)grid));
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, MAC)");
    }

    PYTHON()
    void massMomentumConservingAdvectWater(const FlagGrid *flags_n, const FlagGrid *flags_n_plus_one, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative, Grid<Real> &phi)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvectWater<Grid<Real>>(flags_n->getParent(), *flags_n, *flags_n_plus_one, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5), phi, NONE);
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags_n->getParent(), *flags_n, *flags_n_plus_one, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), true, phi);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, MAC)");
    }

    KERNEL(points)
    void knAdvectParticlesForward(BasicParticleSystem &particles, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i gs)
    {
        Vec3 pos = particles.getPos(idx);
        pos = customTrace(pos, vel, dt, flags, gs);
        particles.setPos(idx, pos);
    }

    PYTHON()
    void advectParticlesForward(BasicParticleSystem *particles, const MACGrid *vel, const FlagGrid *flags)
    {
        knAdvectParticlesForward(*particles, *vel, vel->getParent()->getDt(), *flags, vel->getParent()->getGridSize());
    }

    PYTHON()
    void simpleSLAdvection(const FlagGrid *flags, const MACGrid *vel, Grid<Real> *grid)
    {
        Manta::FluidSolver *parent = flags->getParent();
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Vec3 offset = Vec3(0.5, 0.5, 0.5);

        Grid<Real> newGrid(parent);

        knAdvectGammaCum<Real>(*vel, *grid, newGrid, dt, gridSize, offset, *flags);

        grid->swap(newGrid);
    }
}