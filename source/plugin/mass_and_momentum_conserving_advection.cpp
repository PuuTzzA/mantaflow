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
#include <functional>

namespace std
{
    template <>
    struct hash<Manta::Vector3D<int>>
    {
        std::size_t operator()(const Manta::Vector3D<int> &v) const noexcept
        {
            std::size_t h1 = std::hash<int>{}(v.x);
            std::size_t h2 = std::hash<int>{}(v.y);
            std::size_t h3 = std::hash<int>{}(v.z);
            // Combine hashes (this is a common technique)
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

namespace Manta
{

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
            return (flags.isFluid(i, j, k)) && !flags.isObstacle(i, j, k);
        }
    }

    inline bool isSampleableFluid(IndexInt i, IndexInt j, IndexInt k, const FlagGrid &flags, MACGridComponent component)
    {
        Vec3i gs = flags.getParent()->getGridSize();
        bool inBounds;
        switch (component)
        {
        case MAC_X:
            inBounds = 1 <= i && i < gs.x && 0 <= j && j < gs.y && 0 <= k && k < gs.z;
            return inBounds && (flags.isFluid(i, j, k) || flags.isFluid(i - 1, j, k) || flags.isOutflow(i, j, k) || flags.isOutflow(i - 1, j, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i - 1, j, k));
        case MAC_Y:
            inBounds = 0 <= i && i < gs.x && 1 <= j && j < gs.y && 0 <= k && k < gs.z;
            return inBounds && (flags.isFluid(i, j, k) || flags.isFluid(i, j - 1, k) || flags.isOutflow(i, j, k) || flags.isOutflow(i, j - 1, k)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j - 1, k));
        case MAC_Z:
            inBounds = 0 <= i && i < gs.x && 0 <= j && j < gs.y && 1 <= k && k < gs.z;
            return inBounds && (flags.isFluid(i, j, k) || flags.isFluid(i, j, k - 1) || flags.isOutflow(i, j, k) || flags.isOutflow(i, j, k - 1)) && !(flags.isObstacle(i, j, k) || flags.isObstacle(i, j, k - 1));
        default:
            inBounds = 0 <= i && i < gs.x && 0 <= j && j < gs.y && 0 <= k && k < gs.z;
            return inBounds && (flags.isFluid(i, j, k) || flags.isOutflow(i, j, k)) && !flags.isObstacle(i, j, k);
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

    // UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION
    // UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION
    // UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION
    // UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION
    // UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION UNIFIED MASS_MOMENTUM_CONSERVING_ADVECTION

    inline Vec3 rungeKutta4(Vec3 pos, Real dt, const MACGrid &vel)
    {
        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        return pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    }

    /// @brief normal trilinear (bilinear) interpolation
    bool getInterpolationStencilWithWeights(std::vector<std::tuple<Vec3i, Real>> &result, Vec3 pos, const FlagGrid &flags, Vec3 offset, MACGridComponent component)
    {
        result.clear();

        pos -= offset;

        const int i = std::floor(pos[0]);
        const int j = std::floor(pos[1]);
        const int k = std::floor(pos[2]);

        const Real fx = pos[0] - i;
        const Real fy = pos[1] - j;
        const Real fz = pos[2] - k;

        Real w000, w100, w010, w110;
        Real w001, w101, w011, w111;

        w000 = w100 = w010 = w110 = 0.;
        w001 = w101 = w011 = w111 = 0.;

        if (isSampleableFluid(i, j, k, flags, component))
        {
            w000 = (1 - fx) * (1 - fy) * (1 - fz);
        }
        if (isSampleableFluid(i + 1, j, k, flags, component))
        {
            w100 = fx * (1 - fy) * (1 - fz);
        }
        if (isSampleableFluid(i, j + 1, k, flags, component))
        {
            w010 = (1 - fx) * fy * (1 - fz);
        }
        if (isSampleableFluid(i + 1, j + 1, k, flags, component))
        {
            w110 = fx * fy * (1 - fz);
        }
        if (flags.is3D())
        {
            if (isSampleableFluid(i, j, k + 1, flags, component))
            {
                w001 = (1 - fx) * (1 - fy) * fz;
            }
            if (isSampleableFluid(i + 1, j, k + 1, flags, component))
            {
                w101 = fx * (1 - fy) * fz;
            }
            if (isSampleableFluid(i, j + 1, k + 1, flags, component))
            {
                w011 = (1 - fx) * fy * fz;
            }
            if (isSampleableFluid(i + 1, j + 1, k + 1, flags, component))
            {
                w111 = fx * fy * fz;
            }
        }

        Real tot = w000 + w010 + w100 + w110;
        if (flags.is3D())
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
        if (flags.is3D())
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

    inline Real cubicInterpolationWeight(Real s, int idx)
    {
        if (idx == -1)
            return -0.333333 * s + 0.5 * s * s - 0.166666 * s * s * s;
        if (idx == 0)
            return 1.0 - s * s + 0.5 * (s * s * s - s);
        if (idx == 1)
            return s + 0.5 * (s * s - s * s * s);
        if (idx == 2)
            return 0.166666 * (s * s * s - s);
        return 0.0;
    }

    /// @brief better interpolation to reduce numerical diffusion alla Catmull-Rom
    bool getInterpolationStencilWithWeighCatmullRom(std::vector<std::tuple<Vec3i, Real>> &result, Vec3 pos, const FlagGrid &flags, Vec3 offset, MACGridComponent component)
    {
        result.clear();
        pos -= offset;

        const int i = std::floor(pos[0]);
        const int j = std::floor(pos[1]);
        const int k = std::floor(pos[2]);

        const Real fx = pos[0] - i;
        const Real fy = pos[1] - j;
        const Real fz = pos[2] - k;

        const Real EPSILON = 1e-5;
        Real totalWeight = 0.;

        // Loop over 4x4x4 neighborhood
        for (int dz = -1; dz <= (flags.is3D() ? 2 : -1); ++dz)
        {
            Real wz = flags.is3D() ? cubicInterpolationWeight(fz, dz) : 1.;
            int zk = flags.is3D() ? k + dz : k;
            for (int dy = -1; dy <= 2; ++dy)
            {
                Real wy = cubicInterpolationWeight(fy, dy);
                int yj = j + dy;
                for (int dx = -1; dx <= 2; ++dx)
                {
                    Real wx = cubicInterpolationWeight(fx, dx);
                    int xi = i + dx;

                    Real w = wx * wy * wz;

                    if (std::abs(w) > EPSILON && isSampleableFluid(xi, yj, zk, flags, component))
                    {
                        result.push_back({Vec3i{xi, yj, zk}, w});
                        totalWeight += w;
                    }
                }
            }
        }

        if (totalWeight < EPSILON)
            return false;

        // Normalize
        for (auto &entry : result)
        {
            std::get<1>(entry) /= totalWeight;
        }

        return true;
    }

    std::vector<std::tuple<Vec3i, Real>> traceBack(Vec3 pos, Real dt, const MACGrid &vel, const FlagGrid &flags, Vec3 offset, MACGridComponent component, bool doFluid, InterpolationType interpolationType)
    {
        pos += offset;
        Vec3 newPos = rungeKutta4(pos, -dt, vel);

        std::vector<std::tuple<Vec3i, Real>> resultVec{};
        resultVec.reserve(4);

        if (interpolationType == TRILIENAR)
        {
            if (getInterpolationStencilWithWeights(resultVec, newPos, flags, offset, component))
            {
                return resultVec;
            }
        }
        else if (interpolationType == CATMULL_ROM)
        {
            if (getInterpolationStencilWithWeighCatmullRom(resultVec, newPos, flags, offset, component))
            {
                return resultVec;
            }
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

        if (component != NONE && doFluid)
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

            if (interpolationType == TRILIENAR)
            {
                if (getInterpolationStencilWithWeights(testResultVec, current, flags, offset, component))
                {
                    resultVec = testResultVec;
                }
                else
                {
                    return resultVec;
                }
            }
            else if (interpolationType == CATMULL_ROM)
            {
                if (getInterpolationStencilWithWeighCatmullRom(testResultVec, current, flags, offset, component))
                {
                    resultVec = testResultVec;
                }
                else
                {
                    return resultVec;
                }
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
    void knAdvectTraditional(const FlagGrid &flags, const MACGrid &vel, const Grid<T> &oldGrid, Grid<T> &newGrid, Vec3 &offset, Real dt, MACGridComponent component, bool doFluid, InterpolationType interpolationType)
    {
        if (isValidFluid(i, j, k, flags, component))
        {
            newGrid(i, j, k) = 0;
            auto neighboursAndWeights = traceBack(Vec3(i, j, k), dt, vel, flags, offset, component, doFluid, interpolationType);

            for (const auto &[n, w] : neighboursAndWeights)
            {
                newGrid(i, j, k) += w * oldGrid(n);
            }
        }
        else
        {
            // newGrid(i, j, k) = 1.;
        }
    }

    inline void insertIntoWeights(Sparse2DMap<Real> &map, Reverse2dMap &rmap, Vec3i cellI, Vec3i cellJ, Real value)
    {
        map[cellI][cellJ] = value;
        rmap[cellJ].insert(cellI);
    }

    inline void addToWeights(Sparse2DMap<Real> &map, Reverse2dMap &rmap, Vec3i cellI, Vec3i cellJ, Real value)
    {
        map[cellI][cellJ] += value;
        rmap[cellJ].insert(cellI);
    }

    inline void scaleWeightBy(Sparse2DMap<Real> &map, Vec3i cellI, Vec3i cellJ, Real factor)
    {
        map[cellI][cellJ] *= factor;
    }

    void recalculateGamma(Grid<Real> &gamma, const Sparse2DMap<Real> &weights)
    {
        Grid<Real> newGrid(gamma.getParent());
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                newGrid(cellJ) += value;
            }
        }
        gamma.swap(newGrid);
    }

    void recalculateBeta(Grid<Real> &beta, const Sparse2DMap<Real> &weights)
    {
        Grid<Real> newGrid(beta.getParent());
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, value] : innerMap)
            {
                newGrid(cellI) += value;
            }
        }
        beta.swap(newGrid);
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectUnified(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, const Grid<Real> *phi, MACGridComponent component)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-6;
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();

        // Step 0: Advect Gamma with the same tratitional sceme
        Grid<Real> newGammaCum(parent);
        knAdvectTraditional(flags_n_plus_one, vel, gammaCumulative, newGammaCum, offset, dt, component, phi, TRILIENAR);
        // knAdvectGammaCum(vel, grid, newGammaCum, dt, gridSize, offset, flags_n_plus_one);
        gammaCumulative.swap(newGammaCum);

        Sparse2DMap<Real> weights;
        Reverse2dMap reverseWeights;
        Grid<Real> beta(parent);
        Grid<Real> gamma(parent);
        Grid<Real> newGrid(parent);

        int bnd = 1;
        // Step 1: Backwards step
        FOR_IJK_BND(grid, bnd)
        {
            if (!isValidFluid(i, j, k, flags_n_plus_one, component))
            {
                continue;
            }

            auto neighboursAndWeights = traceBack(Vec3(i, j, k), dt, vel, flags_n, offset, component, phi, TRILIENAR);

            if (neighboursAndWeights.empty()) // Find the nearest surface point and dump the excess momentum there
            {
                /* auto surfaceNeighboursAndWeights = getClosestSurfacePoint(Vec3(i, j, k), phi, offset, flags_n, component);
                for (const auto &[n, w] : surfaceNeighboursAndWeights)
                {
                    insertIntoWeights(weights, reverseWeights, n, Vec3i(i, j, k), gridSize, w);
                    beta(i, j, k) += w;
                } */
                /* insertIntoWeights(weights, reverseWeights, Vec3i(i, j, k), Vec3i(i, j, k), gridSize, 1.0);
                beta(i, j, k) += 1; */
            }
            else
            {
                for (const auto &[n, w] : neighboursAndWeights)
                {
                    insertIntoWeights(weights, reverseWeights, n, Vec3i(i, j, k), w);
                    beta(n) += w;
                }
            }
        }

        // Step 2: Forwards Step
        FOR_IJK_BND(grid, bnd)
        {
            if (!isValidFluid(i, j, k, flags_n, component))
            {
                continue;
            }
            if (beta(i, j, k) < 1)
            {
                Real amountToDistribute = 1 - beta(i, j, k);

                auto neighboursAndWeights = traceForward(Vec3(i, j, k), dt, vel, flags_n_plus_one, offset, component);

                if (neighboursAndWeights.empty() && phi) // Find the nearest surface point and dump the excess momentum there
                {
                    std::cout << "ASLDkfj" << std::endl;
                    auto surfaceNeighboursAndWeights = getClosestSurfacePoint(Vec3(i, j, k), *phi, offset, flags_n_plus_one, component);
                    if (surfaceNeighboursAndWeights.empty()) // Now really just dump it into the same cell
                    {
                        addToWeights(weights, reverseWeights, Vec3i(i, j, k), Vec3i(i, j, k), amountToDistribute);
                    }
                    else
                    {
                        for (const auto &[n, w] : surfaceNeighboursAndWeights)
                        {
                            addToWeights(weights, reverseWeights, Vec3i(i, j, k), n, w * amountToDistribute);
                        }
                    }
                }
                else
                {
                    for (const auto &[n, w] : neighboursAndWeights)
                    {
                        addToWeights(weights, reverseWeights, Vec3i(i, j, k), n, w * amountToDistribute);
                    }
                }
            }
        }

        // Step 3: clamp gamma
        recalculateGamma(gamma, weights);
        FOR_IJK_BND(grid, bnd)
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

            for (Vec3i cellI : reverseWeights[Vec3i(i, j, k)])
            {
                scaleWeightBy(weights, cellI, Vec3i(i, j, k), factor);
            }
        }

        // Step 4: clamp beta
        recalculateBeta(beta, weights);
        FOR_IJK_BND(grid, bnd)
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
            for (auto &[_, value] : weights[Vec3i(i, j, k)])
            {
                value *= factor;
            }
        }

        // Step 5: intermediate Result
        for (const auto &[cellI, innerMap] : weights)
        {
            for (const auto &[cellJ, weight] : innerMap)
            {
                newGrid(cellJ) += weight * grid(cellI);
            }
        }

        // Step 6: Diffuse Gamma with per-axis sweeps
        recalculateGamma(gamma, weights);
        for (int _ = 0; _ < 10; _++)
        {
            std::array<Vec3i, 3> dirs{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

            for (auto &d : dirs)
            {
                Grid<Real> deltaGrid(parent);
                Grid<Real> deltaGamma(parent);

                FOR_IJK_BND(grid, bnd)
                {
                    if (!isValidFluid(i, j, k, flags_n_plus_one, component) || !isValidFluid(i + d.x, j + d.y, k + d.z, flags_n_plus_one, component))
                    {
                        continue;
                    }

                    Vec3i idx_moved = Vec3i(i + d.x, j + d.y, k + d.z);

                    Real gammaToMove = (gamma(idx_moved) - gamma(i, j, k)) / 2.;
                    Real denominator = gamma(idx_moved);

                    if (abs(denominator) < EPSILON)
                    {
                        continue;
                    }

                    Real fraction_to_move = gammaToMove / denominator;

                    fraction_to_move = Manta::clamp(fraction_to_move, static_cast<Real>(-0.5), static_cast<Real>(0.5));

                    Real phiToMove = newGrid(idx_moved) * fraction_to_move;

                    if (!std::isfinite(phiToMove) || std::isnan(phiToMove))
                    {
                        continue;
                    }

                    gammaToMove = gamma(idx_moved) * fraction_to_move;

                    deltaGamma(idx_moved) -= gammaToMove;
                    deltaGamma(i, j, k) += gammaToMove;
                    deltaGrid(idx_moved) -= phiToMove;
                    deltaGrid(i, j, k) += phiToMove;
                }

                gamma.add(deltaGamma);
                newGrid.add(deltaGrid);
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

    void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const FlagGrid &flags_n_plus_one, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative, bool water, const Grid<Real> *phi)
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
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velX, gammaX, offsetX, nullptr, MAC_X);
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velY, gammaY, offsetY, nullptr, MAC_Y);
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velZ, gammaZ, offsetZ, nullptr, MAC_Z);
        }
        else
        {
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velX, gammaX, offsetX, phi, MAC_X);
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velY, gammaY, offsetY, phi, MAC_Y);
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(parent, flags, flags_n_plus_one, vel, velZ, gammaZ, offsetZ, phi, MAC_Z);
        }

        knGrids2MAC(grid, velX, velY, velZ, flags);
        knGrids2MAC(gammaCumulative, gammaX, gammaY, gammaZ, flags);
    }

    PYTHON()
    void massMomentumConservingAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(flags->getParent(), *flags, *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative), false, nullptr);
        }
        else if (grid->getType() & GridBase::TypeVec3)
        {
            // fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
        }
        else
            errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, MAC)");
    }

    PYTHON()
    void massMomentumConservingAdvectWater(const FlagGrid *flags_n, const FlagGrid *flags_n_plus_one, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative, Grid<Real> *phi)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            fnMassMomentumConservingAdvectUnified<Grid<Real>>(flags_n->getParent(), *flags_n, *flags_n_plus_one, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5), phi, NONE);
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

    // OTHER_ADVECTION_FUNCTIONS OTHER_ADVECTION_FUNCITONS OTHER_ADVECTION_FUNCTIONS
    // OTHER_ADVECTION_FUNCTIONS OTHER_ADVECTION_FUNCITONS OTHER_ADVECTION_FUNCTIONS
    // OTHER_ADVECTION_FUNCTIONS OTHER_ADVECTION_FUNCITONS OTHER_ADVECTION_FUNCTIONS
    // OTHER_ADVECTION_FUNCTIONS OTHER_ADVECTION_FUNCITONS OTHER_ADVECTION_FUNCTIONS
    // OTHER_ADVECTION_FUNCTIONS OTHER_ADVECTION_FUNCITONS OTHER_ADVECTION_FUNCTIONS

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

    void fnSimpleSLAdcetMAC(const FlagGrid &flags, const MACGrid &vel, MACGrid &grid, FluidSolver *parent, Real dt, InterpolationType interpolationType)
    {
        Grid<Real> gridX(parent);
        Grid<Real> gridY(parent);
        Grid<Real> gridZ(parent);

        Grid<Real> newGridX(parent);
        Grid<Real> newGridY(parent);
        Grid<Real> newGridZ(parent);

        knMAC2Grids(grid, gridX, gridY, gridZ);

        /* Vec3 offsetX = Vec3(0.5, 0.0, 0.0);
        Vec3 offsetY = Vec3(0.0, 0.5, 0.0);
        Vec3 offsetZ = Vec3(0.0, 0.0, 0.5); */

        Vec3 offsetX = Vec3(0.0, 0.5, 0.5);
        Vec3 offsetY = Vec3(0.5, 0.0, 0.5);
        Vec3 offsetZ = Vec3(0.5, 0.5, 0.0);

        knAdvectTraditional<Real>(flags, vel, gridX, newGridX, offsetX, dt, MAC_X, false, interpolationType);
        knAdvectTraditional<Real>(flags, vel, gridY, newGridY, offsetY, dt, MAC_Y, false, interpolationType);
        knAdvectTraditional<Real>(flags, vel, gridZ, newGridZ, offsetZ, dt, MAC_Z, false, interpolationType);

        knGrids2MAC(grid, newGridX, newGridY, newGridZ, flags);
    }

    PYTHON()
    void simpleSLAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, int interpolationType)
    {
        Manta::FluidSolver *parent = flags->getParent();
        Real dt = parent->getDt();

        if (grid->getType() & GridBase::TypeReal)
        {
            Grid<Real> newGrid(parent);
            Vec3 offset = Vec3(0.5, 0.5, 0.5);
            knAdvectTraditional<Real>(*flags, *vel, *((Grid<Real> *)grid), newGrid, offset, dt, NONE, false, (InterpolationType)interpolationType);
            ((Grid<Real> *)grid)->swap(newGrid);
        }
        else if (grid->getType() & GridBase::TypeMAC)
        {
            fnSimpleSLAdcetMAC(*flags, *vel, *((MACGrid *)grid), parent, dt, (InterpolationType)interpolationType);
        }
        else
            errMsg("simpleSLAdvect: Grid Type is not supported (only Real, MAC)");
    }
}