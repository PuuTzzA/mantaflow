#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
#include "mass_and_momentum_conserving_advection.h"
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
            throw std::runtime_error("trace starting from obstacle!");
        }

        Vec3 nextPos = RK4(pos, dt, vel);

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

    Vec3 customTraceWaterBack(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs, Vec3 &offset, MACGridComponent component)
    {
        Vec3 nextPos = RK4(pos, dt, vel);

        if (isFluid(nextPos.x, nextPos.y, nextPos.z, flags, component))
        {
            return nextPos;
        }

        if (component == MAC_X) // MAC Grid X-component
        {
            Vec3 nextPos1 = RK4(pos - Vec3(0.5, 0, 0), dt, vel);
            bool firstInside = flags.isFluid(nextPos1);

            Vec3 nextPos2 = RK4(pos + Vec3(0.5, 0, 0), dt, vel);
            bool secondInside = flags.isFluid(nextPos2);

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

        if (component == MAC_Y) // MAC Grid Y-component
        {
            Vec3 nextPos1 = RK4(pos - Vec3(0, 0.5, 0), dt, vel);
            bool firstInside = flags.isFluid(nextPos1);

            Vec3 nextPos2 = RK4(pos + Vec3(0, 0.5, 0), dt, vel);
            bool secondInside = flags.isFluid(nextPos2);

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

            if (isFluid(currentTestPoint.x, currentTestPoint.y, currentTestPoint.z, flags, component))
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

    Vec3 customTraceWaterForwardOld(Vec3 pos, const MACGrid &vel, const Grid<Real> &phi, Real dt, const FlagGrid &flags, Vec3i &gs, MACGridComponent component)
    {
        if (flags.isObstacle(pos))
        {
            std::cout << "STARTED IN OBSTACLE DIO CAN" << std::endl;
        }

        Vec3 nextPos = RK4(pos, dt, vel);

        nextPos.x = Manta::clamp(nextPos.x, 0.f, static_cast<Real>(gs.x));
        nextPos.y = Manta::clamp(nextPos.y, 0.f, static_cast<Real>(gs.y));
        nextPos.z = Manta::clamp(nextPos.z, 0.f, static_cast<Real>(gs.z));

        // if (isFluid(nextPos.x, nextPos.y, nextPos.z, flags, component))
        if (isFluid(nextPos.x, nextPos.y, nextPos.z, flags, component))
        {
            return nextPos;
        }
        /*         if (flags.isObstacle(nextPos))
                {
                    std::cout << "next pos in obstacle" << std::endl;
                } */

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

            if (isFluid(currentTestPoint.x, currentTestPoint.y, currentTestPoint.z, flags, component))
            // if (flags.isFluid(currentTestPoint))
            {
                // std::cout << "changed last known fluid pos" << std::endl;
                lastKnownFluidPos = currentTestPoint;
            }
            else
            {
                break;
            }
        }
        return lastKnownFluidPos;
    }

    Vec3 getLevelSetGradient(const Grid<Real> &phi, const Vec3 &pos)
    {
        // Standard central differencing for the gradient.
        // Assumes dx=1. If not, divide by dx.
        Real dx = 1;

        Real Gx = (phi.getInterpolated(pos + Vec3(dx, 0, 0)) - phi.getInterpolated(pos - Vec3(dx, 0, 0))) / (2.0 * dx);
        Real Gy = (phi.getInterpolated(pos + Vec3(0, dx, 0)) - phi.getInterpolated(pos - Vec3(0, dx, 0))) / (2.0 * dx);
        Real Gz = (phi.getInterpolated(pos + Vec3(0, 0, dx)) - phi.getInterpolated(pos - Vec3(0, 0, dx))) / (2.0 * dx);

        return Vec3(Gx, Gy, Gz);
    }

    Vec3 customTraceWaterForward(Vec3 pos, const MACGrid &vel, const Grid<Real> &phi, Real dt, const FlagGrid &flags, Vec3i &gs, MACGridComponent component)
    {
        if (flags.isObstacle(pos))
        {
            std::cout << "STARTED IN OBSTACLE DIO CAN" << std::endl;
        }

        Vec3 nextPos = RK4(pos, dt, vel);

        nextPos.x = Manta::clamp(nextPos.x, 0.f, static_cast<Real>(gs.x));
        nextPos.y = Manta::clamp(nextPos.y, 0.f, static_cast<Real>(gs.y));
        nextPos.z = Manta::clamp(nextPos.z, 0.f, static_cast<Real>(gs.z));

        Real phiVal = phi.getInterpolatedHi(nextPos, 2);

        // Check if we are outside. We only need to project if we are.
        if (phi.getInterpolatedHi(nextPos, 2) > 0)
        {
            Vec3 current_pos = nextPos;
            const int max_iter = 5; // 3-5 iterations is usually enough

            for (int i = 0; i < max_iter; ++i)
            {
                Real phi_val = phi.getInterpolated(current_pos);
                Vec3 gradient = getLevelSetGradient(phi, current_pos);

                Real grad_mag_sq = normSquare(gradient);
                if (grad_mag_sq < 1e-8)
                {
                    // Gradient is unstable, can't continue the projection.
                    // Fallback: Use the original position before the trace.
                    current_pos = pos;
                    break;
                }

                // Newton's method step to find the zero-crossing
                current_pos -= gradient * (phi_val / grad_mag_sq);
            }

            // Now `current_pos` is very close to the phi=0 surface.
            // We still need to nudge it inwards.
            Vec3 final_gradient = getLevelSetGradient(phi, current_pos);
            if (normSquare(final_gradient) > 1e-8)
            {
                normalize(final_gradient);
                Real safety_nudge = 0.1;
                nextPos = current_pos - safety_nudge * final_gradient;
            }
            else
            {
                // If even here the gradient is bad, fallback to the pre-trace position.
                nextPos = pos;
            }
        }

        return nextPos;

        // if (isFluid(nextPos.x, nextPos.y, nextPos.z, flags, component))
        if (isFluid(nextPos.x, nextPos.y, nextPos.z, flags, component))
        {
            return nextPos;
        }
        return pos;

        /*         if (flags.isObstacle(nextPos))
                {
                    std::cout << "next pos in obstacle" << std::endl;
                } */

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

            if (isFluid(currentTestPoint.x, currentTestPoint.y, currentTestPoint.z, flags, component))
            // if (flags.isFluid(currentTestPoint))
            {
                // std::cout << "changed last known fluid pos" << std::endl;
                lastKnownFluidPos = currentTestPoint;
            }
            else
            {
                break;
            }
        }
        return lastKnownFluidPos;
    }

    // Helper function to get the grid cell size, crucial for robust projection.
    Real getDx(const Grid<Real> &grid)
    {
        return grid.getParent()->getDx();
    }

    Vec3 customTraceWaterForward_V4(Vec3 pos, const MACGrid &vel, const Grid<Real> &phi, Real dt, Vec3i gs)
    {
        Vec3 nextPos = RK4(pos, dt, vel);

        // Clamp to grid boundaries first to avoid out-of-bounds access.
        // The -1.001f is a common trick to prevent landing exactly on the boundary.
        nextPos.x = Manta::clamp(nextPos.x, 0.f, (Real)gs.x - 1.001f);
        nextPos.y = Manta::clamp(nextPos.y, 0.f, (Real)gs.y - 1.001f);
        nextPos.z = Manta::clamp(nextPos.z, 0.f, (Real)gs.z - 1.001f);

        // Check if we are outside the new fluid domain (phi > 0 means outside)
        if (phi.getInterpolated(nextPos) > 0)
        {
            Vec3 current_pos = nextPos;
            const int max_iter = 5;
            bool projection_succeeded = true;

            for (int i = 0; i < max_iter; ++i)
            {
                Real phi_val = phi.getInterpolated(current_pos);
                if (phi_val <= 0)
                { // We've successfully projected inside
                    break;
                }

                Vec3 gradient = getLevelSetGradient(phi, current_pos);
                Real grad_mag_sq = normSquare(gradient);

                if (grad_mag_sq < 1e-9)
                {
                    // Gradient is unstable, abort the iterative process.
                    projection_succeeded = false;
                    break;
                }

                // --- DAMPENED NEWTON'S STEP ---
                Vec3 correction = gradient * (phi_val / grad_mag_sq);

                // Clamp the magnitude of the correction to prevent huge jumps.
                // A max step of 1.0 grid cell is a reasonable limit.
                if (normSquare(correction) > 1.0)
                {
                    normalize(correction); // Now length is 1.0
                }
                current_pos -= correction;
            }

            if (projection_succeeded && phi.getInterpolated(current_pos) <= 0)
            {
                // Projection succeeded, and we are inside.
                nextPos = current_pos;
            }
            else
            {
                // If projection failed or still outside, we have a problem.
                // Returning the original position is the source of your bug.
                // A better, but not perfect, fallback is to find the closest point
                // on the line segment between pos and nextPos that is inside phi.
                // For now, let's just mark it as invalid to handle it in the main loop.
                // A special value like (-1,-1,-1) can signal failure.
                return Vec3(-1.0);
            }
        }

        // Final clamp for safety before returning.
        nextPos.x = Manta::clamp(nextPos.x, 0.f, (Real)gs.x - 1.001f);
        nextPos.y = Manta::clamp(nextPos.y, 0.f, (Real)gs.y - 1.001f);
        nextPos.z = Manta::clamp(nextPos.z, 0.f, (Real)gs.z - 1.001f);

        return nextPos;
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

    std::vector<std::tuple<Vec3i, Real>> getInterpolationstencilAndWeightsWater(const FlagGrid &flags, Vec3 x, Vec3i &gs, Vec3 &offset, MACGridComponent component)
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

        if (isValidWater(i, j, k, flags, gs, component))
        {
            w000 = (1 - fx) * (1 - fy) * (1 - fz);
        }
        if (isValidWater(i + 1, j, k, flags, gs, component))
        {
            w100 = fx * (1 - fy) * (1 - fz);
        }
        if (isValidWater(i, j + 1, k, flags, gs, component))
        {
            w010 = (1 - fx) * fy * (1 - fz);
        }
        if (isValidWater(i + 1, j + 1, k, flags, gs, component))
        {
            w110 = fx * fy * (1 - fz);
        }
        if (gs.z > 1)
        {
            if (isValidWater(i, j, k + 1, flags, gs, component))
            {
                w001 = (1 - fx) * (1 - fy) * fz;
            }
            if (isValidWater(i + 1, j, k + 1, flags, gs, component))
            {
                w101 = fx * (1 - fy) * fz;
            }
            if (isValidWater(i, j + 1, k + 1, flags, gs, component))
            {
                w011 = (1 - fx) * fy * fz;
            }
            if (isValidWater(i + 1, j + 1, k + 1, flags, gs, component))
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

    KERNEL(bnd = 0)
    template <class T>
    void knAdvectGammaCumWater(const MACGrid &vel, Grid<T> &grid, Grid<T> &newGrid, float dt, Vec3i gridSize, Vec3 &offset, const FlagGrid &flags, const FlagGrid &flags2, MACGridComponent component)
    {
        if (!isFluid(i, j, k, flags2, component))
        {
            // newGrid(i, j, k) = 1;
        }
        else
        {
            Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
            newPos = customTraceWaterBack(newPos, vel, -dt, flags, gridSize, offset, component);

            auto neighboursAndWeights = getInterpolationstencilAndWeightsWater(flags, newPos, gridSize, offset, component);
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
    void fnMassMomentumConservingAdvect(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        // For testing of the "normal" advection step that is used in this function
        /* Grid<T> testGrid(parent);
        advectGammaCum<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags);
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

                    if (!flags.isFluid(x, y, k) || !flags.isFluid(x + 1, y, k))
                    {
                        continue;
                    }

                    IndexInt cellI = x * gridSize[1] + y;
                    IndexInt cellI_1 = (x + 1) * gridSize[1] + y;

                    Real gammaAvg = (gamma[cellI_1] - gamma[cellI]) / 2.;
                    T phiToMove = newGrid(x + 1, y, k) * (gammaAvg / gamma[cellI_1]);

                    gamma[cellI_1] -= gammaAvg;
                    gamma[cellI] += gammaAvg;

                    newGrid(x + 1, y, k) -= phiToMove;
                    newGrid(x, y, k) += phiToMove;
                }
            }

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

                    gamma[cellI_1] -= gammaAvg;
                    gamma[cellI] += gammaAvg;

                    newGrid(x, y + 1, k) -= phiToMove;
                    newGrid(x, y, k) += phiToMove;
                }
            }
        }

        knSetNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }

    inline bool isValidWater(const FlagGrid &flags, const FlagGrid &flags2, int i, int j, int k)
    {
        return flags.isFluid(i, j, k) || flags2.isFluid(i, j, k) || flags.isFluid(i - 1, j, k) || flags2.isFluid(i - 1, j, k) || flags.isFluid(i, j - 1, k) || flags2.isFluid(i, j - 1, k) || flags.isFluid(i + 1, j, k) || flags2.isFluid(i + 1, j, k) || flags.isFluid(i, j - 1, k) || flags2.isFluid(i, j + 1, k);
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWaterOld(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, const Grid<Real> &phi, MACGridComponent component = NONE)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;

        // For testing of the "normal" advection step that is used in this function
        /* Grid<T> testGrid(parent);
        knAdvectGammaCumWater<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags_n, flags_n_plus_one, component);
        grid.swap(testGrid);
        return; */

        // Advect the cummulative Gamma the same way as later the rest
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();
        Grid<Real> newGammaCum(parent);
        knAdvectGammaCumWater<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags_n, flags_n_plus_one, component);
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

                if (!isFluid(i, j, k, flags_n_plus_one, component))
                {
                    continue;
                }

                IndexInt cellJ = i * gridSize[1] + j;

                Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                newPos = customTraceWaterBack(newPos, vel, -dt, flags_n, gridSize, offset, component);

                if (newPos.x == -1. && newPos.y == -1.)
                {
                    continue;
                }

                auto neighboursAndWeights = getInterpolationstencilAndWeightsWater(flags_n, newPos, gridSize, offset, component);

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

                if (!isFluid(i, j, k, flags_n, component))
                {
                    continue;
                }

                IndexInt cellI = i * gridSize[1] + j;

                if (beta[cellI] < 1 - EPSILON)
                {
                    Vec3 startPosForward = Vec3(i + offset[0], j + offset[1], k + offset[2]);
                    // Use the corrected forward tracer
                    Vec3 posForward = customTraceWaterForward_V4(startPosForward, vel, phi, dt, gridSize);

                    // Check if the forward trace failed
                    if (posForward.x < 0.0)
                    { // Using the (-1,-1,-1) failure signal
                        // FALLBACK 1: Trace failed. To conserve mass/momentum,
                        // add the remaining amount back to the original cell.
                        weights[cellI][cellI] += (1. - beta[cellI]);
                        reverseWeights[cellI].insert(cellI); // Ensure reverse map is updated
                        continue;                            // Move to the next cell
                    }

                    Real amountToDistribute = 1. - beta[cellI];
                    auto neighboursAndWeights = getInterpolationstencilAndWeightsWater(flags_n_plus_one, posForward, gridSize, offset, component);

                    if (neighboursAndWeights.empty())
                    {
                        // FALLBACK 2: Stencil is empty, even with a valid trace.
                        // This means the projected point is in an "island" of air cells.
                        // Again, to conserve, add the remaining amount back to the source.
                        weights[cellI][cellI] += amountToDistribute;
                        reverseWeights[cellI].insert(cellI); // Ensure reverse map is updated
                        continue;                            // Move to the next cell
                    }

                    // Original logic for distributing momentum
                    for (const auto &[n, w] : neighboursAndWeights)
                    {
                        IndexInt cellJ = n[0] * gridSize[1] + n[1];
                        weights[cellI][cellJ] += w * amountToDistribute;
                        reverseWeights[cellJ].insert(cellI);
                    }
                }
            }
        }

        for (int __ = 0; __ < 0; __++)
        {
            // Step 3: Clamp gamma to the cumulative gamma
            recalculateGamma(gamma, weights);
            for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
            {
                for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
                {
                    int k = 0;

                    if (!isFluid(i, j, k, flags_n_plus_one, component))
                    {
                        continue;
                    }

                    IndexInt cellJ = i * gridSize[1] + j;

                    if (gamma[cellJ] < EPSILON)
                        continue; // avoid division by 0

                    Real factor = gammaCumulative(i, j, k) / gamma[cellJ];
                    // factor = 1 / gamma[cellJ];

                    if (!std::isfinite(factor))
                    {
                        continue;
                    }
                    factor = Manta::clamp(factor, (Real)0.1, (Real)10.0);

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
                int i = cellI / gridSize[1];
                int j = cellI % gridSize[1];
                int k = 0;

                if (!isFluid(i, j, k, flags_n, component))
                {
                    continue;
                }

                if (beta[cellI] < EPSILON)
                    continue; // avoid division by 0

                Real factor = 1 / beta[cellI];

                if (!std::isfinite(factor))
                    continue;

                for (auto &[_, value] : weights[cellI])
                {
                    value *= factor;
                }
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
            for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
            {
                for (IndexInt i = bnd; i < gridSize[0] - bnd - 1; i++)
                {
                    int k = 0;

                    if (!(isFluid(i, j, k, flags_n_plus_one, component) && isFluid(i + 1, j, k, flags_n_plus_one, component)))
                    {
                        continue;
                    }

                    IndexInt cellI = i * gridSize[1] + j;
                    IndexInt cellI_1 = (i + 1) * gridSize[1] + j;

                    Real gammaAvg = (gamma[cellI_1] - gamma[cellI]) / 2.;
                    T phiToMove = newGrid(i + 1, j, k) * (gammaAvg / gamma[cellI_1]);

                    if (!std::isfinite(phiToMove) || std::isnan(phiToMove))
                    {
                        continue;
                    }

                    gamma[cellI_1] -= gammaAvg;
                    gamma[cellI] += gammaAvg;

                    newGrid(i + 1, j, k) -= phiToMove;
                    newGrid(i, j, k) += phiToMove;
                }
            }

            // Y-Dimension
            for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
            {
                for (IndexInt j = bnd; j < gridSize[1] - bnd - 1; j++)
                {
                    int k = 0;

                    if (!(isFluid(i, j, k, flags_n_plus_one, component) && isFluid(i, j + 1, k, flags_n_plus_one, component)))
                    {
                        continue;
                    }

                    IndexInt cellI = i * gridSize[1] + j;
                    IndexInt cellI_1 = i * gridSize[1] + j + 1;

                    Real gammaAvg = (gamma[cellI_1] - gamma[cellI]) / 2.;
                    T phiToMove = newGrid(i, j + 1, k) * (gammaAvg / gamma[cellI_1]);

                    if (!std::isfinite(phiToMove) || std::isnan(phiToMove))
                    {
                        continue;
                    }

                    gamma[cellI_1] -= gammaAvg;
                    gamma[cellI] += gammaAvg;

                    newGrid(i, j + 1, k) -= phiToMove;
                    newGrid(i, j, k) += phiToMove;
                }
            }
        }

        knSetNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }

    // Completely new Try. gang gang
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

    inline Vec3 rungeKutta4(Vec3 pos, Real dt, const MACGrid &vel)
    {
        Vec3 k1 = vel.getInterpolatedHi(pos, 2);
        Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
        Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
        Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

        return pos + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
    }

    std::vector<std::tuple<Vec3i, Real>> getInterpolationStencilWithWeights(Vec3 pos, const FlagGrid &flags)
    {

    }

    std::vector<std::tuple<Vec3i, Real>> traceBack(Vec3 pos, Real dt, const MACGrid &vel, const FlagGrid &flags, MACGridComponent component)
    {
        pos = rungeKutta4(pos, dt, vel);
    }

    KERNEL()
    template <class T>
    void knAdvectTraditional(const FlagGrid &flags, const MACGrid &vel, const Grid<T> &oldGrid, Grid<T> &newGrid, Vec3 &offset, Real dt, MACGridComponent component)
    {
        if (isValidFluid(i, j, k, flags, component))
        {

            auto neighboursAndWeights =
        }
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWaterOld(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, const Grid<Real> &phi, MACGridComponent component = NONE)
    {
        typedef typename GridType::BASETYPE T;
        const Real EPSILON = 1e-5;
        Real dt = parent->getDt();
        Vec3i gridSize = parent->getGridSize();

        // Step 0: Advect Gamma with the same tratitional sceme
        Grid<Real> newGammaCum(parent);
    }

    // End of completely new Try

    // agealdfjöaslkfdjölaskdjfölaksjdf ölasjkd flkasjldfökjas
    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code

    // Helper function to compute bilinear interpolation weights
    // f: fractional position within the base cell [0,1]^3
    // ii, jj: integer offsets (0 or 1) for the corner
    static inline Real getBilinearWeight(const Vec3 &f, int ii, int jj)
    {
        Real w = 1.0;
        w *= (ii == 0) ? (1.0 - f.x) : f.x;
        w *= (jj == 0) ? (1.0 - f.y) : f.y;
        return w;
    }

    template <class GridType>
    void fnMassMomentumConservingAdvectWater(FluidSolver *parent,
                                             const FlagGrid &flags_n,
                                             const FlagGrid &flags_n_plus_one,
                                             const MACGrid &vel,
                                             GridType &grid,
                                             Grid<Real> &gammaCumulative,
                                             Vec3 offset,
                                             const Grid<Real> &phi, // not strictly used, flags are sufficient
                                             MACGridComponent component)
    {
        // =================================================================================================
        // 1. PREPARATION
        // =================================================================================================
        const Real dt = parent->getDt();
        const Vec3i gridSize = parent->getGridSize();
        const int bnd = 1; // Boundary width to ignore in loops

        // Temporary grids for the new advected quantity, and for tracking weights
        GridType grid_new(parent, 0.0);
        Grid<Real> beta(parent, 0.0);
        Grid<Real> gamma_advected(parent);

        // =================================================================================================
        // 2. ADVECTION of CUMULATIVE GAMMA (y)
        // Advect the cumulative gamma field from the previous step using standard semi-Lagrangian advection.
        // This provides the initial state for the gamma field at t_{n+1}.
        // =================================================================================================
        FOR_IJK(grid)
        {
            const Vec3i p(i, j, k);
            if (flags_n_plus_one.isFluid(p))
            {
                Vec3 pos = toVec3(p) + offset;
                Vec3 pos_start = vel.getInterpolated(pos) * -dt + pos;
                gamma_advected(p) = gammaCumulative.getInterpolated(pos_start);
            }
            else
            {
                gamma_advected(p) = 1.0; // Default for non-fluid cells
            }
        }

        // =================================================================================================
        // 3. BACKWARD PASS (Main Advection)
        // Advect 'grid' from t_n to t_{n+1}. Source points are filtered to be inside the fluid domain (flags_n).
        // Interpolation weights are renormalized to ensure we only gather from valid fluid cells.
        // We compute 'beta', the sum of outgoing weights from each source cell.
        // =================================================================================================
        FOR_IJK(grid)
        {
            const Vec3i p_dest(i, j, k);
            if (!flags_n_plus_one.isFluid(p_dest))
                continue;

            Vec3 pos_dest = toVec3(p_dest) + offset;
            Vec3 pos_start = vel.getInterpolated(pos_dest) * -dt + pos_dest;

            Vec3i s_base = toVec3i(pos_start - offset);
            Vec3 f = (pos_start - offset) - toVec3(s_base);

            // Collect source cells and their weights that are within the fluid domain
            Real total_weight = 0;
            Vec3i sources_s[4];
            Real sources_w[4];
            int num_valid_sources = 0;

            for (int jj = 0; jj <= 1; ++jj)
            {
                for (int ii = 0; ii <= 1; ++ii)
                {
                    Vec3i s_curr = s_base + Vec3i(ii, jj, 0);
                    if (grid.isInBounds(s_curr, bnd) && flags_n.isFluid(s_curr))
                    {
                        sources_s[num_valid_sources] = s_curr;
                        sources_w[num_valid_sources] = getBilinearWeight(f, ii, jj);
                        total_weight += sources_w[num_valid_sources];
                        num_valid_sources++;
                    }
                }
            }

            if (total_weight < 1e-7)
            {
                // Fallback for isolated cells: use standard (less accurate) advection
                grid_new(p_dest) = grid.getInterpolated(pos_start);
                continue;
            }

            // Apply normalized weights to compute the new grid value and update beta
            for (int k_idx = 0; k_idx < num_valid_sources; ++k_idx)
            {
                const Real w_norm = sources_w[k_idx] / total_weight;
                const Vec3i s_curr = sources_s[k_idx];
                grid_new(p_dest) += grid(s_curr) * w_norm;
                beta(s_curr) += w_norm;
            }
        }

        // =================================================================================================
        // 4. FORWARD PASS (Conservation Correction)
        // For source cells where beta < 1, some of their quantity was not advected.
        // We advect this residual amount forward in time and add it to the result.
        // This ensures mass/momentum conservation (all of 'grid' is moved).
        // =================================================================================================
        GridType R(parent, 0.0);
        Grid<Real> R_gamma(parent, 0.0);
        FOR_IJK(grid)
        {
            const Vec3i p_src(i, j, k);
            if (flags_n.isFluid(p_src) && beta(p_src) < 0.9999)
            {
                const Real missing_frac = 1.0 - beta(p_src);
                if (missing_frac > 0)
                {
                    R(p_src) = grid(p_src) * missing_frac;
                    R_gamma(p_src) = gammaCumulative(p_src) * missing_frac;
                }
            }
        }

        FOR_IJK(grid)
        {
            const Vec3i p_src(i, j, k);
            if (R_gamma(p_src) == 0.0)
                continue; // no residual to advect

            Vec3 pos_src = toVec3(p_src) + offset;
            Vec3 pos_end = vel.getInterpolated(pos_src) * dt + pos_src;

            // Distribute residual R and R_gamma to valid fluid destination cells
            Vec3i d_base = toVec3i(pos_end - offset);
            Vec3 f = (pos_end - offset) - toVec3(d_base);

            Real total_weight = 0;
            Vec3i dests_d[4];
            Real dests_w[4];
            int num_valid_dests = 0;

            for (int jj = 0; jj <= 1; ++jj)
            {
                for (int ii = 0; ii <= 1; ++ii)
                {
                    Vec3i d_curr = d_base + Vec3i(ii, jj, 0);
                    if (grid.isInBounds(d_curr, bnd) && flags_n_plus_one.isFluid(d_curr))
                    {
                        dests_d[num_valid_dests] = d_curr;
                        dests_w[num_valid_dests] = getBilinearWeight(f, ii, jj);
                        total_weight += dests_w[num_valid_dests];
                        num_valid_dests++;
                    }
                }
            }

            if (total_weight < 1e-7)
                continue; // No valid fluid destination

            for (int k_idx = 0; k_idx < num_valid_dests; ++k_idx)
            {
                const Real w_norm = dests_w[k_idx] / total_weight;
                const Vec3i d_curr = dests_d[k_idx];
                grid_new(d_curr) += R(p_src) * w_norm;
                gamma_advected(d_curr) += R_gamma(p_src) * w_norm;
            }
        }

        // =================================================================================================
        // 5. DIFFUSION PASS (Incompressibility Correction)
        // Iteratively smooth the 'gamma_advected' field towards a constant value (1.0).
        // For each diffusion step on gamma, a corresponding amount of the quantity 'grid_new'
        // is moved between cells to maintain consistency. This enforces incompressibility.
        // =================================================================================================
        const int num_diff_iter = 5; // Typically 1-7 iterations are sufficient
        for (int iter = 0; iter < num_diff_iter; ++iter)
        {
            // X-sweep
            for (int j = bnd; j < gridSize.y - bnd; ++j)
            {
                for (int i = bnd; i < gridSize.x - 1 - bnd; ++i)
                {
                    if (flags_n_plus_one.isFluid(i, j, 0) && flags_n_plus_one.isFluid(i + 1, j, 0))
                    {
                        Real &g1 = gamma_advected(i, j, 0);
                        Real &g2 = gamma_advected(i + 1, j, 0);

                        if (g2 < 1e-7)
                            continue;

                        const Real avg = (g1 + g2) * 0.5;
                        const auto delta_q = grid_new(i + 1, j, 0) * ((g2 - g1) * 0.5 / g2);

                        grid_new(i, j, 0) += delta_q;
                        grid_new(i + 1, j, 0) -= delta_q;

                        g1 = avg;
                        g2 = avg;
                    }
                }
            }
            // Y-sweep
            for (int i = bnd; i < gridSize.x - bnd; ++i)
            {
                for (int j = bnd; j < gridSize.y - 1 - bnd; ++j)
                {
                    if (flags_n_plus_one.isFluid(i, j, 0) && flags_n_plus_one.isFluid(i, j + 1, 0))
                    {
                        Real &g1 = gamma_advected(i, j, 0);
                        Real &g2 = gamma_advected(i, j + 1, 0);

                        if (g2 < 1e-7)
                            continue;

                        const Real avg = (g1 + g2) * 0.5;
                        const auto delta_q = grid_new(i, j + 1, 0) * ((g2 - g1) * 0.5 / g2);

                        grid_new(i, j, 0) += delta_q;
                        grid_new(i, j + 1, 0) -= delta_q;

                        g1 = avg;
                        g2 = avg;
                    }
                }
            }
        }

        // =================================================================================================
        // 6. FINALIZE
        // Swap the temporary grids with the output grids.
        // =================================================================================================
        grid.swap(grid_new);
        gammaCumulative.swap(gamma_advected);
    }

    // Explicit instantiations if needed for linking
    // template void fnMassMomentumConservingAdvectWater<MACGrid>(...);
    // template void fnMassMomentumConservingAdvectWater<Grid<Real>>(...);

    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code
    // gemini code

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
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velX, gammaX, offsetX);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velY, gammaY, offsetY);
            fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velZ, gammaZ, offsetZ);
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
            fnMassMomentumConservingAdvect<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
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
        advectGammaCum<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags);
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

        knSetNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
        grid.swap(newGrid);
    }
}