#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
#include "mass_and_momentum_conserving_advection.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <cstdlib>

namespace Manta
{
#define skipDelted(particles)                        \
    if (particles[idx].flag & ParticleBase::PDELETE) \
    {                                                \
        return;                                      \
    }

#define DISCRETIZATION 2
#define POSITIVE_CUTOFF 0.f
#define NEGATIVE_CUTOFF -3.f
#define MAX_PARTICLES_PER_CELL 10
#define MIN_PARTICLES_PER_CELL 3
#define ESCAPE_CONDITION .1 // a particle counts as escaped, if it is further than ESCAPE_CONDITION * radius on the wrong side

    KERNEL(points)
    void knSetParticleRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        skipDelted(particles);
        radii[idx] = std::abs(Manta::clamp(phi.getInterpolatedHi(particles[idx].pos, 2), NEGATIVE_CUTOFF, POSITIVE_CUTOFF));
    }

    PYTHON()
    void sampleLevelsetBorderWithParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, Real randomness = 0.5)
    {
        Real step = 1. / DISCRETIZATION;
        FOR_IJK_BND(flags, 0)
        {
            if (flags.isObstacle(i, j, k))
            {
                continue;
            }
            if (0 < phi(i, j, k) && phi(i, j, k) < POSITIVE_CUTOFF)
            {
                for (int di = 0; di < DISCRETIZATION; di++)
                {
                    for (int dj = 0; dj < DISCRETIZATION; dj++)
                    {
                        Real randi = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;
                        Real randj = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;

                        Vec3 pos = Vec3(i + (0.5 + di) * step + randi, j + (0.5 + dj) * step + randj, 0.5);
                        particles.addBuffered(pos, phi(i, j, k) < ParticleBase::POUTSIDE);
                    }
                }
            }
            else if (NEGATIVE_CUTOFF < phi(i, j, k) && phi(i, j, k) <= 0)
            {
                for (int di = 0; di < DISCRETIZATION; di++)
                {
                    for (int dj = 0; dj < DISCRETIZATION; dj++)
                    {
                        Real randi = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;
                        Real randj = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;

                        Vec3 pos = Vec3(i + (0.5 + di) * step + randi, j + (0.5 + dj) * step + randj, 0.5);
                        particles.addBuffered(pos, ParticleBase::PINSIDE);
                    }
                }
            }
        }
        particles.insertBufferedParticles();
        knSetParticleRadii(particles, radii, phi);
    }

    KERNEL(points)
    void testXX(BasicParticleSystem &particles, Grid<Real> &g)
    {
        skipDelted(particles);
        Vec3 pos = particles[idx].pos;
        g(std::floor(pos.x), std::floor(pos.y), 0) = particles[idx].flag == 0 ? 0. : 1.;
        // g(std::floor(pos.x), std::floor(pos.y), 0) = 0.;
    }

    PYTHON()
    void testSeedParticles(const LevelsetGrid &phi, const FlagGrid &flags, Grid<Real> &g, BasicParticleSystem &particles, Real cutoff = 3., Real randomness = 0.05)
    {
        // testXX(particles, g);

        for (IndexInt idx = 0, max = particles.size(); idx < max; idx++)
        {
            Vec3 pos = particles[idx].pos;
            g(std::floor(pos.x), std::floor(pos.y), 0) = particles[idx].flag & ParticleBase::POUTSIDE ? 0.5 : 1.;
        }
    }

    inline Real particleSphere(int flag, Real radius, Vec3 particlePosition, Vec3 x)
    {
        Real sp = flag & ParticleBase::PINSIDE ? -1. : 1.;
        return sp * (radius - norm(x - particlePosition));
    }

    KERNEL()
    void knCorrectOmegas(Grid<Real> &originalPhi, Grid<Real> &omega, std::vector<IndexInt> escapedParticles, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, bool useMin)
    {
        omega(i, j, k) = useMin ? std::numeric_limits<Real>::infinity() : -std::numeric_limits<Real>::infinity();

        for (auto idx : escapedParticles)
        {
            if (useMin)
            {
                omega(i, j, k) = std::min(originalPhi(i, j, k), particleSphere(particles[idx].flag, radii[idx], particles[idx].pos, Vec3(i + 0.5, j + 0.5, k + 0.5)));
            }
            else
            {
                omega(i, j, k) = std::max(originalPhi(i, j, k), particleSphere(particles[idx].flag, radii[idx], particles[idx].pos, Vec3(i + 0.5, j + 0.5, k + 0.5)));
            }
        }
    }

    KERNEL()
    void knCombineOmegas(Grid<Real> &omegaPlus, Grid<Real> &omegaMinus, Grid<Real> &phi, const FlagGrid &flags)
    {
        if (flags.isObstacle(i, j, k))
        {
            return;
        }
        phi(i, j, k) = std::abs(omegaPlus(i, j, k)) <= std::abs(omegaMinus(i, j, k)) ? omegaPlus(i, j, k) : omegaMinus(i, j, k);
    }

    PYTHON()
    void correctErrorsWithParticles(Grid<Real> &phi, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const FlagGrid &flags)
    {
        // Step 1: find escaped particles
        std::vector<IndexInt> positiveEscaped{}; // escaped particles with flag POUTSIDE, but are escapeConditio * radius inside
        std::vector<IndexInt> negativeEscaped{}; // escaped particles with flag PINSIDE, but are escapeConditio * radius outside

        for (IndexInt idx = 0, max = particles.size(); idx < max; idx++)
        {
            Vec3 pos = particles[idx].pos;
            Real phiAtPos = phi.getInterpolatedHi(pos, 2);

            if (particles[idx].flag & ParticleBase::PDELETE)
            {
                continue; // delted
            }

            if (particles[idx].flag & (phiAtPos < 0 ? ParticleBase::PINSIDE : ParticleBase::POUTSIDE))
            {
                continue; // on the right side
            }

            if (std::abs(phiAtPos) > ESCAPE_CONDITION * radii[idx])
            {
                continue; // particle still not far enough on wrong side
            }

            particles[idx].flag |= ParticleBase::PESCAPED;
            if (particles[idx].flag & ParticleBase::PINSIDE)
            {
                // Particles that should be inside but are outside
                negativeEscaped.push_back(idx);
            }
            else
            {
                // Particles that should be outside but are inside
                positiveEscaped.push_back(idx);
            }
        }

        std::cout << "POSITIVE ESCAPED: " << positiveEscaped.size() << std::endl;
        std::cout << "NEGATIVE ESCAPED: " << negativeEscaped.size() << std::endl;

        /*         for (auto p : negativeEscaped)
                {
                    std::cout << "pos: " << (particles[p].pos) << "; radius: " << radii[p] << "; stuff: " << static_cast<bool>(particles[p].flag & ParticleBase::PINSIDE) << std::endl;
                } */

        // Step 2: correct omega+ and omega-
        Grid<Real> omegaPlus(phi.getParent());
        Grid<Real> omegaMinus(phi.getParent());

        // positiveEscaped.clear();

        if (positiveEscaped.size() > 0)
        {
            knCorrectOmegas(phi, omegaPlus, positiveEscaped, particles, radii, false);
        }
        if (negativeEscaped.size() > 0)
        {
            knCorrectOmegas(phi, omegaMinus, negativeEscaped, particles, radii, true);
        }

        if (positiveEscaped.size() > 0 && negativeEscaped.size() == 0)
        {
            knCombineOmegas(omegaPlus, phi, phi, flags);
        }
        else if (negativeEscaped.size() > 0 && positiveEscaped.size() == 0)
        {
            knCombineOmegas(phi, omegaMinus, phi, flags);
        }
        else if (positiveEscaped.size() > 0 && negativeEscaped.size() > 0)
        {
            knCombineOmegas(omegaPlus, omegaMinus, phi, flags);
        }
    }

    Real calcGradMagnitude(Grid<Real> &phi, int i, int j, int k, Real sng, const FlagGrid &flags)
    {
        // TODO: UNDERSTAND AND LOOK OVER THIS ALGORITHM (this one was made by google aistudio)

        // --- Calculate directional derivatives, respecting obstacles ---
        Real dx = 1.;
        // Backward difference in x (dPhi/dx-)
        Real dPhiDx_minus = 0.0;
        if (!flags.isObstacle(i - 1, j, k))
        {
            dPhiDx_minus = (phi(i, j, k) - phi(i - 1, j, k)) / dx;
        }

        // Forward difference in x (dPhi/dx+)
        Real dPhiDx_plus = 0.0;
        if (!flags.isObstacle(i + 1, j, k))
        {
            dPhiDx_plus = (phi(i + 1, j, k) - phi(i, j, k)) / dx;
        }

        // Backward difference in y (dPhi/dy-)
        Real dPhiDy_minus = 0.0;
        if (!flags.isObstacle(i, j - 1, k))
        {
            dPhiDy_minus = (phi(i, j, k) - phi(i, j - 1, k)) / dx;
        }

        // Forward difference in y (dPhi/dy+)
        Real dPhiDy_plus = 0.0;
        if (!flags.isObstacle(i, j + 1, k))
        {
            dPhiDy_plus = (phi(i, j + 1, k) - phi(i, j, k)) / dx;
        }

        // --- The Godunov scheme part is identical to before ---
        // It now operates on derivatives that are zeroed out at obstacle boundaries.

        Real grad_x_sq = 0.0;
        Real grad_y_sq = 0.0;

        if (sng > 0.0)
        {
            grad_x_sq = std::max(std::pow(std::max(dPhiDx_minus, 0.0f), 2), std::pow(std::min(dPhiDx_plus, 0.0f), 2));
            grad_y_sq = std::max(std::pow(std::max(dPhiDy_minus, 0.0f), 2), std::pow(std::min(dPhiDy_plus, 0.0f), 2));
        }
        else if (sng < 0.0)
        {
            grad_x_sq = std::max(std::pow(std::min(dPhiDx_minus, 0.0f), 2), std::pow(std::max(dPhiDx_plus, 0.0f), 2));
            grad_y_sq = std::max(std::pow(std::min(dPhiDy_minus, 0.0f), 2), std::pow(std::max(dPhiDy_plus, 0.0f), 2));
        }
        // If s is exactly 0, the gradient remains 0.

        return std::sqrt(grad_x_sq + grad_y_sq);
    }

    KERNEL(bnd = 1)
    void knReinitializeLevelset(Grid<Real> &phi, Grid<Real> &phiOld, Real dt_pseudo, const FlagGrid &flags)
    {
        if (flags.isObstacle(i, j, k))
        {
            return;
        }

        Real p = phiOld(i, j, k);
        Real sng = p / std::sqrt(p * p + 1 * 1); // smooth sng function

        phi(i, j, k) = phiOld(i, j, k) - dt_pseudo * sng * (calcGradMagnitude(phiOld, i, j, k, sng, flags) - 1);
    }

    PYTHON()
    void reinitializeLevelset(Grid<Real> &phi, const FlagGrid &flags)
    {
        Grid<Real> phiOld(phi.getParent());
        Real dt_pseudo = 0.5;
        for (int __ = 0; __ < 10; __++)
        {
            phiOld.copyFrom(phi);
            knReinitializeLevelset(phi, phiOld, dt_pseudo, flags);
        }
    }

    PYTHON()
    void reinitializeRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        knSetParticleRadii(particles, radii, phi);
    }

    PYTHON()
    //! you always need to run reinitializeRadii after this
    void reseedParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, Real randomness = 0.5)
    {
        Grid<Real> counterGrid(phi.getParent());

        int deletedParticles = 0;
        int addedParticles = 0;

        // Delete Out of bounds particles or overcrowded particles
        for (IndexInt idx = 0, max = particles.size(); idx < max; idx++)
        {
            if (particles[idx].flag & ParticleBase::PESCAPED)
            {
                continue; // don't delete escaped particles
            }
            int i = std::floor(particles[idx].pos[0]);
            int j = std::floor(particles[idx].pos[1]);
            int k = std::floor(particles[idx].pos[2]);

            if (flags.isObstacle(i, j, k))
            {
                continue;
            }

            if (!(NEGATIVE_CUTOFF < phi(i, j, k) && phi(i, j, k) < POSITIVE_CUTOFF))
            {
                particles.kill(idx);
                ++deletedParticles;
            }

            if (++counterGrid(i, j, k) > MAX_PARTICLES_PER_CELL)
            {
                particles.kill(idx);
                ++deletedParticles;
            }
        }

        // Insert new particles in cells where there are not enough particles
        FOR_IJK(phi)
        {
            if (flags.isObstacle(i, j, k))
            {
                continue;
            }
            if (!(NEGATIVE_CUTOFF < phi(i, j, k) && phi(i, j, k) < POSITIVE_CUTOFF))
            {
                continue;
            }
            for (int _ = 0; _ < MIN_PARTICLES_PER_CELL - counterGrid(i, j, k); ++_)
            {
                ++addedParticles;
                Real randi = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * randomness;
                Real randj = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * randomness;

                Vec3 pos = Vec3(i + 0.5 + randi, j + 0.5 + randj, 0.5);
                particles.addBuffered(pos, phi(i, j, k) < 0 ? ParticleBase::PINSIDE : ParticleBase::POUTSIDE);
            }
        }
        particles.insertBufferedParticles();

        std::cout << "Deleted " << deletedParticles << " particles." << std::endl;
        std::cout << "Created " << addedParticles << " particles." << std::endl;
    }

    KERNEL()
    void knSetFlagsFromParticleLevelset(const Grid<Real> &phi, FlagGrid &flags, Real level)
    {
        if (flags.isFluid(i, j, k) || flags.isEmpty(i, j, k))
        {
            flags(i, j, k) = (flags(i, j, k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid; // clears the fluid flags

            if (phi(i, j, k) < level)
            {
                flags(i, j, k) = (flags(i, j, k) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
            }
        }
    }

    PYTHON()
    void setFlagsFromParticleLevelset(const Grid<Real> &phi, FlagGrid &flags, Real level)
    {
        knSetFlagsFromParticleLevelset(phi, flags, level);
    }

    KERNEL()
    void knFluidLSFillHelper(Grid<Real> &grid, const FlagGrid &flags, const Grid<Real> &phi, Real level)
    {
        if (!flags.isObstacle(i, j, k) && phi(i, j, k) < level)
        {
            grid(i, j, k) = 1.;
        }
        else
        {
            grid(i, j, k) = 0.;
        }
    }

    PYTHON()
    void fillLevelsetWithOnes(GridBase *grid, const FlagGrid *flags, const Grid<Real> *phi, Real level)
    {
        if (grid->getType() & GridBase::TypeReal)
        {
            knFluidLSFillHelper(*((Grid<Real> *)grid), *flags, *phi, level);
            std::cout << "filled fluid based on phi in Grid<Real> with ones" << std::endl;
        }
        else
        {
            throw std::runtime_error("fill fluid with ones not implemented for this grid type!");
        }
    }

    KERNEL()
    void knFillSolvedFlags(const FlagGrid &flags, Grid<int> &solvedGrid, Grid<Real> &velComponent, MACGridComponent component)
    {
        const Real D_INF = std::numeric_limits<Real>::max();
        bool isSolid;
        bool isFluid;

        // direction: 0=x, 1=y, 3=z
        switch (component)
        {
        case MAC_X:
            isSolid = flags.isObstacle(i - 1, j, k) || flags.isObstacle(i, j, k);
            isFluid = flags.isFluid(i - 1, j, k) || flags.isFluid(i, j, k);
            break;
        case MAC_Y:
            isSolid = flags.isObstacle(i, j - 1, k) || flags.isObstacle(i, j, k);
            isFluid = flags.isFluid(i, j - 1, k) || flags.isFluid(i, j, k);
            break;
        case MAC_Z:
            throw std::runtime_error("not implemented yet! (getPhiAtMACCell with MAC_Z)");
        }

        if (isSolid || isFluid)
        {
            solvedGrid(i, j, k) = true;
            if (isSolid)
            {
                velComponent(i, j, k) = 0.;
            }
        }
        else
        {
            solvedGrid(i, j, k) = false;
            velComponent(i, j, k) = D_INF;
        }
    }

    inline Real getPhiAtMACCell(IndexInt i, IndexInt j, IndexInt k, const Grid<Real> &phi, MACGridComponent component)
    {
        switch (component)
        {
        case MAC_X:
            return 0.5 * (phi(i, j, k) + phi(i - 1, j, k));
        case MAC_Y:
            return 0.5 * (phi(i, j, k) + phi(i, j - 1, k));
        case MAC_Z:
            throw std::runtime_error("not implemented yet! (getPhiAtMACCell with MAC_Z)");
        default:
            throw std::runtime_error("should never happen");
        }
    }

    KERNEL()
    void knRemoveInfinities(Grid<Real> &grid, int _)
    {
        grid(i, j, k) = grid(i, j, k) == std::numeric_limits<Real>::max() ? 0. : grid(i, j, k);
    }

    void extrapolateComponentFSM(Grid<Real> &velComponent, const Grid<Real> &phi, Grid<int> &solvedMask, int sweepIterations, Vec3i gs, MACGridComponent component)
    {
        const Real D_INF = std::numeric_limits<Real>::max();

        for (int iter = 0; iter < sweepIterations; iter++)
        {
            // Sweep 1 (i++, j++)
            for (IndexInt j = 1; j < gs[1] - 1; ++j)
            {
                for (IndexInt i = 1; i < gs[0] - 1; ++i)
                {
                    IndexInt k = 0;
                    if (solvedMask(i, j, k))
                    {
                        continue;
                    }

                    Real sumVal = 0.;
                    Real sumW = 0.;

                    Real phiHere = getPhiAtMACCell(i, j, k, phi, component);

                    if (velComponent(i - 1, j, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i - 1, j, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i - 1, j, k);
                            sumW += w;
                        }
                    }

                    if (velComponent(i, j - 1, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i, j - 1, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i, j - 1, k);
                            sumW += w;
                        }
                    }

                    if (sumW > 1e-6)
                    {
                        velComponent(i, j, k) = std::min(sumVal / sumW, velComponent(i, j, k));
                    }
                }
            }

            // Sweep 2 (i--, j++)
            for (IndexInt j = 1; j < gs[1] - 1; ++j)
            {
                for (IndexInt i = gs[0] - 2; i >= 1; --i)
                {
                    IndexInt k = 0;
                    if (solvedMask(i, j, k))
                    {
                        continue;
                    }

                    Real sumVal = 0.;
                    Real sumW = 0.;

                    Real phiHere = getPhiAtMACCell(i, j, k, phi, component);

                    if (velComponent(i + 1, j, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i + 1, j, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i + 1, j, k);
                            sumW += w;
                        }
                    }

                    if (velComponent(i, j - 1, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i, j - 1, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i, j - 1, k);
                            sumW += w;
                        }
                    }

                    if (sumW > 1e-6)
                    {
                        velComponent(i, j, k) = std::min(sumVal / sumW, velComponent(i, j, k));
                    }
                }
            }

            // Sweep 3 (i++, j--)
            for (IndexInt j = gs[1] - 2; j >= 1; --j)
            {
                for (IndexInt i = 1; i < gs[0] - 1; ++i)
                {
                    IndexInt k = 0;
                    if (solvedMask(i, j, k))
                    {
                        continue;
                    }

                    Real sumVal = 0.;
                    Real sumW = 0.;

                    Real phiHere = getPhiAtMACCell(i, j, k, phi, component);

                    if (velComponent(i - 1, j, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i - 1, j, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i - 1, j, k);
                            sumW += w;
                        }
                    }

                    if (velComponent(i, j + 1, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i, j + 1, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i, j + 1, k);
                            sumW += w;
                        }
                    }

                    if (sumW > 1e-6)
                    {
                        velComponent(i, j, k) = std::min(sumVal / sumW, velComponent(i, j, k));
                    }
                }
            }

            // Sweep 4 (i--, j--)
            for (IndexInt j = gs[1] - 2; j >= 1; --j)
            {
                for (IndexInt i = gs[0] - 2; i >= 1; --i)
                {
                    IndexInt k = 0;
                    if (solvedMask(i, j, k))
                    {
                        continue;
                    }

                    Real sumVal = 0.;
                    Real sumW = 0.;

                    Real phiHere = getPhiAtMACCell(i, j, k, phi, component);

                    if (velComponent(i + 1, j, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i + 1, j, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i + 1, j, k);
                            sumW += w;
                        }
                    }

                    if (velComponent(i, j + 1, k) < D_INF)
                    {
                        Real phiNeighbor = getPhiAtMACCell(i, j + 1, k, phi, component);
                        Real w = std::max(0.f, phiHere - phiNeighbor);
                        if (w > 1e-6)
                        {
                            sumVal += w * velComponent(i, j + 1, k);
                            sumW += w;
                        }
                    }

                    if (sumW > 1e-6)
                    {
                        velComponent(i, j, k) = std::min(sumVal / sumW, velComponent(i, j, k));
                    }
                }
            }
        }

        knRemoveInfinities(velComponent, 0);
    }

    PYTHON()
    void extrapolateVelFSM(const Grid<Real> &phi, const FlagGrid &flags, MACGrid &vel, int steps)
    {
        FluidSolver *parent = phi.getParent();
        Grid<Real> velX(parent);
        Grid<Real> velY(parent);
        Grid<Real> velZ(parent);

        knMAC2Grids(vel, velX, velY, velZ);

        Vec3i gs = parent->getGridSize();
        Real dx = 1.;

        Grid<int> solvedFlags(parent);

        // X-component
        knFillSolvedFlags(flags, solvedFlags, velX, MAC_X);
        extrapolateComponentFSM(velX, phi, solvedFlags, steps, gs, MAC_X);

        // Y-component
        knFillSolvedFlags(flags, solvedFlags, velY, MAC_Y);
        extrapolateComponentFSM(velY, phi, solvedFlags, steps, gs, MAC_Y);

        knGrids2MAC(vel, velX, velY, velZ, flags);
    }
}