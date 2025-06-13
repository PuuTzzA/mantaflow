#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
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

    KERNEL(points)
    void knSetParticleRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        skipDelted(particles);
        radii[idx] = std::abs(phi.getInterpolatedHi(particles[idx].pos, 2));
    }

    PYTHON()
    void sampleLevelsetBorderWithParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii,
                                           Real cutoff = 2., int discretization = 2, Real randomness = 0.5)
    {
        Real step = 1. / discretization;
        FOR_IJK_BND(flags, 0)
        {
            if (flags.isObstacle(i, j, k))
            {
                continue;
            }
            if (std::abs(phi(i, j, k)) < cutoff)
            {
                for (int di = 0; di < discretization; di++)
                {
                    for (int dj = 0; dj < discretization; dj++)
                    {
                        Real randi = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;
                        Real randj = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;

                        Vec3 pos = Vec3(i + (0.5 + di) * step + randi, j + (0.5 + dj) * step + randj, 0.5);
                        particles.addBuffered(pos, phi(i, j, k) < 0 ? ParticleBase::PINSIDE : ParticleBase::POUTSIDE);
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
        static Real escapeCondition = 1.;
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

            if (std::abs(phiAtPos) > escapeCondition * radii[idx])
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
    void reseedParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii,
                         Real cutoff = 2., Real randomness = 0.5)
    {
        const static int minParticles = 2;
        const static int maxParticles = 8;

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

            if (std::abs(phi(i, j, k)) > cutoff || flags.isObstacle(i, j, k))
            {
                particles.kill(idx);
                ++deletedParticles;
            }

            if (++counterGrid(i, j, k) > maxParticles)
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
            if (std::abs(phi(i, j, k)) > cutoff)
            {
                continue;
            }
            for (int _ = 0; _ < minParticles + 1 - counterGrid(i, j, k); ++_)
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
}