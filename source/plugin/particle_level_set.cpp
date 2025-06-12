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
    void knCombineOmegas(Grid<Real> &omegaPlus, Grid<Real> &omegaMinus, Grid<Real> &phi)
    {
        phi(i, j, k) = std::abs(omegaPlus(i, j, k)) <= std::abs(omegaMinus(i, j, k)) ? omegaPlus(i, j, k) : omegaMinus(i, j, k);
    }

    PYTHON()
    void correctErrorsWithParticles(Grid<Real> &phi, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii)
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
            knCombineOmegas(omegaPlus, phi, phi);
        }
        else if (negativeEscaped.size() > 0 && positiveEscaped.size() == 0)
        {
            knCombineOmegas(phi, omegaMinus, phi);
        }
        else if (positiveEscaped.size() > 0 && negativeEscaped.size() > 0)
        {
            knCombineOmegas(omegaPlus, omegaMinus, phi);
        }
    }

    Real calcGradMagnitude(Grid<Real> &phi, int i, int j, int k, Real sng)
    {
        Real dx = 1.;

        // Get directional derivatives (forward and backward differences)
        Real dPhiDx_minus = (phi(i, j, k) - phi(i - 1, j, k)) / dx;
        Real dPhiDx_plus = (phi(i + 1, j, k) - phi(i, j, k)) / dx;
        Real dPhiDy_minus = (phi(i, j, k) - phi(i, j - 1, k)) / dx;
        Real dPhiDy_plus = (phi(i, j + 1, k) - phi(i, j, k)) / dx;

        Real grad_x_sq = 0.0;
        Real grad_y_sq = 0.0;

        // Choose the upwind derivative based on the sign 's'
        if (sng > 0.0)
        {
            // Information flows outwards from the interface. We look inwards for the gradient.
            grad_x_sq = std::max(std::pow(std::max(dPhiDx_minus, 0.0f), 2), std::pow(std::min(dPhiDx_plus, 0.0f), 2));
            grad_y_sq = std::max(std::pow(std::max(dPhiDy_minus, 0.0f), 2), std::pow(std::min(dPhiDy_plus, 0.0f), 2));
        }
        else if (sng < 0.0)
        {
            // Information flows inwards towards the center. We look outwards for the gradient.
            grad_x_sq = std::max(std::pow(std::min(dPhiDx_minus, 0.0f), 2), std::pow(std::max(dPhiDx_plus, 0.0f), 2));
            grad_y_sq = std::max(std::pow(std::min(dPhiDy_minus, 0.0f), 2), std::pow(std::max(dPhiDy_plus, 0.0f), 2));
        }
        // If s is exactly 0, the gradient is 0, as handled by initialization.

        return std::sqrt(grad_x_sq + grad_y_sq);
    }

    KERNEL(bnd = 1)
    void knReinitializeLevelset(Grid<Real> &phi, Grid<Real> &phiOld, Real dt_pseudo)
    {
        Real p = phiOld(i, j, k);
        Real sng = p / std::sqrt(p * p + 1 * 1);

        phi(i, j, k) = phiOld(i, j, k) - dt_pseudo * sng * (calcGradMagnitude(phiOld, i, j, k, sng) - 1);
    }

    PYTHON()
    void reinitializeLevelset(Grid<Real> &phi)
    {
        Grid<Real> phiOld(phi.getParent());
        Real dt_pseudo = 0.5;
        for (int __ = 0; __ < 10; __++)
        {
            phiOld.copyFrom(phi);
            knReinitializeLevelset(phi, phiOld, dt_pseudo);
        }
    }

    PYTHON()
    void reinitializeRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        knSetParticleRadii(particles, radii, phi);
    }

    PYTHON()
    void reseedParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii,
                         Real cutoff = 2., int discretization = 2, Real randomness = 0.5)
    {
        
    }
}