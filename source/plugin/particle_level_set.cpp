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
    /*     KERNEL()
        void knSampleLevelSet(const LevelsetGrid &phi, const FlagGrid &flags, BasicParticleSystem &particles, Real cutoff = 3., Real randomness = 0.05)
        {
            if (!flags.isObstacle(i, j, k))
            {
            }
        } */

    PYTHON()
    void sampleLevelsetBorderWithParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, Real cutoff = 3., int discretization = 2, Real randomness = 0.5)
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

                        Vec3 pos = Vec3(i + (0.5 + di) * step + randi, j + (0.5 + dj) * step + randj, 0);
                        particles.addBuffered(pos, phi(i, j, k) < 0 ? ParticleBase::PINSIDE : ParticleBase::POUTSIDE);
                    }
                }
            }
        }
        particles.insertBufferedParticles();
    }

    KERNEL(points)
    void testXX(BasicParticleSystem &particles, Grid<Real> &g)
    {
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

}