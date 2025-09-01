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
#include <random>
#include <queue>
#include <assert.h>

namespace Manta
{
#define skipDelted(particles)                        \
    if (particles[idx].flag & ParticleBase::PDELETE) \
    {                                                \
        return;                                      \
    }

#define POSITIVE_SEED_CUTOFF 1.0f  // 3 * dx
#define NEGATIVE_SEED_CUTOFF -3.0f // 3 * dx
#define MIN_RADIUS 0.1f            // .1 * min(dx, dy, dz)
#define MAX_RADIUS 0.5f            // .5 * min(dx, dy, dz)

#define DISCRETIZATION 4 // particles per spatial dimension

#define ESCAPE_CONDITION 1 // a particle counts as escaped, if it is further than ESCAPE_CONDITION * radius on the wrong side

#define RANDOM_SEED 1

    // LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET
    // LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET
    // LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET
    // LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET
    // LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET LEVEL SET

    Real calcGradMagnitude(Grid<Real> &phi, int i, int j, int k, Real sng, const FlagGrid &flags)
    {
        Real dx = 1.0;

        // Backward differences
        Real dPhiDx_minus = 0.0, dPhiDy_minus = 0.0, dPhiDz_minus = 0.0;
        if (!flags.isObstacle(i - 1, j, k))
            dPhiDx_minus = (phi(i, j, k) - phi(i - 1, j, k)) / dx;
        if (!flags.isObstacle(i, j - 1, k))
            dPhiDy_minus = (phi(i, j, k) - phi(i, j - 1, k)) / dx;
        if (phi.is3D() && !flags.isObstacle(i, j, k - 1))
            dPhiDz_minus = (phi(i, j, k) - phi(i, j, k - 1)) / dx;

        // Forward differences
        Real dPhiDx_plus = 0.0, dPhiDy_plus = 0.0, dPhiDz_plus = 0.0;
        if (!flags.isObstacle(i + 1, j, k))
            dPhiDx_plus = (phi(i + 1, j, k) - phi(i, j, k)) / dx;
        if (!flags.isObstacle(i, j + 1, k))
            dPhiDy_plus = (phi(i, j + 1, k) - phi(i, j, k)) / dx;
        if (phi.is3D() && !flags.isObstacle(i, j, k + 1))
            dPhiDz_plus = (phi(i, j, k + 1) - phi(i, j, k)) / dx;

        // Godunov scheme
        Real grad_x_sq = 0.0, grad_y_sq = 0.0, grad_z_sq = 0.0;

        if (sng > 0.0)
        {
            grad_x_sq = std::max(std::pow(std::max(dPhiDx_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::min(dPhiDx_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
            grad_y_sq = std::max(std::pow(std::max(dPhiDy_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::min(dPhiDy_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
            grad_z_sq = std::max(std::pow(std::max(dPhiDz_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::min(dPhiDz_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
        }
        else if (sng < 0.0)
        {
            grad_x_sq = std::max(std::pow(std::min(dPhiDx_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::max(dPhiDx_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
            grad_y_sq = std::max(std::pow(std::min(dPhiDy_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::max(dPhiDy_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
            grad_z_sq = std::max(std::pow(std::min(dPhiDz_minus, static_cast<Real>(0.0)), static_cast<Real>(2)),
                                 std::pow(std::max(dPhiDz_plus, static_cast<Real>(0.0)), static_cast<Real>(2)));
        }
        // else if sng == 0.0: leave all grad_*_sq at zero.

        return std::sqrt(grad_x_sq + grad_y_sq + grad_z_sq);
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

        Real gradientMagnitude = calcGradMagnitude(phiOld, i, j, k, sng, flags); // use this bc some kind of hyperbolic stuff and the cental difference is unstable

        phi(i, j, k) = phiOld(i, j, k) - dt_pseudo * sng * (gradientMagnitude - 1);
    }

    PYTHON()
    void reinitializeLevelset(Grid<Real> &phi, const FlagGrid &flags)
    {
        Grid<Real> phiOld(phi.getParent());
        Real dt_pseudo = 0.3;
        for (int __ = 0; __ < 10; __++)
        {
            phiOld.copyFrom(phi);
            knReinitializeLevelset(phi, phiOld, dt_pseudo, flags);
        }
    }

    // ---------- helpers: safe sampling & WENO5 1D derivative ----------
    inline Real safePhiSample(const Grid<Real> &phi, const FlagGrid &flags,
                              int i, int j, int k, int di, int dj, int dk)
    {
        const int I = i + di, J = j + dj, K = k + dk;
        if (!phi.isInBounds(Vec3i(I, J, K), 0) || flags.isObstacle(I, J, K))
            return phi(i, j, k);
        return phi(I, J, K);
    }

    // Compute WENO5 left/right-biased derivatives at cell center
    inline void weno5Derivatives1D(const Grid<Real> &phi, const FlagGrid &flags,
                                   int i, int j, int k, Real dx,
                                   int ax, // 0=x, 1=y, 2=z
                                   Real &Dm, Real &Dp)
    {
        auto S = [&](int s) -> Real
        {
            int di = (ax == 0) ? s : 0;
            int dj = (ax == 1) ? s : 0;
            int dk = (ax == 2) ? s : 0;
            return safePhiSample(phi, flags, i, j, k, di, dj, dk);
        };

        const Real v_m3 = S(-3), v_m2 = S(-2), v_m1 = S(-1), v_0 = S(0),
                   v_p1 = S(+1), v_p2 = S(+2), v_p3 = S(+3);

        bool poorStencil =
            (v_m3 == v_0 && v_m2 == v_0 && v_m1 == v_0) ||
            (v_p3 == v_0 && v_p2 == v_0 && v_p1 == v_0);

        if (poorStencil || dx <= 0)
        {
            Real back = (v_0 - v_m1) / (dx > 0 ? dx : 1.0);
            Real forw = (v_p1 - v_0) / (dx > 0 ? dx : 1.0);
            Dm = back;
            Dp = forw;
            return;
        }

        auto weno_minus = [&](Real vm2, Real vm1, Real v0, Real vp1, Real vp2) -> Real
        {
            Real p0 = (2.0 * vm2 - 7.0 * vm1 + 11.0 * v0) / 6.0;
            Real p1 = (-1.0 * vm1 + 5.0 * v0 + 2.0 * vp1) / 6.0;
            Real p2 = (2.0 * v0 + 5.0 * vp1 - 1.0 * vp2) / 6.0;

            Real b0 = (13.0 / 12.0) * ((vm2 - 2.0 * vm1 + v0) * (vm2 - 2.0 * vm1 + v0)) +
                      (1.0 / 4.0) * ((3.0 * vm2 - 4.0 * vm1 + v0) * (3.0 * vm2 - 4.0 * vm1 + v0));
            Real b1 = (13.0 / 12.0) * ((vm1 - 2.0 * v0 + vp1) * (vm1 - 2.0 * v0 + vp1)) +
                      (1.0 / 4.0) * ((vm1 - vp1) * (vm1 - vp1));
            Real b2 = (13.0 / 12.0) * ((v0 - 2.0 * vp1 + vp2) * (v0 - 2.0 * vp1 + vp2)) +
                      (1.0 / 4.0) * ((v0 - 4.0 * vp1 + 3.0 * vp2) * (v0 - 4.0 * vp1 + 3.0 * vp2));

            const Real eps = (Real)1e-6;
            const Real g0 = 0.1, g1 = 0.6, g2 = 0.3;
            Real a0 = g0 / ((eps + b0) * (eps + b0));
            Real a1 = g1 / ((eps + b1) * (eps + b1));
            Real a2 = g2 / ((eps + b2) * (eps + b2));
            Real sum = a0 + a1 + a2;

            return (a0 * p0 + a1 * p1 + a2 * p2) / sum;
        };

        auto weno_plus = [&](Real vm2, Real vm1, Real v0, Real vp1, Real vp2) -> Real
        {
            Real q0 = (-1.0 * vm2 + 5.0 * vm1 + 2.0 * v0) / 6.0;
            Real q1 = (2.0 * vm1 + 5.0 * v0 - 1.0 * vp1) / 6.0;
            Real q2 = (11.0 * v0 - 7.0 * vp1 + 2.0 * vp2) / 6.0;

            Real b0 = (13.0 / 12.0) * ((vm2 - 2.0 * vm1 + v0) * (vm2 - 2.0 * vm1 + v0)) +
                      (1.0 / 4.0) * ((vm2 - 4.0 * vm1 + 3.0 * v0) * (vm2 - 4.0 * vm1 + 3.0 * v0));
            Real b1 = (13.0 / 12.0) * ((vm1 - 2.0 * v0 + vp1) * (vm1 - 2.0 * v0 + vp1)) +
                      (1.0 / 4.0) * ((vm1 - vp1) * (vm1 - vp1));
            Real b2 = (13.0 / 12.0) * ((v0 - 2.0 * vp1 + vp2) * (v0 - 2.0 * vp1 + vp2)) +
                      (1.0 / 4.0) * ((3.0 * v0 - 4.0 * vp1 + vp2) * (3.0 * v0 - 4.0 * vp1 + vp2));

            const Real eps = (Real)1e-6;
            const Real g0 = 0.1, g1 = 0.6, g2 = 0.3;
            Real a0 = g0 / ((eps + b0) * (eps + b0));
            Real a1 = g1 / ((eps + b1) * (eps + b1));
            Real a2 = g2 / ((eps + b2) * (eps + b2));
            Real sum = a0 + a1 + a2;

            return (a0 * q0 + a1 * q1 + a2 * q2) / sum;
        };

        Real phi_imh_minus = weno_minus(v_m3, v_m2, v_m1, v_0, v_p1);
        Real phi_iph_minus = weno_minus(v_m2, v_m1, v_0, v_p1, v_p2);
        Real phi_imh_plus = weno_plus(v_p2, v_p1, v_0, v_m1, v_m2);
        Real phi_iph_plus = weno_plus(v_p3, v_p2, v_p1, v_0, v_m1);

        Dm = (phi_iph_minus - phi_imh_minus) / dx;
        Dp = (phi_iph_plus - phi_imh_plus) / dx;
    }

    inline void weno5DerivativesXYZ(const Grid<Real> &phi, const FlagGrid &flags,
                                    int i, int j, int k, Real dx,
                                    Real &Dmx, Real &Dpx,
                                    Real &Dmy, Real &Dpy,
                                    Real &Dmz, Real &Dpz)
    {
        weno5Derivatives1D(phi, flags, i, j, k, dx, 0, Dmx, Dpx);
        weno5Derivatives1D(phi, flags, i, j, k, dx, 1, Dmy, Dpy);
        if (phi.is3D())
            weno5Derivatives1D(phi, flags, i, j, k, dx, 2, Dmz, Dpz);
        else
        {
            Dmz = Dpz = 0.0;
        }
    }

    inline Real gradMagGodunovWENO5(const Grid<Real> &phi, int i, int j, int k,
                                    Real sng, const FlagGrid &flags, Real dx)
    {
        Real Dmx, Dpx, Dmy, Dpy, Dmz, Dpz;
        weno5DerivativesXYZ(phi, flags, i, j, k, dx, Dmx, Dpx, Dmy, Dpy, Dmz, Dpz);

        Real gx2 = 0, gy2 = 0, gz2 = 0;
        if (sng > 0.0)
        {
            gx2 = std::max((std::max(Dmx, (Real)0.0)) * (std::max(Dmx, (Real)0.0)),
                           (std::min(Dpx, (Real)0.0)) * (std::min(Dpx, (Real)0.0)));
            gy2 = std::max((std::max(Dmy, (Real)0.0)) * (std::max(Dmy, (Real)0.0)),
                           (std::min(Dpy, (Real)0.0)) * (std::min(Dpy, (Real)0.0)));
            gz2 = std::max((std::max(Dmz, (Real)0.0)) * (std::max(Dmz, (Real)0.0)),
                           (std::min(Dpz, (Real)0.0)) * (std::min(Dpz, (Real)0.0)));
        }
        else if (sng < 0.0)
        {
            gx2 = std::max((std::min(Dmx, (Real)0.0)) * (std::min(Dmx, (Real)0.0)),
                           (std::max(Dpx, (Real)0.0)) * (std::max(Dpx, (Real)0.0)));
            gy2 = std::max((std::min(Dmy, (Real)0.0)) * (std::min(Dmy, (Real)0.0)),
                           (std::max(Dpy, (Real)0.0)) * (std::max(Dpy, (Real)0.0)));
            gz2 = std::max((std::min(Dmz, (Real)0.0)) * (std::min(Dmz, (Real)0.0)),
                           (std::max(Dpz, (Real)0.0)) * (std::max(Dpz, (Real)0.0)));
        }

        return std::sqrt(gx2 + gy2 + gz2);
    }

    // ---------- kernel ----------
    KERNEL(bnd = 3) // WENO5 needs ±3 halo
    void knReinitializeLevelsetWENO(Grid<Real> &phi, Grid<Real> &phiOld,
                                    Real dt_pseudo, const FlagGrid &flags)
    {
        if (flags.isObstacle(i, j, k))
            return;

        const Real dx = 1.0;
        const Real p0 = phiOld(i, j, k);

        const Real eps = dx;
        const Real sng = p0 / std::sqrt(p0 * p0 + eps * eps);

        const Real gradMag = gradMagGodunovWENO5(phiOld, i, j, k, sng, flags, dx);

        phi(i, j, k) = phiOld(i, j, k) - dt_pseudo * sng * (gradMag - (Real)1.0);
    }

    PYTHON()
    void reinitializeLevelsetWENO(Grid<Real> &phi, const FlagGrid &flags)
    {
        Grid<Real> phiOld(phi.getParent());
        const Real dt_pseudo = 0.3;
        const int iters = 10;

        for (int n = 0; n < iters; ++n)
        {
            phiOld.copyFrom(phi);
            knReinitializeLevelset(phi, phiOld, dt_pseudo, flags);
        }
    }

    // PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES
    // PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES
    // PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES
    // PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES
    // PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES PARTICLES

    void addParticlesToCell(int i, int j, int k, Real phi, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, bool is3D, Real randomness = 0.5)
    {
        Real step = 1. / DISCRETIZATION;

        int alsdi = 0;

        for (int dk = 0; dk < (is3D ? DISCRETIZATION : 1); dk++)
        {
            for (int di = 0; di < DISCRETIZATION; di++)
            {
                for (int dj = 0; dj < DISCRETIZATION; dj++)
                {
                    Real randi = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;
                    Real randj = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;
                    Real randk = (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) - 0.5) * step * randomness;

                    Vec3 pos = Vec3(i + (0.5 + di) * step + randi, j + (0.5 + dj) * step + randj, k + (0.5 + dk) * step + randk);

                    particles.addBuffered(pos, phi > 0 ? ParticleBase::POUTSIDE : ParticleBase::PINSIDE);
                    alsdi++;
                }
            }
        }
    }

    KERNEL(points)
    void knSetTargetPhi(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii)
    {
        skipDelted(particles);

        Real randomNum = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

        randomNum = randomNum * (std::abs(particles.isInside(idx) ? NEGATIVE_SEED_CUTOFF : POSITIVE_SEED_CUTOFF) - MIN_RADIUS) + MIN_RADIUS;
        radii[idx] = (particles.isInside(idx) ? -1. : 1.) * randomNum;
    }

    Vec3 calcNormal(const Grid<Real> &phi, int i, int j, int k, const FlagGrid &flags, Real eps = (Real)1e-12)
    {
        // Compute normalized gradient (unit normal) of phi in 3D
        // Uses central differences where possible, else one-sided differences.

        Vec3i gs = flags.getParent()->getGridSize();
        Real dx = 1.0;

        auto valid = [&](int ii, int jj, int kk)
        {
            bool inBounds = 0 <= i && i < gs.x && 0 <= j && j < gs.y && 0 <= k && k < gs.z;
            return inBounds && !flags.isObstacle(ii, jj, kk);
        };

        Real dphix = 0.0, dphiy = 0.0, dphiz = 0.0;

        // ----- X derivative -----
        if (valid(i - 1, j, k) && valid(i + 1, j, k))
        {
            dphix = (phi(i + 1, j, k) - phi(i - 1, j, k)) * 0.5 / dx; // central
        }
        else if (valid(i + 1, j, k))
        { // forward
            dphix = (phi(i + 1, j, k) - phi(i, j, k)) / dx;
        }
        else if (valid(i - 1, j, k))
        { // backward
            dphix = (phi(i, j, k) - phi(i - 1, j, k)) / dx;
        }

        // ----- Y derivative -----
        if (valid(i, j - 1, k) && valid(i, j + 1, k))
        {
            dphiy = (phi(i, j + 1, k) - phi(i, j - 1, k)) * 0.5 / dx; // central
        }
        else if (valid(i, j + 1, k))
        {
            dphiy = (phi(i, j + 1, k) - phi(i, j, k)) / dx;
        }
        else if (valid(i, j - 1, k))
        {
            dphiy = (phi(i, j, k) - phi(i, j - 1, k)) / dx;
        }

        if (phi.is3D())
        {
            // ----- Z derivative -----
            if (valid(i, j, k - 1) && valid(i, j, k + 1))
            {
                dphiz = (phi(i, j, k + 1) - phi(i, j, k - 1)) * 0.5 / dx; // central
            }
            else if (valid(i, j, k + 1))
            {
                dphiz = (phi(i, j, k + 1) - phi(i, j, k)) / dx;
            }
            else if (valid(i, j, k - 1))
            {
                dphiz = (phi(i, j, k) - phi(i, j, k - 1)) / dx;
            }
        }

        Real len = std::sqrt(dphix * dphix + dphiy * dphiy + dphiz * dphiz) + eps;
        return Vec3(dphix / len, dphiy / len, dphiz / len);
    }

    KERNEL(points)
    void knAttractionStep(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi, const FlagGrid &flags)
    {
        skipDelted(particles);

        Vec3 pos = particles.getPos(idx);

        Real lambda = 1;
        Vec3 newPos = pos + lambda * (radii[idx] - phi.getInterpolatedHi(pos, 2)) * calcNormal(phi, std::floor(pos.x), std::floor(pos.y), std::floor(pos.z), flags);

        particles.setPos(idx, newPos);
    }

    KERNEL(points)
    void knSetParticleRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        skipDelted(particles);

        Real phiIdx = std::abs(phi.getInterpolatedHi(particles[idx].pos, 2));

        radii[idx] = Manta::clamp(std::abs(phiIdx), MIN_RADIUS, MAX_RADIUS);
    }

    PYTHON()
    void reinitializeRadii(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, const Grid<Real> &phi)
    {
        knSetParticleRadii(particles, radii, phi);
    }

    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET
    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET
    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET
    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET
    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET
    // PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET PARTICLE LEVEL SET

    inline Real particleSphere(BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, IndexInt idx, Vec3 x)
    {
        return (particles.isInside(idx) ? -1.0 : 1.0) * (radii[idx] - norm(x - particles[idx].pos));
    }

    KERNEL()
    void knUpdateOmega(Grid<Real> &omega, std::vector<IndexInt> &escapedParticles, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii, bool useMin)
    {
        for (auto idx : escapedParticles)
        {
            if (useMin)
            {
                omega(i, j, k) = std::min(omega(i, j, k), particleSphere(particles, radii, idx, Vec3(i + 0.5, j + 0.5, k + 0.5)));
            }
            else
            {
                omega(i, j, k) = std::max(omega(i, j, k), particleSphere(particles, radii, idx, Vec3(i + 0.5, j + 0.5, k + 0.5)));
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
            if (particles.isDeleted(idx))
            {
                continue;
            }

            Vec3 pos = particles[idx].pos;
            Real phiAtPos = phi.getInterpolatedHi(pos, 2);

            bool isInside = (phiAtPos <= 0);
            if ((isInside && particles.isInside(idx)) || (!isInside && particles.isOutside(idx)))
            {
                particles[idx].flag &= ~ParticleBase::PESCAPED;
                continue; // on the right side
            }

            if (std::abs(phiAtPos) < ESCAPE_CONDITION * radii[idx])
            {
                continue; // particle still not far enough on wrong side
            }

            particles[idx].flag |= ParticleBase::PESCAPED;
            if (particles.isInside(idx))
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

        // std::cout << "POSITIVE ESCAPED: " << positiveEscaped.size() << std::endl;
        // std::cout << "NEGATIVE ESCAPED: " << negativeEscaped.size() << std::endl;

        // Step 2: correct omega+ and omega-
        if (positiveEscaped.size() == 0 && negativeEscaped.size() == 0)
        {
            return;
        }

        Grid<Real> omegaPlus(phi.getParent());
        Grid<Real> omegaMinus(phi.getParent());
        omegaPlus.copyFrom(phi);
        omegaMinus.copyFrom(phi);

        knUpdateOmega(omegaPlus, positiveEscaped, particles, radii, false);
        knUpdateOmega(omegaMinus, negativeEscaped, particles, radii, true);

        // Step 3: combine omage+ and omega-
        knCombineOmegas(omegaPlus, omegaMinus, phi, flags);
    }

    struct ParticleHeapEntry
    {
        Real metric;
        IndexInt idx;
        bool operator<(const ParticleHeapEntry &other) const
        {
            return metric < other.metric; // reverse: largest on top
        }
    };

    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON
    // PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON PYTHON

    PYTHON()
    void sampleLevelsetBorderWithParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii)
    {
        srand(RANDOM_SEED);

        FOR_IJK_BND(flags, 0)
        {
            if (flags.isObstacle(i, j, k))
            {
                continue;
            }
            if (NEGATIVE_SEED_CUTOFF < phi(i, j, k) && phi(i, j, k) < POSITIVE_SEED_CUTOFF)
            {
                addParticlesToCell(i, j, k, phi(i, j, k), particles, radii, phi.is3D());
            }
        }

        particles.insertBufferedParticles();

        for (IndexInt idx = 0, max = particles.size(); idx < particles.size(); idx++) // do this for a repeatable random seed, with threads not random
        {
            if (particles[idx].flag & ParticleBase::PDELETE)
            {
                continue;
            }

            Real randomNum = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

            randomNum = randomNum * (std::abs(particles.isInside(idx) ? NEGATIVE_SEED_CUTOFF : POSITIVE_SEED_CUTOFF) - MIN_RADIUS) + MIN_RADIUS;
            radii[idx] = (particles.isInside(idx) ? -1. : 1.) * randomNum;
        }
        // knSetTargetPhi(particles, radii); // store phi_target in radii for now, later the actual radii

        for (int __ = 0; __ < 10; __++)
        {
            knAttractionStep(particles, radii, phi, flags);
        }

        knSetParticleRadii(particles, radii, phi);
    }

    PYTHON()
    void advectParticleLevelSet(Grid<Real> *phi, BasicParticleSystem *particles, ParticleDataImpl<Real> *radii, const MACGrid *vel, const FlagGrid *flags)
    {
        // Time integration
        simpleSLAdvect(flags, vel, phi, 2, true, NORMAL);
        advectParticlesForward(particles, vel, flags);

        // Error correction
        correctErrorsWithParticles(*phi, *particles, *radii, *flags);

        // Reinitialization
        reinitializeLevelset(*phi, *flags);
        correctErrorsWithParticles(*phi, *particles, *radii, *flags);

        // Adjust Radii
        knSetParticleRadii(*particles, *radii, *phi);
    }

    PYTHON()
    void reseedParticles(const Grid<Real> &phi, const FlagGrid &flags, BasicParticleSystem &particles, ParticleDataImpl<Real> &radii)
    {
        Vec3i gs = phi.getParent()->getGridSize();
        std::vector<std::vector<IndexInt>> particlesInEachCell(gs.x * gs.y * gs.z);
        IndexInt strideZ = gs.z > 1 ? (gs.x * gs.y) : 0;
        IndexInt particlesPerCell = DISCRETIZATION * DISCRETIZATION * (phi.is3D() ? DISCRETIZATION : 1);

        int deletedParticles = 0;
        int addedParticles = 0;

        auto index = [&](IndexInt x, IndexInt y, IndexInt z)
        {
            return (IndexInt)x + (IndexInt)gs.x * y + (IndexInt)strideZ * z;
        };

        // Step 1: delete particles from not-important regions and count not-escaped particles in border-regions
        for (IndexInt idx = 0, max = particles.size(); idx < particles.size(); idx++)
        {
            if (particles.isDeleted(idx))
            {
                continue;
            }

            if (particles.isEscaped(idx))
            {
                continue; // don't delete escaped particles
            }

            int i = std::floor(particles[idx].pos[0]);
            int j = std::floor(particles[idx].pos[1]);
            int k = std::floor(particles[idx].pos[2]);

            // Negative particle too far inside the fluid
            if (particles.isInside(idx) && phi.getInterpolatedHi(particles.getPos(idx), 2) < NEGATIVE_SEED_CUTOFF)
            {
                particles.kill(idx);
                deletedParticles++;
                continue;
            }

            // Positive particle too far outside the fluid
            if (particles.isOutside(idx) && phi.getInterpolatedHi(particles.getPos(idx), 2) > POSITIVE_SEED_CUTOFF)
            {
                particles.kill(idx);
                deletedParticles++;
                continue;
            }

            particlesInEachCell[index(i, j, k)].push_back(idx);
        }

        // Step 2: Delete particles from cells that have to many.
        FOR_IJK(phi)
        {
            std::vector<IndexInt> &particleIndexList = particlesInEachCell[index(i, j, k)];

            if ((int)particleIndexList.size() > particlesPerCell)
            {
                std::priority_queue<ParticleHeapEntry> heap;

                // Build initial heap of desired size
                for (IndexInt idx : particleIndexList)
                {
                    Real metric = (particles.isOutside(idx) ? 1.0 : -1.0) * phi.getInterpolatedHi(particles.getPos(idx), 2) - radii[idx];

                    heap.push({metric, idx});
                }

                // Kill the "worst" perfoming particles, i.e. top-of-heap particles
                while ((int)heap.size() > particlesPerCell)
                {
                    deletedParticles++;
                    particles.kill(heap.top().idx);
                    heap.pop();
                    deletedParticles++;
                }
            }
        }

        // Step 3: insert particles to cells that don't have enough
        FOR_IJK(phi)
        {
            if (flags.isObstacle(i, j, k))
            {
                continue;
            }

            if (NEGATIVE_SEED_CUTOFF <= phi(i, j, k) && phi(i, j, k) <= POSITIVE_SEED_CUTOFF && particlesInEachCell[index(i, j, k)].size() < particlesPerCell)
            {
                for (int _ = 0; _ < particlesPerCell - particlesInEachCell[index(i, j, k)].size(); _++)
                {
                    Real randi = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
                    Real randj = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
                    Real randz = phi.is3D() ? static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX) : 0.0;

                    Vec3 pos = Vec3(i + randi, j + randj, k + randz);

                    particles.addBuffered(pos, phi(i, j, k) > 0 ? ParticleBase::POUTSIDE : ParticleBase::PINSIDE);
                    addedParticles++;
                }
            }
        }
        particles.insertBufferedParticles();
        knSetParticleRadii(particles, radii, phi);

        std::cout << "Deleted " << deletedParticles << " particles." << std::endl;
        std::cout << "Created " << addedParticles << " particles." << std::endl;
    }

    // OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS
    // OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS
    // OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS
    // OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS
    // OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS OTHER FUNCTIONS

    KERNEL()
    void knSetFlagsFromParticleLevelset(const Grid<Real> &phi, FlagGrid &flags, Real level)
    {
        if (!flags.isObstacle(i, j, k))
        {
            flags(i, j, k) = 0;

            if (phi(i, j, k) < level)
            {
                flags(i, j, k) = FlagGrid::TypeFluid;
            }
            else
            {
                flags(i, j, k) = FlagGrid::TypeEmpty;
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
    void knMarkPhiFromFlagGrid(Grid<Real> &phi, const FlagGrid &flags)
    {
        phi(i, j, k) = flags.isFluid(i, j, k) ? -1. : (flags.isObstacle(i, j, k) ? 0. : 1.);
    }

    PYTHON()
    void markPhiFromFlagGrid(Grid<Real> &phi, const FlagGrid &flags)
    {
        knMarkPhiFromFlagGrid(phi, flags);
    }

    // VELOCITY EXTRAPOLATION VELOCITY EXTRAPOLATION VELOCTY EXTRAPOLAION
    // VELOCITY EXTRAPOLATION VELOCITY EXTRAPOLATION VELOCTY EXTRAPOLAION
    // VELOCITY EXTRAPOLATION VELOCITY EXTRAPOLATION VELOCTY EXTRAPOLAION
    // VELOCITY EXTRAPOLATION VELOCITY EXTRAPOLATION VELOCTY EXTRAPOLAION
    // VELOCITY EXTRAPOLATION VELOCITY EXTRAPOLATION VELOCTY EXTRAPOLAION

    // FAST SWEEPING MARCHING EXTRAPOLATION FAST SWEEPING MARCHING EXTRAPOLATION
    // FAST SWEEPING MARCHING EXTRAPOLATION FAST SWEEPING MARCHING EXTRAPOLATION
    // FAST SWEEPING MARCHING EXTRAPOLATION FAST SWEEPING MARCHING EXTRAPOLATION
    // FAST SWEEPING MARCHING EXTRAPOLATION FAST SWEEPING MARCHING EXTRAPOLATION
    // FAST SWEEPING MARCHING EXTRAPOLATION FAST SWEEPING MARCHING EXTRAPOLATION

    KERNEL()
    void knFillLocked(const FlagGrid &flags, Grid<int> &locked, Grid<Real> &velComponent, MACGridComponent component)
    {
        const Real D_INF = std::numeric_limits<Real>::max();
        bool isSolid;
        bool isFluid;

        // direction: 0=x, 1=y, 3=z
        switch (component)
        {
        case MAC_X:
            isSolid = flags.isObstacle(i - 1, j, k) || flags.isObstacle(i, j, k);
            isFluid = flags.isFluid(i - 1, j, k) && flags.isFluid(i, j, k);
            break;
        case MAC_Y:
            isSolid = flags.isObstacle(i, j - 1, k) || flags.isObstacle(i, j, k);
            isFluid = flags.isFluid(i, j - 1, k) && flags.isFluid(i, j, k);
            break;
        case MAC_Z:
            throw std::runtime_error("not implemented yet! (getPhiAtMACCell with MAC_Z)");
        }

        if (isSolid || isFluid)
        {
            locked(i, j, k) = true;
            if (isSolid)
            {
                velComponent(i, j, k) = 0.;
            }
        }
        else
        {
            locked(i, j, k) = false;
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
        /* if (i == 1 || j == 1){
            grid(i , j, k) = 10;
            return;
        } */
        grid(i, j, k) = grid(i, j, k) == std::numeric_limits<Real>::max() ? 0. : grid(i, j, k);
    }

    void extrapolateComponentFSM(Grid<Real> &velComponent, const Grid<Real> &phi, Grid<int> &locked, int sweepIterations, Vec3i gs, MACGridComponent component)
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
                    if (locked(i, j, k))
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

                    // std::cout << "sum: " << sumW << std::endl;

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
                    if (locked(i, j, k))
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
                    if (locked(i, j, k))
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
                    if (locked(i, j, k))
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

        Grid<int> locked(parent);

        // X-component
        knFillLocked(flags, locked, velX, MAC_X);
        extrapolateComponentFSM(velX, phi, locked, steps, gs, MAC_X);

        // Y-component
        knFillLocked(flags, locked, velY, MAC_Y);
        extrapolateComponentFSM(velY, phi, locked, steps, gs, MAC_Y);

        knGrids2MAC(vel, velX, velY, velZ);
    }

    // MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION
    // MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION
    // MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION
    // MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION
    // MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION MISCELLANEOUS FUNCTION

    KERNEL()
    void knSetAllEmptyFlagsToLiquid(FlagGrid &flags, int _)
    {
        if (flags.isEmpty(i, j, k))
        {
            flags(i, j, k) = (flags(i, j, k) | FlagGrid::TypeFluid) & ~FlagGrid::TypeEmpty;
        }
    }

    PYTHON()
    void setAllEmptyFlagsToLiquid(FlagGrid &flags)
    {
        knSetAllEmptyFlagsToLiquid(flags, 1);
    }

    KERNEL()
    void knSetFlagsFromDensity(FlagGrid &flags, const Grid<Real> &density, Real level)
    {
        if (!flags.isObstacle(i, j, k))
        {
            flags(i, j, k) = 0;

            if (density(i, j, k) > level)
            {
                flags(i, j, k) = FlagGrid::TypeFluid;
            }
            else
            {
                flags(i, j, k) = FlagGrid::TypeEmpty;
            }
        }
    }

    PYTHON()
    void setFlagsFromDensity(FlagGrid &flags, const Grid<Real> &density, Real level)
    {
        knSetFlagsFromDensity(flags, density, level);
    }

    KERNEL()
    void knVisualizeFlags(const FlagGrid &flags, Grid<Real> &grid, const FlagGrid *flags_n_plus_one, bool onlyFluid)
    {
        if (onlyFluid)
        {
            if (flags.isFluid(i, j, k) || (flags_n_plus_one && flags_n_plus_one->isFluid(i, j, k)))
            {
                grid(i, j, k) = 1;
            }
            else
            {
                grid(i, j, k) = 0;
            }
            return;
        }

        if (flags_n_plus_one)
        {
            if (flags.isFluid(i, j, k))
            {
                grid(i, j, k) = -1;
            }
            if (flags_n_plus_one->isFluid(i, j, k) && !flags.isFluid(i, j, k))
            {
                grid(i, j, k) = -100;
            }
            else if (flags.isEmpty(i, j, k))
            {
                grid(i, j, k) = 0;
            }
            else if (flags.isObstacle(i, j, k))
            {
                grid(i, j, k) = 1;
            }
            return;
        }
        if (flags.isEmpty(i, j, k))
        {
            grid(i, j, k) = 0;
        }
        else if (flags.isObstacle(i, j, k))
        {
            grid(i, j, k) = 1;
        }
        else if (flags.isFluid(i, j, k))
        {
            grid(i, j, k) = -1;
        }
        else
        {
            grid(i, j, k) = -100;
        }
    }

    PYTHON()
    void visualizeFlags(const FlagGrid &flags, Grid<Real> &grid, const FlagGrid *flags_n_plus_one = nullptr, bool onlyFluid = false)
    {
        knVisualizeFlags(flags, grid, flags_n_plus_one, onlyFluid);
    }
}