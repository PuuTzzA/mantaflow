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
    enum MACGridComponent
    {
        MAC_X,
        MAC_Y,
        MAC_Z,
        NONE
    };

    enum InterpolationType
    {
        TRILIENAR,
        CATMULL_ROM
    };

    //! Semi-Lagrange interpolation kernel
    void
    fillWithOnes(GridBase *grid);

    void fillFluidWithOnes(GridBase *grid, const FlagGrid *flags);

    bool isValid(int i, int j, int k, const FlagGrid &flags, Vec3i &gs);

    inline Vec3 RK4(Vec3 pos, Real dt, const MACGrid &vel);

    template <class GridType>
    void fnMassMomentumConservingAdvectUnified(FluidSolver *parent, const FlagGrid &flags_n, const FlagGrid &flags_n_plus_one, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset, const Grid<Real> *phi = nullptr, MACGridComponent component = NONE);

    KERNEL()
    void knMAC2Grids(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
    {
        Vec3 data = vel(i, j, k);
        velX(i, j, k) = data.x;
        velY(i, j, k) = data.y;
        velZ(i, j, k) = data.z;
    }

    KERNEL()
    void knGrids2MAC(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ, const FlagGrid &flags)
    {
        vel(i, j, k) = Vec3(velX(i, j, k), velY(i, j, k), velZ(i, j, k));
    }

    void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const FlagGrid &flags_n_plus_one, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative, bool water);

    void massMomentumConservingAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative);

    void massMomentumConservingAdvectWater(const FlagGrid *flags_n, const FlagGrid *flags_n_plus_one, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative);

    std::string calculateMass(const Grid<Real> *grid);

    void advectParticlesForward(BasicParticleSystem *particles, const MACGrid *vel, const FlagGrid *flags);

    void simpleSLAdvection(const FlagGrid *flags, const MACGrid *vel, Grid<Real> *grid);

}