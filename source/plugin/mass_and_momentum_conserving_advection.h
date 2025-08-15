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
        LINEAR,
        CUBIC_CONVOLUTIONAL,    // standart cubic interpolation (less accurate than CUBIC_POLYNOMIAL)
        CUBIC_POLYNOMIAL,       // polynomial interpolation with Lagrange basis functions
        MONOTONE_CUBIC_HERMITE, // Monotone cubic Hermite interpolation (nach Fritsch und Carlson)
    };

    enum TargetCellType
    {
        NOT_OBSTACLE,
        FLUID_ISH,
        FLUID_STRICT
    };

    inline std::string toString(MACGridComponent component)
    {
        switch (component)
        {
        case MAC_X:
            return "MAC_X";
        case MAC_Y:
            return "MAC_Y";
        case MAC_Z:
            return "MAC_Z";
        case NONE:
            return "NONE";
        default:
            return "???";
        }
    }

    inline std::string toString(InterpolationType component)
    {
        switch (component)
        {
        case LINEAR:
            return "Linear";
        case CUBIC_CONVOLUTIONAL:
            return "Cubic Convolutional";
        case CUBIC_POLYNOMIAL:
            return "Cubic Polynomial";
        case MONOTONE_CUBIC_HERMITE:
            return "Monotone Cubic Hermite";
        default:
            return "???";
        }
    }

    void fillWithOnes(GridBase *grid);

    void fillFluidWithOnes(GridBase *grid, const FlagGrid *flags);

    bool isValidFluid(IndexInt i, IndexInt j, IndexInt k, const FlagGrid &flags, MACGridComponent component);

    KERNEL()
    void knMAC2Grids(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
    {
        Vec3 data = vel(i, j, k);
        velX(i, j, k) = data.x;
        velY(i, j, k) = data.y;
        velZ(i, j, k) = data.z;
    }

    KERNEL()
    void knGrids2MAC(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
    {
        vel(i, j, k) = Vec3(velX(i, j, k), velY(i, j, k), velZ(i, j, k));
    }

    void simpleSLAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, int interpolationType, bool all);

    void advectParticlesForward(BasicParticleSystem *particles, const MACGrid *vel, const FlagGrid *flags);
}