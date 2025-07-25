#include "grid.h"
#include <cmath>
#include <unordered_map>
#include <string>
#include <functional>

namespace Manta
{
    using velFunction = std::function<Vec3(IndexInt, IndexInt, IndexInt)>;

#define M 0.001
#define U0 1
#define W 3

    KERNEL()
    void knSetVelocityField(MACGrid &vel, const FlagGrid &flags, velFunction function)
    {
        if (!flags.isObstacle(i, j, k))
        {
            vel(i, j, k) = function(i, j, k);
        }
    }

    PYTHON()
    void setVelocityField(MACGrid &vel, const FlagGrid &flags, std::string functionName, Real t = 0.)
    {
        std::unordered_map<std::string, velFunction> functionMap;

        functionMap["diagonal_motion"] = [](IndexInt i, IndexInt j, IndexInt k)
        { return Vec3(U0); };

        functionMap["zalesak_rotation"] = [&flags](IndexInt i, IndexInt j, IndexInt k)
        {
            Vec3i _gs = flags.getParent()->getGridSize();
            Vec3 gs = Vec3(Real(_gs.x), Real(_gs.y), 0.);

            Real u = -U0 * M_PI * ((j + 0.5) / gs.x - 0.5);
            Real v = U0 * M_PI * ((i + 0.5) / gs.y - 0.5);
            Real w = 0.;
            return Vec3(u, v, w);
        };

        functionMap["shear_flow"] = [&flags](IndexInt i, IndexInt j, IndexInt k)
        {
            Vec3i _gs = flags.getParent()->getGridSize();
            Vec3 gs = Vec3(Real(_gs.x), Real(_gs.y), 0.);

            Real u = -U0 * M_PI * std::cos(M_PI * (i / gs.x - 0.5)) * std::sin(M_PI * ((j + 0.5) / gs.x - 0.5));
            Real v = U0 * M_PI * std::sin(M_PI * ((i + 0.5) / gs.y - 0.5)) * std::cos(M_PI * (j / gs.y - 0.5));
            Real w = 0.;
            return Vec3(u, v, w);
        };

        functionMap["deformation_field"] = [&flags](IndexInt i, IndexInt j, IndexInt k)
        {
            Vec3i _gs = flags.getParent()->getGridSize();
            Vec3 gs = Vec3(Real(_gs.x), Real(_gs.y), 0.);

            Real u = -U0 * std::sin(4. * M_PI * (i / gs.x + 0.5)) * std::sin(4. * M_PI * ((j + 0.5) / gs.x + 0.5));
            Real v = -U0 * std::cos(4. * M_PI * ((i + 0.5) / gs.x + 0.5)) * std::cos(4. * M_PI * (j / gs.x + 0.5));
            Real w = 0.;
            return Vec3(u, v, w);
        };

        knSetVelocityField(vel, flags, functionMap[functionName]);
    }
}