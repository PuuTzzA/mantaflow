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
#include <sstream>
#include <iomanip>

using namespace std;
namespace Manta
{
    PYTHON()
    std::string calculateTotalMassOfGrid(const Grid<Real> *grid)
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

    inline bool doCurlHere(IndexInt i, IndexInt j, IndexInt k, const FlagGrid &flags)
    {
        bool inbounds = i >= 1 && j >= 1 && k >= 1 && i < flags.getSize().x - 1 && j < flags.getSize().y - 1 && k < flags.getSize().z - 1;
        return inbounds && !flags.isObstacle(i, j, k);
    }

    KERNEL()
    void knCalculateCurl(const MACGrid &vel, Grid<Real> &curl, const FlagGrid &flags)
    {
        if (!curl.is3D())
        {

            curl(i, j, k) = (vel.getCentered(i + 1, j, k).y - vel.getCentered(i - 1, j, k).y) / 2;
            curl(i, j, k) -= (vel.getCentered(i, j + 1, k).y - vel.getCentered(i, j - 1, k).y) / 2;
        }
        else
        {
            if (!doCurlHere(i, j, k, flags))
            {
                curl(i, j, k) = 0.;
                return;
            }

            Vec3 curlVec(0.);

            curlVec.x = (vel.getCentered(i, j + 1, k).z - vel.getCentered(i, j - 1, k).z) / 2;
            curlVec.x -= (vel.getCentered(i, j, k + 1).y - vel.getCentered(i, j, k - 1).y) / 2;

            curlVec.y = (vel.getCentered(i, j, k + 1).x - vel.getCentered(i, j, k - 1).x) / 2;
            curlVec.y -= (vel.getCentered(i + 1, j, k).z - vel.getCentered(i - 1, j, k).z) / 2;

            curlVec.z = (vel.getCentered(i + 1, j, k).y - vel.getCentered(i - 1, j, k).y) / 2;
            curlVec.z -= (vel.getCentered(i, j + 1, k).x - vel.getCentered(i, j - 1, k).x) / 2;

            curl(i, j, k) = normalize(curlVec);
        }
    }

    PYTHON()
    void calculateCurl(const MACGrid &vel, Grid<Real> &curl, const FlagGrid &flags)
    {
        knCalculateCurl(vel, curl, flags);
    }

    PYTHON()
    std::string realGridStats(const Grid<Real> &grid, const FlagGrid &flags)
    {
        double sum{0};
        Real min = std::numeric_limits<Real>::infinity();
        Real max = -std::numeric_limits<Real>::infinity();

        double sumOnlyFluid{0};
        IndexInt numFluidCells = 0;
        Real minOnlyFluid = std::numeric_limits<Real>::infinity();
        Real maxOnlyFluid = -std::numeric_limits<Real>::infinity();

        FOR_IJK(grid)
        {
            Real val = grid(i, j, k);
            sum += val;
            min = std::min(min, val);
            max = std::max(max, val);

            if (flags.isFluid(i, j, k))
            {
                numFluidCells++;
                sumOnlyFluid += val;
                minOnlyFluid = std::min(minOnlyFluid, val);
                maxOnlyFluid = std::max(maxOnlyFluid, val);
            }
        }

        Vec3i gs = grid.getParent()->getGridSize();
        const int n = gs.x * gs.y * gs.z;

        double mean = sum / n;
        double meanOnlyFluid = sumOnlyFluid / numFluidCells;

        double squaredDif{0};
        double squaredDifOnlyFluid{0};

        FOR_IJK(grid)
        {
            Real val = grid(i, j, k);

            squaredDif += (val - mean) * (val - mean);

            if (flags.isFluid(i, j, k))
            {
                squaredDifOnlyFluid += (val - meanOnlyFluid) * (val - meanOnlyFluid);
            }
        }

        double std = std::sqrt(squaredDif / n);
        double stdOnlyFluid = std::sqrt(squaredDifOnlyFluid / numFluidCells);

        std::ostringstream json;
        json.setf(std::ios::fixed);
        json << std::setprecision(6); // control decimal places

        json << '{'
             << "\"cells\":" << n << ','
             << "\"sum\":" << sum << ','
             << "\"min\":" << min << ','
             << "\"max\":" << max << ','
             << "\"mean\":" << mean << ','
             << "\"std\":" << std << ','
             << "\"fluidCells\":" << numFluidCells << ','
             << "\"sumFluid\":" << sumOnlyFluid << ','
             << "\"minFluid\":" << (numFluidCells ? minOnlyFluid : 0.0) << ','
             << "\"maxFluid\":" << (numFluidCells ? maxOnlyFluid : 0.0) << ','
             << "\"meanFluid\":" << meanOnlyFluid << ','
             << "\"stdFluid\":" << stdOnlyFluid
             << '}';

        return json.str(); // ready to log, save, or send
    }

    //! Kernel: Compute max norm of vec grid
    KERNEL(ijk, reduce = max)
    returns(Real maxVal = -std::numeric_limits<Real>::max())
        Real knGetMaxVal(const Grid<Real> &grid, const FlagGrid &flags)
    {
        maxVal = std::max(maxVal, grid(i, j, k));
    }

    PYTHON()
    Real getMaxVal(const Grid<Real> &grid, const FlagGrid &flags)
    {
        return knGetMaxVal(grid, flags);
    }

    KERNEL(bnd = 2)
    void knComputeVelocityMagnitude(Grid<Real> &dest, const MACGrid &vel)
    {
        dest(i, j, k) = norm(vel.getCentered(i, j, k));
    }

    PYTHON()
    void computeVelocityMagnitude(Grid<Real> &dest, const MACGrid &vel)
    {
        knComputeVelocityMagnitude(dest, vel);
    }

    KERNEL()
    void knStoreVelocityMagnitude(Grid<Real> &dest, const MACGrid &vel)
    {
        dest(i, j, k) = norm(vel(i, j, k));
    }

    PYTHON()
    void storeVelocityMagnitude(Grid<Real> &dest, const MACGrid &vel)
    {
        knStoreVelocityMagnitude(dest, vel);
    }

    PYTHON()
    double calculateRelativeError(Grid<Real> &phi0, Grid<Real> &phin)
    {
        double sum0squared{0.};
        double differenceSquared{0.};

        FOR_IJK(phi0)
        {
            sum0squared += phi0(i, j, k) * phi0(i, j, k);
            differenceSquared += (phin(i, j, k) - phi0(i, j, k)) * (phin(i, j, k) - phi0(i, j, k));
        }

        return std::sqrt(differenceSquared / sum0squared);
    }
}