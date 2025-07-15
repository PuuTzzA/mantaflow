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

    KERNEL()
    void knCalculateCurl(const MACGrid &vel, Grid<Real> &curl)
    {
        curl(i, j, k) = (vel.getCentered(i + 1, j, k).y - vel.getCentered(i - 1, j, k).y) / 2;
        curl(i, j, k) -= (vel.getCentered(i, j + 1, k).y - vel.getCentered(i, j - 1, k).y) / 2;
    }

    PYTHON()
    void calculateCurl(const MACGrid &vel, Grid<Real> &curl)
    {
        knCalculateCurl(vel, curl);
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
}