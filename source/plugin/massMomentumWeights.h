#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
#include <unordered_map>
#include <unordered_set>

namespace Manta
{
    class MassMomentumWeights
    {
    public:
        MassMomentumWeights(Vec3i gridSize);

        void insert(Vec3i cellI, Vec3i cellJ, Real value);

        void add(Vec3i cellI, Vec3i cellJ, Real value);

        void scale(Vec3i cellI, Vec3i cellJ, Real factor);

        void scaleAllWeightsAt(Vec3i cellI, Real factor);

        void scaleAllReverseWeightsAt(Vec3i cellJ, Real factor);

        void calculateBeta(Grid<Real> &beta);

        void calculateGamma(Grid<Real> &gamma);

        void calculateIntermediateResult(Grid<Real> &dest, Grid<Real> &src, Grid<Real> &min, Grid<Real> &max);

        void distributeLostMass(Grid<Real>&grid, Grid<Real> &lostMass, Grid<Real>& min, Grid<Real>& max);

        inline IndexInt index(Vec3i pos) const
        {
            return (IndexInt)pos.x + (IndexInt)gs.x * pos.y + (IndexInt)strideZ * pos.z;
        };

    protected:
        Vec3i gs;
        IndexInt strideZ;
        std::vector<std::vector<std::pair<Vec3i, Real>>> weights;
        std::vector<std::vector<Vec3i>> reverseWeights;

        Real &getOrInsert(std::vector<std::pair<Vec3i, Real>> &vec, Vec3i key, Real defaultValue = 0.0);
    };
}
