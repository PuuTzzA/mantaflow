#include "vectorbase.h"
#include "grid.h"
#include "particle.h"
#include "kernel.h"
#include "massMomentumWeights.h"
#include <unordered_map>
#include <unordered_set>

namespace std
{
    template <>
    struct hash<Manta::Vector3D<int>>
    {
        std::size_t operator()(const Manta::Vector3D<int> &v) const noexcept
        {
            std::size_t h1 = std::hash<int>{}(v.x);
            std::size_t h2 = std::hash<int>{}(v.y);
            std::size_t h3 = std::hash<int>{}(v.z);
            // Combine hashes (this is a common technique)
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

namespace Manta
{
    MassMomentumWeights::MassMomentumWeights(Vec3i gridSize) : gs{gridSize}, weights(gridSize.x * gridSize.y * gridSize.z), reverseWeights(gridSize.x * gridSize.y * gridSize.z)
    {
        int amount = gs.z > 1 ? 10 : 5;
        for (int i = 0; i < gs.x * gs.y * gs.z; i++)
        {
            weights[i].reserve(amount);
            reverseWeights[i].reserve(amount);
        }

        strideZ = gs.z > 1 ? (gs.x * gs.y) : 0;
    };

    Real &MassMomentumWeights::getOrInsert(std::vector<std::pair<Vec3i, Real>> &vec, Vec3i key, Real defaultValue)
    {
        for (auto &[k, v] : vec)
        {
            if (k == key)
                return v;
        }
        vec.emplace_back(key, defaultValue);
        return vec.back().second;
    }

    void MassMomentumWeights::insert(Vec3i cellI, Vec3i cellJ, Real value)
    {
        getOrInsert(weights[index(cellI)], cellJ) = value;
        std::vector<Vec3i> &rev = reverseWeights[index(cellJ)];
        if (std::find(rev.begin(), rev.end(), cellI) == rev.end())
        {
            rev.push_back(cellI);
        }
    }

    void MassMomentumWeights::add(Vec3i cellI, Vec3i cellJ, Real value)
    {
        getOrInsert(weights[index(cellI)], cellJ) += value;
        std::vector<Vec3i> &rev = reverseWeights[index(cellJ)];
        if (std::find(rev.begin(), rev.end(), cellI) == rev.end())
        {
            rev.push_back(cellI);
        }
    }

    void MassMomentumWeights::scale(Vec3i cellI, Vec3i cellJ, Real factor)
    {
        for (auto &[_cellJ, w] : weights[index(cellI)])
        {
            if (_cellJ == cellJ)
            {
                w *= factor;
                return;
            }
        }
    }

    void MassMomentumWeights::scaleAllWeightsAt(Vec3i cellI, Real factor)
    {
        for (auto &[_, value] : weights[index(cellI)])
        {
            value *= factor;
        }
    }

    void MassMomentumWeights::scaleAllReverseWeightsAt(Vec3i cellJ, Real factor)
    {
        for (Vec3i cellI : reverseWeights[index(cellJ)])
        {
            scale(cellI, cellJ, factor);
        }
    }

    void MassMomentumWeights::calculateBeta(Grid<Real> &beta)
    {
        Grid<Real> newBeta(beta.getParent());

        FOR_IJK(beta)
        {
            Vec3i cellI = Vec3i(i, j, k);
            for (const auto &[cellJ, value] : weights[index(cellI)])
            {
                newBeta(cellI) += value;
            }
        }

        beta.swap(newBeta);
    }

    void MassMomentumWeights::calculateGamma(Grid<Real> &gamma)
    {
        Grid<Real> newGamma(gamma.getParent());

        FOR_IJK(gamma)
        {
            Vec3i cellI = Vec3i(i, j, k);
            for (const auto &[cellJ, value] : weights[index(cellI)])
            {
                newGamma(cellJ) += value;
            }
        }

        gamma.swap(newGamma);
    }

    void MassMomentumWeights::calculateIntermediateResult(Grid<Real> &dest, Grid<Real> &src, Grid<Real> &min, Grid<Real> &max)
    {
        FOR_IJK(dest)
        {
            Vec3i cellI = Vec3i(i, j, k);

            if (weights[index(cellI)].empty())
            {
                continue;
            }

            for (const auto &[cellJ, w] : weights[index(cellI)])
            {
                dest(cellJ) += src(cellI) * w;
                min(cellJ) = std::min(min(cellJ), src(cellI));
                max(cellJ) = std::max(max(cellJ), src(cellI));
            }
        }
    }

#define EPSILON 1e-9

    void MassMomentumWeights::distributeLostMass(Grid<Real> &grid, Grid<Real> &lostMass, Grid<Real> &min, Grid<Real> &max, Real &subractedMass)
    {
        for (int _ = 0; _ < 4; _++)
        {
            FOR_IJK(grid)
            {
                Vec3i cellJ = Vec3i(i, j, k);
                std::vector<Vec3i> &reverseWeightsCellJ = reverseWeights[index(cellJ)];

                Real weight = 1.0 / (Real)reverseWeightsCellJ.size();
                Real massToMove = lostMass(cellJ);

                if (std::abs(massToMove) < EPSILON)
                {
                    continue;
                }

                for (auto it = reverseWeightsCellJ.begin(); it != reverseWeightsCellJ.end();)
                {
                    Vec3i cellI = (*it);

                    bool deleted = false;

                    // find the correct cell (bc it is only a vector)
                    for (auto &[_cellJ, w] : weights[index(cellI)])
                    {
                        if (_cellJ == cellJ) // found correct weight from cellI -> cellJ
                        {
                            if (min(cellI) <= grid(cellI) - massToMove * weight && grid(cellI) - massToMove * weight <= max(cellI))
                            {
                                grid(cellI) -= massToMove * weight;
                                lostMass(cellJ) -= massToMove * weight;
                                break;
                            }

                            Real cellIBefore = grid(cellI);
                            grid(cellI) = Manta::clamp(grid(cellI) - massToMove * weight, min(cellI), max(cellI));
                            Real massMoved = cellIBefore - grid(cellI);
                            lostMass(cellJ) -= massMoved;

                            // if (std::abs(massMoved - massToMove * weight) > 1e-6) // can't fit more mass
                            //{
                            reverseWeightsCellJ.erase(it);
                            deleted = true;
                            //}
                            break;
                        }
                    }

                    if (!deleted)
                    {
                        it++;
                    }
                }
            }
        }
    }
}
