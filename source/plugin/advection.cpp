/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2015 Tobias Pfaff, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugins for advection
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <array>

using namespace std;

namespace Manta
{
	template <typename T>
	using Sparse2DMap = std::unordered_map<IndexInt, std::unordered_map<IndexInt, T>>;
	using Reverse2dMap = std::unordered_map<IndexInt, std::unordered_set<IndexInt>>;

	//! Semi-Lagrange interpolation kernel
	KERNEL(bnd = 1)
	template <class T>
	void SemiLagrange(const FlagGrid &flags, const MACGrid &vel, Grid<T> &dst, const Grid<T> &src, Real dt, bool isLevelset, int orderSpace, int orderTrace)
	{
		if (orderTrace == 1)
		{
			// traceback position
			Vec3 pos = Vec3(i + 0.5f, j + 0.5f, k + 0.5f) - vel.getCentered(i, j, k) * dt;
			dst(i, j, k) = src.getInterpolatedHi(pos, orderSpace);
		}
		else if (orderTrace == 2)
		{
			// backtracing using explicit midpoint
			Vec3 p0 = Vec3(i + 0.5f, j + 0.5f, k + 0.5f);
			Vec3 p1 = p0 - vel.getCentered(i, j, k) * dt * 0.5;
			Vec3 p2 = p0 - vel.getInterpolated(p1) * dt;
			dst(i, j, k) = src.getInterpolatedHi(p2, orderSpace);
		}
		else
		{
			assertMsg(false, "Unknown backtracing order " << orderTrace);
		}
	}

	//! Semi-Lagrange interpolation kernel for MAC grids
	KERNEL(bnd = 1)
	void SemiLagrangeMAC(const FlagGrid &flags, const MACGrid &vel, MACGrid &dst, const MACGrid &src, Real dt, int orderSpace, int orderTrace)
	{
		if (orderTrace == 1)
		{
			// get currect velocity at MAC position
			// no need to shift xpos etc. as lookup field is also shifted
			Vec3 xpos = Vec3(i + 0.5f, j + 0.5f, k + 0.5f) - vel.getAtMACX(i, j, k) * dt;
			Real vx = src.getInterpolatedComponentHi<0>(xpos, orderSpace);
			Vec3 ypos = Vec3(i + 0.5f, j + 0.5f, k + 0.5f) - vel.getAtMACY(i, j, k) * dt;
			Real vy = src.getInterpolatedComponentHi<1>(ypos, orderSpace);
			Vec3 zpos = Vec3(i + 0.5f, j + 0.5f, k + 0.5f) - vel.getAtMACZ(i, j, k) * dt;
			Real vz = src.getInterpolatedComponentHi<2>(zpos, orderSpace);

			dst(i, j, k) = Vec3(vx, vy, vz);
		}
		else if (orderTrace == 2)
		{
			Vec3 p0 = Vec3(i + 0.5, j + 0.5, k + 0.5);
			Vec3 xp0 = Vec3(i, j + 0.5f, k + 0.5f);
			Vec3 xp1 = xp0 - src.getAtMACX(i, j, k) * dt * 0.5;
			Vec3 xp2 = p0 - src.getInterpolated(xp1) * dt;
			Real vx = src.getInterpolatedComponentHi<0>(xp2, orderSpace);
			Vec3 yp0 = Vec3(i + 0.5f, j, k + 0.5f);
			Vec3 yp1 = yp0 - src.getAtMACY(i, j, k) * dt * 0.5;
			Vec3 yp2 = p0 - src.getInterpolated(yp1) * dt;
			Real vy = src.getInterpolatedComponentHi<1>(yp2, orderSpace);
			Vec3 zp0 = Vec3(i + 0.5f, j + 0.5f, k);
			Vec3 zp1 = zp0 - src.getAtMACZ(i, j, k) * dt * 0.5;
			Vec3 zp2 = p0 - src.getInterpolated(zp1) * dt;
			Real vz = src.getInterpolatedComponentHi<2>(zp2, orderSpace);

			dst(i, j, k) = Vec3(vx, vy, vz);
		}
		else
		{
			assertMsg(false, "Unknown backtracing order " << orderTrace);
		}
	}

	//! Kernel: Correct based on forward and backward SL steps (for both centered & mac grids)
	KERNEL(idx)
	template <class T>
	void MacCormackCorrect(const FlagGrid &flags, Grid<T> &dst, const Grid<T> &old, const Grid<T> &fwd, const Grid<T> &bwd,
						   Real strength, bool isLevelSet, bool isMAC = false)
	{
		dst[idx] = fwd[idx];

		if (flags.isFluid(idx))
		{
			// only correct inside fluid region; note, strenth of correction can be modified here
			dst[idx] += strength * 0.5 * (old[idx] - bwd[idx]);
		}
	}

	//! Kernel: Correct based on forward and backward SL steps (for both centered & mac grids)
	KERNEL()
	template <class T>
	void MacCormackCorrectMAC(const FlagGrid &flags, Grid<T> &dst, const Grid<T> &old, const Grid<T> &fwd, const Grid<T> &bwd,
							  Real strength, bool isLevelSet, bool isMAC = false)
	{
		bool skip[3] = {false, false, false};

		if (!flags.isFluid(i, j, k))
			skip[0] = skip[1] = skip[2] = true;
		if (isMAC)
		{
			if ((i > 0) && (!flags.isFluid(i - 1, j, k)))
				skip[0] = true;
			if ((j > 0) && (!flags.isFluid(i, j - 1, k)))
				skip[1] = true;
			if ((k > 0) && (!flags.isFluid(i, j, k - 1)))
				skip[2] = true;
		}

		for (int c = 0; c < 3; ++c)
		{
			if (skip[c])
			{
				dst(i, j, k)[c] = fwd(i, j, k)[c];
			}
			else
			{
				// perform actual correction with given strength
				dst(i, j, k)[c] = fwd(i, j, k)[c] + strength * 0.5 * (old(i, j, k)[c] - bwd(i, j, k)[c]);
			}
		}
	}

	// Helper to collect min/max in a template
	template <class T>
	inline void getMinMax(T &minv, T &maxv, const T &val)
	{
		if (val < minv)
			minv = val;
		if (val > maxv)
			maxv = val;
	}
	template <>
	inline void getMinMax<Vec3>(Vec3 &minv, Vec3 &maxv, const Vec3 &val)
	{
		getMinMax(minv.x, maxv.x, val.x);
		getMinMax(minv.y, maxv.y, val.y);
		getMinMax(minv.z, maxv.z, val.z);
	}

	//! detect out of bounds value
	template <class T>
	inline bool cmpMinMax(T &minv, T &maxv, const T &val)
	{
		if (val < minv)
			return true;
		if (val > maxv)
			return true;
		return false;
	}
	template <>
	inline bool cmpMinMax<Vec3>(Vec3 &minv, Vec3 &maxv, const Vec3 &val)
	{
		return (cmpMinMax(minv.x, maxv.x, val.x) | cmpMinMax(minv.y, maxv.y, val.y) | cmpMinMax(minv.z, maxv.z, val.z));
	}

#define checkFlag(x, y, z) (flags((x), (y), (z)) & (FlagGrid::TypeFluid | FlagGrid::TypeEmpty))

	//! Helper function for clamping non-mac grids (those have specialized per component version below)
	//  Note - 2 clamp modes, a sharper one (default, clampMode 1, also uses backward step),
	//         and a softer version (clampMode 2) that is recommended in Andy's paper
	template <class T>
	inline T doClampComponent(const Vec3i &gridSize, const FlagGrid &flags, T dst, const Grid<T> &orig, const T fwd, const Vec3 &pos, const Vec3 &vel, const int clampMode)
	{
		T minv(std::numeric_limits<Real>::max()), maxv(-std::numeric_limits<Real>::max());
		bool haveFl = false;

		// forward (and optionally) backward
		Vec3i positions[2];
		int numPos = 1;
		positions[0] = toVec3i(pos - vel);
		if (clampMode == 1)
		{
			numPos = 2;
			positions[1] = toVec3i(pos + vel);
		}

		for (int l = 0; l < numPos; ++l)
		{
			Vec3i &currPos = positions[l];

			// clamp lookup to grid
			const int i0 = clamp(currPos.x, 0, gridSize.x - 1); // note! gridsize already has -1 from call
			const int j0 = clamp(currPos.y, 0, gridSize.y - 1);
			const int k0 = clamp(currPos.z, 0, (orig.is3D() ? (gridSize.z - 1) : 1));
			const int i1 = i0 + 1, j1 = j0 + 1, k1 = (orig.is3D() ? (k0 + 1) : k0);

			// find min/max around source pos
			if (checkFlag(i0, j0, k0))
			{
				getMinMax(minv, maxv, orig(i0, j0, k0));
				haveFl = true;
			}
			if (checkFlag(i1, j0, k0))
			{
				getMinMax(minv, maxv, orig(i1, j0, k0));
				haveFl = true;
			}
			if (checkFlag(i0, j1, k0))
			{
				getMinMax(minv, maxv, orig(i0, j1, k0));
				haveFl = true;
			}
			if (checkFlag(i1, j1, k0))
			{
				getMinMax(minv, maxv, orig(i1, j1, k0));
				haveFl = true;
			}

			if (orig.is3D())
			{
				if (checkFlag(i0, j0, k1))
				{
					getMinMax(minv, maxv, orig(i0, j0, k1));
					haveFl = true;
				}
				if (checkFlag(i1, j0, k1))
				{
					getMinMax(minv, maxv, orig(i1, j0, k1));
					haveFl = true;
				}
				if (checkFlag(i0, j1, k1))
				{
					getMinMax(minv, maxv, orig(i0, j1, k1));
					haveFl = true;
				}
				if (checkFlag(i1, j1, k1))
				{
					getMinMax(minv, maxv, orig(i1, j1, k1));
					haveFl = true;
				}
			}
		}

		if (!haveFl)
			return fwd;
		if (clampMode == 1)
		{
			dst = clamp(dst, minv, maxv); // hard clamp
		}
		else
		{
			if (cmpMinMax(minv, maxv, dst))
				dst = fwd; // recommended in paper, "softer"
		}
		return dst;
	}

	//! Helper function for clamping MAC grids, slight differences in flag checks
	//  similar to scalar version, just uses single component c of vec3 values
	//  for symmetry, reverts to first order near boundaries for clampMode 2
	template <int c>
	inline Real doClampComponentMAC(const FlagGrid &flags, const Vec3i &gridSize, Real dst, const MACGrid &orig, Real fwd, const Vec3 &pos, const Vec3 &vel, const int clampMode)
	{
		Real minv = std::numeric_limits<Real>::max(), maxv = -std::numeric_limits<Real>::max();
		// bool haveFl = false;

		// forward (and optionally) backward
		Vec3i positions[2];
		int numPos = 1;
		positions[0] = toVec3i(pos - vel);
		if (clampMode == 1)
		{
			numPos = 2;
			positions[1] = toVec3i(pos + vel);
		}

		Vec3i oPos = toVec3i(pos);
		Vec3i nbPos = oPos;
		nbPos[c] -= 1;
		if (clampMode == 2 && (!(checkFlag(oPos.x, oPos.y, oPos.z) && checkFlag(nbPos.x, nbPos.y, nbPos.z))))
			return fwd; // replaces haveFl check

		for (int l = 0; l < numPos; ++l)
		{
			Vec3i &currPos = positions[l];

			const int i0 = clamp(currPos.x, 0, gridSize.x - 1); // note! gridsize already has -1 from call
			const int j0 = clamp(currPos.y, 0, gridSize.y - 1); // but we need a clamp to -2 for the +1 offset below
			const int k0 = clamp(currPos.z, 0, (orig.is3D() ? (gridSize.z - 1) : 0));
			const int i1 = i0 + 1, j1 = j0 + 1, k1 = (orig.is3D() ? (k0 + 1) : k0);

			// find min/max around source pos
			getMinMax(minv, maxv, orig(i0, j0, k0)[c]);
			getMinMax(minv, maxv, orig(i1, j0, k0)[c]);
			getMinMax(minv, maxv, orig(i0, j1, k0)[c]);
			getMinMax(minv, maxv, orig(i1, j1, k0)[c]);

			if (orig.is3D())
			{
				getMinMax(minv, maxv, orig(i0, j0, k1)[c]);
				getMinMax(minv, maxv, orig(i1, j0, k1)[c]);
				getMinMax(minv, maxv, orig(i0, j1, k1)[c]);
				getMinMax(minv, maxv, orig(i1, j1, k1)[c]);
			}
		}

		if (clampMode == 1)
		{
			dst = clamp(dst, minv, maxv); // hard clamp
		}
		else
		{
			if (cmpMinMax(minv, maxv, dst))
				dst = fwd; // recommended in paper, "softer"
		}
		return dst;
	}

#undef checkFlag

	//! Kernel: Clamp obtained value to min/max in source area, and reset values that point out of grid or into boundaries
	//          (note - MAC grids are handled below)
	KERNEL(bnd = 1)
	template <class T>
	void MacCormackClamp(const FlagGrid &flags, const MACGrid &vel, Grid<T> &dst, const Grid<T> &orig, const Grid<T> &fwd, Real dt, const int clampMode)
	{
		T dval = dst(i, j, k);
		Vec3i gridUpper = flags.getSize() - 1;

		dval = doClampComponent<T>(gridUpper, flags, dval, orig, fwd(i, j, k), Vec3(i, j, k), vel.getCentered(i, j, k) * dt, clampMode);

		if (1 && clampMode == 1)
		{
			// lookup forward/backward , round to closest NB
			Vec3i posFwd = toVec3i(Vec3(i, j, k) + Vec3(0.5, 0.5, 0.5) - vel.getCentered(i, j, k) * dt);
			Vec3i posBwd = toVec3i(Vec3(i, j, k) + Vec3(0.5, 0.5, 0.5) + vel.getCentered(i, j, k) * dt);

			// test if lookups point out of grid or into obstacle (note doClampComponent already checks sides, below is needed for valid flags access)
			if (posFwd.x < 0 || posFwd.y < 0 || posFwd.z < 0 ||
				posBwd.x < 0 || posBwd.y < 0 || posBwd.z < 0 ||
				posFwd.x > gridUpper.x || posFwd.y > gridUpper.y || ((posFwd.z > gridUpper.z) && flags.is3D()) ||
				posBwd.x > gridUpper.x || posBwd.y > gridUpper.y || ((posBwd.z > gridUpper.z) && flags.is3D()) ||
				flags.isObstacle(posFwd) || flags.isObstacle(posBwd))
			{
				dval = fwd(i, j, k);
			}
		}
		// clampMode 2 handles flags in doClampComponent call

		dst(i, j, k) = dval;
	}

	//! Kernel: same as MacCormackClamp above, but specialized version for MAC grids
	KERNEL(bnd = 1)
	void MacCormackClampMAC(const FlagGrid &flags, const MACGrid &vel, MACGrid &dst, const MACGrid &orig, const MACGrid &fwd, Real dt, const int clampMode)
	{
		Vec3 pos(i, j, k);
		Vec3 dval = dst(i, j, k);
		Vec3 dfwd = fwd(i, j, k);
		Vec3i gridUpper = flags.getSize() - 1;

		dval.x = doClampComponentMAC<0>(flags, gridUpper, dval.x, orig, dfwd.x, pos, vel.getAtMACX(i, j, k) * dt, clampMode);
		dval.y = doClampComponentMAC<1>(flags, gridUpper, dval.y, orig, dfwd.y, pos, vel.getAtMACY(i, j, k) * dt, clampMode);
		if (flags.is3D())
			dval.z = doClampComponentMAC<2>(flags, gridUpper, dval.z, orig, dfwd.z, pos, vel.getAtMACZ(i, j, k) * dt, clampMode);

		// note - the MAC version currently does not check whether source points were inside an obstacle! (unlike centered version)
		// this would need to be done for each face separately to stay symmetric...

		dst(i, j, k) = dval;
	}

	//! template function for performing SL advection
	//! (Note boundary width only needed for specialization for MAC grids below)
	template <class GridType>
	void fnAdvectSemiLagrange(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &orig, int order, Real strength, int orderSpace, int clampMode, int orderTrace)
	{
		typedef typename GridType::BASETYPE T;

		Real dt = parent->getDt();
		bool levelset = orig.getType() & GridBase::TypeLevelset;

		// forward step
		GridType fwd(parent); // source
		SemiLagrange<T>(flags, vel, fwd, orig, dt, levelset, orderSpace, orderTrace);

		if (order == 1)
		{
			orig.swap(fwd);
		}
		else if (order == 2)
		{						  // MacCormack
			GridType bwd(parent); // destination
			GridType newGrid(parent);

			// bwd <- backwards step
			SemiLagrange<T>(flags, vel, bwd, fwd, -dt, levelset, orderSpace, orderTrace);

			// newGrid <- compute correction
			MacCormackCorrect<T>(flags, newGrid, orig, fwd, bwd, strength, levelset);

			// clamp values
			MacCormackClamp<T>(flags, vel, newGrid, orig, fwd, dt, clampMode);

			orig.swap(newGrid);
		}
	}

	// outflow functions

	//! calculate local propagation velocity for cell (i,j,k)
	Vec3 getBulkVel(const FlagGrid &flags, const MACGrid &vel, int i, int j, int k)
	{
		Vec3 avg = Vec3(0.);
		int count = 0;
		int size = 1; // stencil size
		int nmax = (flags.is3D() ? size : 0);
		// average the neighboring fluid / outflow cell's velocity
		for (int n = -nmax; n <= nmax; n++)
		{
			for (int m = -size; m <= size; m++)
			{
				for (int l = -size; l <= size; l++)
				{
					if (flags.isInBounds(Vec3i(i + l, j + m, k + n)) && (flags.isFluid(i + l, j + m, k + n) || flags.isOutflow(i + l, j + m, k + n)))
					{
						avg += vel(i + l, j + m, k + n);
						count++;
					}
				}
			}
		}
		return count > 0 ? avg / count : avg;
	}

	//! extrapolate normal velocity components into outflow cell
	KERNEL()
	void extrapolateVelConvectiveBC(const FlagGrid &flags, const MACGrid &vel, MACGrid &velDst, const MACGrid &velPrev, Real timeStep)
	{
		if (flags.isOutflow(i, j, k))
		{
			Vec3 bulkVel = getBulkVel(flags, vel, i, j, k);
			int dim = flags.is3D() ? 3 : 2;
			const Vec3i cur = Vec3i(i, j, k);
			Vec3i low, up, flLow, flUp;
			int cnt = 0;
			// iterate over each velocity component x, y, z
			for (int c = 0; c < dim; c++)
			{
				low = up = flLow = flUp = cur;
				Real factor = timeStep * max((Real)1.0, bulkVel[c]); // prevent the extrapolated velocity from exploding when bulk velocity below 1
				low[c] = flLow[c] = cur[c] - 1;
				up[c] = flUp[c] = cur[c] + 1;
				// iterate over bWidth to allow for extrapolation into more distant outflow cells; hard-coded extrapolation distance of two cells
				for (int d = 0; d < 2; d++)
				{
					bool extrapolateFromLower = flags.isInBounds(flLow) && flags.isFluid(flLow);
					bool extrapolateFromUpper = flags.isInBounds(flUp) && flags.isFluid(flUp);
					if (extrapolateFromLower || extrapolateFromUpper)
					{
						if (extrapolateFromLower)
						{
							velDst(i, j, k) += ((vel(i, j, k) - velPrev(i, j, k)) / factor) + vel(low);
							cnt++;
						}
						if (extrapolateFromUpper)
						{
							// check for cells equally far away from two fluid cells -> average value between both sides
							velDst(i, j, k) += ((vel(i, j, k) - velPrev(i, j, k)) / factor) + vel(up);
							cnt++;
						}
						break;
					}
					flLow[c]--;
					flUp[c]++;
				}
			}
			if (cnt > 0)
				velDst(i, j, k) /= cnt;
		}
	}

	//! copy extrapolated velocity components
	KERNEL()
	void copyChangedVels(const FlagGrid &flags, const MACGrid &velDst, MACGrid &vel)
	{
		if (flags.isOutflow(i, j, k))
			vel(i, j, k) = velDst(i, j, k);
	}

	//! extrapolate normal velocity components into open boundary cells (marked as outflow cells)
	void applyOutflowBC(const FlagGrid &flags, MACGrid &vel, const MACGrid &velPrev, double timeStep)
	{
		MACGrid velDst(vel.getParent()); // do not overwrite vel while it is read
		extrapolateVelConvectiveBC(flags, vel, velDst, velPrev, max(1.0, timeStep * 4));
		copyChangedVels(flags, velDst, vel);
	}

	// advection helpers

	//! prevent parts of the surface getting "stuck" in obstacle regions
	KERNEL()
	void knResetPhiInObs(const FlagGrid &flags, Grid<Real> &sdf)
	{
		if (flags.isObstacle(i, j, k) && (sdf(i, j, k) < 0.))
		{
			sdf(i, j, k) = 0.1;
		}
	}
	PYTHON()
	void resetPhiInObs(const FlagGrid &flags, Grid<Real> &sdf) { knResetPhiInObs(flags, sdf); }

	// advection main calls

	//! template function for performing SL advection: specialized version for MAC grids
	template <>
	void fnAdvectSemiLagrange<MACGrid>(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, MACGrid &orig, int order, Real strength, int orderSpace, int clampMode, int orderTrace)
	{
		Real dt = parent->getDt();

		// forward step
		MACGrid fwd(parent);
		SemiLagrangeMAC(flags, vel, fwd, orig, dt, orderSpace, orderTrace);

		if (orderSpace != 1)
		{
			debMsg("Warning higher order for MAC grids not yet implemented...", 1);
		}

		if (order == 1)
		{
			applyOutflowBC(flags, fwd, orig, dt);
			orig.swap(fwd);
		}
		else if (order == 2)
		{ // MacCormack
			MACGrid bwd(parent);
			MACGrid newGrid(parent);

			// bwd <- backwards step
			SemiLagrangeMAC(flags, vel, bwd, fwd, -dt, orderSpace, orderTrace);

			// newGrid <- compute correction
			MacCormackCorrectMAC<Vec3>(flags, newGrid, orig, fwd, bwd, strength, false, true);

			// clamp values
			MacCormackClampMAC(flags, vel, newGrid, orig, fwd, dt, clampMode);

			applyOutflowBC(flags, newGrid, orig, dt);
			orig.swap(newGrid);
		}
	}

	//! Perform semi-lagrangian advection of target Real- or Vec3 grid
	//! Open boundary handling needs information about width of border
	//! Clamping modes: 1 regular clamp leading to more overshoot and sharper results, 2 revert to 1st order slightly smoother less overshoot (enable when 1 gives artifacts)
	PYTHON()
	void advectSemiLagrange(const FlagGrid *flags, const MACGrid *vel, GridBase *grid,
							int order = 1, Real strength = 1.0, int orderSpace = 1, bool openBounds = false, int boundaryWidth = -1, int clampMode = 2, int orderTrace = 1)
	{
		assertMsg(order == 1 || order == 2, "AdvectSemiLagrange: Only order 1 (regular SL) and 2 (MacCormack) supported");
		if ((boundaryWidth != -1) || (openBounds))
		{
			debMsg("Warning: boundaryWidth and openBounds parameters in AdvectSemiLagrange plugin are deprecated (and have no more effect), please remove.", 0);
		}

		// determine type of grid
		if (grid->getType() & GridBase::TypeReal)
		{
			fnAdvectSemiLagrange<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else if (grid->getType() & GridBase::TypeMAC)
		{
			fnAdvectSemiLagrange<MACGrid>(flags->getParent(), *flags, *vel, *((MACGrid *)grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else if (grid->getType() & GridBase::TypeVec3)
		{
			fnAdvectSemiLagrange<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else
			errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");
	}

	// THOMAS mass and momentum conserving advection step

	//! Semi-Lagrange interpolation kernel
	KERNEL(bnd = 1)
	void fillHelper(Grid<Real> &dst)
	{
		dst(i, j, k) = 1;
	}

	KERNEL(bnd = 1)
	void fillHelperMAC(MACGrid &dst)
	{
		dst(i, j, k) = Vec3(1., 1., 1.);
	}

	PYTHON()
	void fillWithOnes(GridBase *grid)
	{
		if (grid->getType() & GridBase::TypeReal)
		{
			fillHelper(*((Grid<Real> *)grid)).run();
			std::cout << "filled Grid<Real> with ones" << std::endl;
		}
		else if (grid->getType() & GridBase::TypeMAC)
		{
			fillHelperMAC(*((MACGrid *)grid)).run();
			std::cout << "filled MACGrid with ones" << std::endl;
		}
	}

	Vec3 customTrace(Vec3 pos, const MACGrid &vel, Real dt, const FlagGrid &flags, Vec3i &gs)
	{
		if (flags.isObstacle(pos))
		{
			throw runtime_error("trace not starting from fluid cell");
		}

		Vec3 k1 = vel.getInterpolatedHi(pos, 2);
		Vec3 k2 = vel.getInterpolatedHi(pos + dt / 2. * k1, 2);
		Vec3 k3 = vel.getInterpolatedHi(pos + dt / 2. * k2, 2);
		Vec3 k4 = vel.getInterpolatedHi(pos + dt * k3, 2);

		Vec3 nextPos = pos + dt * 1. / 6. * (k1 + 2. * k2 + 2. * k3 + k4);

		nextPos.x = Manta::clamp(nextPos.x, static_cast<Real>(0.), static_cast<Real>(gs[0] - 1));
		nextPos.y = Manta::clamp(nextPos.y, static_cast<Real>(0.), static_cast<Real>(gs[1] - 1));
		nextPos.z = Manta::clamp(nextPos.z, static_cast<Real>(0.), static_cast<Real>(gs[2] - 1));

		if (!flags.isObstacle(nextPos))
		{
			return nextPos;
		}

		Vec3 segmentStart = pos;
		Vec3 segmentEnd = nextPos;
		Vec3 lastKnownFluidPos = pos;

		Vec3 direction = segmentEnd - segmentStart;
		Real totalDistance = norm(direction); // Assuming Vec3 has a norm() or length() method

		if (totalDistance < 1e-9)
		{
			if (!flags.isObstacle(segmentStart))
				return segmentStart;

			return pos;
		}

		int numSearchSteps = std::max(2, static_cast<int>(std::ceil(totalDistance / 0.5))); // e.g. step by 0.25 grid units

		for (int i = 1; i <= numSearchSteps; ++i)
		{
			Real t = static_cast<Real>(i) / static_cast<Real>(numSearchSteps);
			Vec3 currentTestPoint = segmentEnd - t * direction;

			if (!flags.isObstacle(currentTestPoint))
			{
				return currentTestPoint;
			}
		}
		return pos;
	}

	std::vector<std::tuple<Vec3i, Real>> getInterpolationstencilAndWeights(const FlagGrid &flags, Vec3 x, Vec3i gs, Vec3 &offset)
	{
		int i = std::floor(x[0] - offset[0]);
		int j = std::floor(x[1] - offset[1]);
		i = Manta::clamp(i, 0, gs[0] - 2);
		j = Manta::clamp(j, 0, gs[1] - 2);

		Real fx = x[0] - offset[0] - i;
		Real fy = x[1] - offset[1] - j;

		Real w00 = 0., w10 = 0., w01 = 0., w11 = 0.;

		if (!flags.isObstacle(i, j, 0))
		{
			w00 = (1 - fx) * (1 - fy);
		}
		if (!flags.isObstacle(i + 1, j, 0))
		{
			w10 = fx * (1 - fy);
		}
		if (!flags.isObstacle(i, j + 1, 0))
		{
			w01 = (1 - fx) * fy;
		}
		if (!flags.isObstacle(i + 1, j + 1, 0))
		{
			w11 = fx * fy;
		}

		Real tot = w00 + w01 + w10 + w11;

		if (tot < 1e-5)
		{
			return {};
		}

		w00 /= tot;
		w10 /= tot;
		w01 /= tot;
		w11 /= tot;

		std::vector<std::tuple<Vec3i, Real>> result{};
		result.reserve(4);

		const Real EPSILON = 1e-5;

		if (w00 > EPSILON)
		{
			result.push_back({Vec3i{i, j, 0}, w00});
		}
		if (w10 > EPSILON)
		{
			result.push_back({Vec3i{i + 1, j, 0}, w10});
		}
		if (w01 > EPSILON)
		{
			result.push_back({Vec3i{i, j + 1, 0}, w01});
		}
		if (w11 > EPSILON)
		{
			result.push_back({Vec3i{i + 1, j + 1, 0}, w11});
		}

		return result;
	}

	KERNEL()
	template <class T>
	void advectGammaCum(const MACGrid &vel, Grid<T> &grid, Grid<T> &newGrid, float dt, Vec3i gridSize, Vec3 &offset, const FlagGrid &flags, Vec3i &gs)
	{
		if (!flags.isFluid(i, j, k))
		{
			newGrid(i, j, k) = 1;
			return;
		}

		Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
		newPos = customTrace(newPos, vel, -dt, flags, gs);

		auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
		for (const auto &[n, w] : neighboursAndWeights)
		{
			newGrid(i, j, k) += w * grid(n);
		}
	}

	KERNEL()
	template <class T>
	void setNewGammaCum(Grid<Real> &gammaCum, std::vector<Real> gamma, Vec3i gridSize)
	{
		gammaCum(i, j, k) = gamma[i * gridSize[1] + j];
	}

	void recalculateBeta(std::vector<Real> &beta, const Sparse2DMap<Real> &weights)
	{
		std::fill(beta.begin(), beta.end(), 0);

		for (const auto &[cellI, innerMap] : weights)
		{
			for (const auto &[cellJ, value] : innerMap)
			{
				beta[cellI] += value;
			}
		}
	}

	void recalculateGamma(std::vector<Real> &gamma, const Sparse2DMap<Real> &weights)
	{
		std::fill(gamma.begin(), gamma.end(), 0);
		for (const auto &[cellI, innerMap] : weights)
		{
			for (const auto &[cellJ, value] : innerMap)
			{
				gamma[cellJ] += value;
			}
		}
	}

	template <class GridType>
	void fnMassMomentumConservingAdvect(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, GridType &grid, Grid<Real> &gammaCumulative, Vec3 offset)
	{
		typedef typename GridType::BASETYPE T;
		const Real EPSILON = 1e-5;

		// For testing of the "normal" advection step that is used in this function
		/* Grid<T> testGrid(parent);
		advectGammaCum<T>(vel, grid, testGrid, parent->getDt(), parent->getGridSize(), offset, flags);
		grid.swap(testGrid);
		return; */

		GridType newGrid(parent); // source

		// Advect the cummulative Gamma the same way as later the rest
		Real dt = parent->getDt();
		Vec3i gridSize = parent->getGridSize();
		Grid<Real> newGammaCum(parent);
		advectGammaCum<Real>(vel, gammaCumulative, newGammaCum, dt, gridSize, offset, flags, gridSize);
		gammaCumulative.swap(newGammaCum);

		// main advection part
		long unsigned numCells = gridSize[0] * gridSize[1] * gridSize[2];

		// weights[k][p] = weight from cell k to cell p (cell indeces k/p = i * gridSize[0] + j)
		Sparse2DMap<Real> weights;
		Reverse2dMap reverseWeights;
		std::vector<Real> beta(numCells, 0.);
		std::vector<Real> gamma(numCells, 0.);

		int bnd = 1;

		// Step 1: backwards step
		for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
		{
			for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
			{
				int k = 0;

				if (flags.isObstacle(i, j, k))
				{
					continue;
				}

				IndexInt cellJ = i * gridSize[1] + j;

				Vec3 newPos = Vec3(i + offset[0], j + offset[1], k + offset[2]);
				newPos = customTrace(newPos, vel, -dt, flags, gridSize);

				auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, newPos, gridSize, offset);
				for (const auto &[n, w] : neighboursAndWeights)
				{
					IndexInt cellI = n[0] * gridSize[1] + n[1];
					weights[cellI][cellJ] = w;
					reverseWeights[cellJ].insert(cellI);
					beta[cellI] += w;

					// newGrid(i, j, k) += interpolationWeight(newPos, n) * grid(n);
				}
			}
		}

		// Step 2: forward step for all beta < 1
		for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
		{
			for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
			{
				int k = 0;
				if (flags.isObstacle(i, j, k))
				{
					continue;
				}

				IndexInt cellI = i * gridSize[1] + j;

				if (beta[cellI] < 1 - EPSILON)
				{
					Vec3 posForward = Vec3(i + offset[0], j + offset[1], k + offset[2]);
					posForward = customTrace(posForward, vel, dt, flags, gridSize);

					Real amountToDistribute = 1. - beta[cellI];

					auto neighboursAndWeights = getInterpolationstencilAndWeights(flags, posForward, gridSize, offset);
					for (const auto &[n, w] : neighboursAndWeights)
					{
						IndexInt cellJ = n[0] * gridSize[1] + n[1];

						weights[cellI][cellJ] += w * amountToDistribute;
						reverseWeights[cellJ].insert(cellI);
					}
				}
			}
		}

		// Step 3: compute Gamma
		recalculateGamma(gamma, weights);

		// Step 4: Clamp gamma to the cumulative gamma
		for (IndexInt i = bnd; i < gridSize[0] - bnd; i++)
		{
			for (IndexInt j = bnd; j < gridSize[1] - bnd; j++)
			{
				int k = 0;

				if (flags.isObstacle(i, j, k))
				{
					continue;
				}

				IndexInt cellJ = i * gridSize[1] + j;

				if (gamma[cellJ] < EPSILON)
					continue; // avoid division by 0

				Real factor = gammaCumulative(i, j, k) / gamma[cellJ];

				for (IndexInt cellI : reverseWeights[cellJ])
				{
					weights[cellI][cellJ] *= factor;
				}
			}
		}

		// Step 5: Clamp beta to 1 for conservation
		recalculateBeta(beta, weights); // should be 1 for all beta, maybe no need to recalculate, just set to 1
		for (IndexInt cellI = 0; cellI < numCells; cellI++)
		{
			if (flags.isObstacle(cellI / gridSize[1], cellI % gridSize[1], 0))
			{
				continue;
			}

			if (beta[cellI] < EPSILON)
				continue; // avoid division by 0

			Real factor = 1 / beta[cellI];

			for (auto &[_, value] : weights[cellI])
			{
				value *= factor;
			}
		}

		// Step 6: calculate the an intermediate result
		for (const auto &[cellI, innerMap] : weights)
		{
			for (const auto &[cellJ, weight] : innerMap)
			{
				int k = 0;

				IndexInt cellI_i = cellI / gridSize[1];
				IndexInt cellI_j = cellI % gridSize[1];

				IndexInt cellJ_i = cellJ / gridSize[1];
				IndexInt cellJ_j = cellJ % gridSize[1];

				newGrid(cellJ_i, cellJ_j, k) += weight * grid(cellI_i, cellI_j, k);
			}
		}

		// Step 7: Diffuse gamma using Gaus seidel Sweep
		recalculateGamma(gamma, weights);

		Real testMaxFlux = 0.1;

		for (int _ = 0; _ < 5; _++)
		{
			// X-Dimension
			for (IndexInt y = bnd; y < gridSize[1] - bnd; y++)
			{
				for (IndexInt x = bnd; x < gridSize[0] - bnd - 1; x++)
				{
					int k = 0;
					IndexInt cellI = x * gridSize[1] + y;
					IndexInt cellI_1 = (x + 1) * gridSize[1] + y;

					if (!flags.isFluid(x, y, k) || !flags.isFluid(x + 1, y, k)) // 1 == Type Fluid
					{
						continue;
					}

					Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;
					fluxGamma = Manta::clamp(fluxGamma, 0.f, testMaxFlux);

					gamma[cellI] += fluxGamma;
					gamma[cellI_1] -= fluxGamma;

					T gammaToMove = newGrid(x + 1, y, k) * (fluxGamma / gamma[cellI_1]);
					newGrid(x, y, k) += gammaToMove;
					newGrid(x + 1, y, k) -= gammaToMove;

					newGrid(x, y, k) = Manta::clamp(newGrid(x, y, k), 0.f, newGrid(x, y, k));
					newGrid(x + 1, y, k) = Manta::clamp(newGrid(x + 1, y, k), 0.f, newGrid(x + 1, y, k));
				}
			}

			// Y-Dimension
			for (IndexInt x = bnd; x < gridSize[0] - bnd; x++)
			{
				for (IndexInt y = bnd; y < gridSize[1] - bnd - 1; y++)
				{
					int k = 0;
					IndexInt cellI = x * gridSize[1] + y;
					IndexInt cellI_1 = cellI + 1;

					if (!flags.isFluid(x, y, k) || !flags.isFluid(x, y + 1, k)) // 1 == Type Fluid
					{
						continue;
					}

					Real fluxGamma = (gamma[cellI_1] - gamma[cellI]) / 2;
					fluxGamma = Manta::clamp(fluxGamma, 0.f, testMaxFlux);

					gamma[cellI] += fluxGamma;
					gamma[cellI_1] -= fluxGamma;

					T gammaToMove = newGrid(x, y + 1, k) * (fluxGamma / gamma[cellI_1]);
					newGrid(x, y, k) += gammaToMove;
					newGrid(x, y + 1, k) -= gammaToMove;

					newGrid(x, y, k) = Manta::clamp(newGrid(x, y, k), 0.f, newGrid(x, y, k));
					newGrid(x, y + 1, k) = Manta::clamp(newGrid(x, y + 1, k), 0.f, newGrid(x, y + 1, k));
				}
			}
		}

		setNewGammaCum<Real>(gammaCumulative, gamma, gridSize);
		grid.swap(newGrid);
		// MassMomentum<T>(flags, vel, grid, newGrid, parent->getDt());
		// grid.swap(newGrid);
	}

	KERNEL()
	void MAC2Grids(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
	{
		Vec3 data = vel.getAtMACnoInterpolation(i, j, k);
		velX(i, j, k) = data.x;
		velY(i, j, k) = data.y;
		velZ(i, j, k) = data.z;
	}

	KERNEL()
	void Grids2MAC(MACGrid &vel, Grid<Real> &velX, Grid<Real> &velY, Grid<Real> &velZ)
	{
		vel.setAtMACnoInterpolation(i, j, k, Vec3(velX(i, j, k), velY(i, j, k), velZ(i, j, k)));
	}

	void fnMassMomentumConservingAdvectMAC(FluidSolver *parent, const FlagGrid &flags, const MACGrid &vel, MACGrid &grid, MACGrid &gammaCumulative)
	{
		const Real EPSILON = 1e-5;
		Real dt = parent->getDt();
		Vec3i gridSize = parent->getGridSize();

		Grid<Real> velX(parent);
		Grid<Real> velY(parent);
		Grid<Real> velZ(parent);

		Grid<Real> gammaX(parent);
		Grid<Real> gammaY(parent);
		Grid<Real> gammaZ(parent);

		MAC2Grids(grid, velX, velY, velZ);
		MAC2Grids(gammaCumulative, gammaX, gammaY, gammaZ);

		Vec3 offsetX = Vec3(0., 0.5, 0.5);
		Vec3 offsetY = Vec3(0.5, 0., 0.5);
		Vec3 offsetZ = Vec3(0.5, 0.5, 0.);

		fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velX, gammaX, offsetX);
		fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velY, gammaY, offsetY);
		fnMassMomentumConservingAdvect<Grid<Real>>(parent, flags, vel, velZ, gammaZ, offsetZ);

		/* 		Grid<Real> newVelX(parent);
				Grid<Real> newVelY(parent);
				Grid<Real> newVelZ(parent);

				advectGammaCum<Real>(vel, velX, newVelX, dt, gridSize, offsetX, flags);
				advectGammaCum<Real>(vel, velY, newVelY, dt, gridSize, offsetY, flags);
				advectGammaCum<Real>(vel, velZ, newVelZ, dt, gridSize, offsetZ, flags);

				Grids2MAC(grid, newVelX, newVelY, newVelZ);
				return;
		 */
		Grids2MAC(grid, velX, velY, velZ);
		Grids2MAC(gammaCumulative, gammaX, gammaY, gammaZ);
		return;
	}

	PYTHON()
	void massMomentumConservingAdvect(const FlagGrid *flags, const MACGrid *vel, GridBase *grid, GridBase *gammaCumulative)
	{
		if (grid->getType() & GridBase::TypeReal)
		{
			fnMassMomentumConservingAdvect<Grid<Real>>(flags->getParent(), *flags, *vel, *((Grid<Real> *)grid), *((Grid<Real> *)gammaCumulative), Vec3(0.5, 0.5, 0.5));
		}
		else if (grid->getType() & GridBase::TypeMAC)
		{
			fnMassMomentumConservingAdvectMAC(flags->getParent(), *flags, *vel, *((MACGrid *)grid), *((MACGrid *)gammaCumulative));
		}
		else if (grid->getType() & GridBase::TypeVec3)
		{
			// fnMassMomentumConservingAdvect<Grid<Vec3>>(flags->getParent(), *flags, *vel, *((Grid<Vec3> *)grid), *((Grid<Real> *)gammaCumulative));
		}
		else
			errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");
	}

	PYTHON()
	std::string calculateMass(const Grid<Real> *grid)
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

} // end namespace DDF