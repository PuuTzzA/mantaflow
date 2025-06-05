#
# Simulation of a flame with smoke (and with adaptive time-stepping)
#
from manta import *

# solver params
dim = 2
res = 30
gs = vec3(res, 1.5 * res, res)
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# buoyancy parameters
smokeDensity = -0.001
smokeTempDiff = 0.1

# set time step range
s.frameLength = 1.2   # length of one frame (in "world time")
s.timestepMin = 0.2   # time step range
s.timestepMax = 2.0
s.cfl         = 3   # maximal velocity per cell
s.timestep    = (s.timestepMax+s.timestepMin)*0.5
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
react = s.create(RealGrid)
fuel = s.create(RealGrid)
heat = s.create(RealGrid)
flame = s.create(RealGrid)
pressure = s.create(RealGrid)
innen0außen1 = s.create(RealGrid)

density_gamma = s.create(RealGrid)
react_gamma = s.create(RealGrid)
fuel_gamma = s.create(RealGrid)
heat_gamma = s.create(RealGrid)
innen0außen1_gamma = s.create(RealGrid)
vel_gamma = s.create(MACGrid)

doOpen = True

# how many frames to calculate 
frames = 100

# noise field
noise = s.create(NoiseField, loadFromFile=True)
noise.posScale = vec3(45)
noise.clamp = True
noise.clampNeg = 1
noise.clampPos = 1
noise.valScale = 1
noise.valOffset = 0.75
noise.timeAnim = 0.2

gravity = vec3(0,-0.0981,0)

# vorticity global is applied to all cells, vorticity flames only to ones with fuel
vortGlobal = 0.1
vortFlames = 0.5

# initialize domain with boundary
bWidth=1
flags.initDomain( boundaryWidth=bWidth )
flags.fillGrid()

obsPos = vec3(0.5, 0.63, 0)
obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
obsSize = 0.11

obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
phiObs = obs.computeLevelset()

setObstacleFlags(flags=flags, phiObs=phiObs, boundaryWidth=4) 

flags.fillGrid()

obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

if doOpen:
	setOpenBound( flags, bWidth,'yY',FlagOutflow|FlagEmpty )

if (GUI):
	gui = Gui()
	gui.show(True)
	#gui.nextRealGrid()
	#gui.nextRealGrid()
	#gui.nextRealGrid()
	#gui.nextRealGrid()
	#gui.nextRealGrid()
	#gui.nextRealGrid()
	##gui.pause()

# source: cube in center of domain (x, y), standing on bottom of the domain
boxSize = vec3(res/8, 0.05*res, res/8)
boxCenter = gs*vec3(0.5, 0.15, 0.5)
sourceBox = s.create( Box, center=boxCenter, size=boxSize )

test_start = 0
test_min = 1000000
test_max = -1000000

# main loop
while s.frame < frames:
	maxvel = vel.getMax()
	s.adaptTimestep( maxvel )
	mantaMsg('\nFrame %i, time-step size %f' % (s.frame, s.timestep))

	if s.frame < 1:
		densityInflow( flags=flags, density=innen0außen1, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		fillWithOnes( grid=density_gamma )
		fillWithOnes( grid=heat_gamma )
		fillWithOnes( grid=fuel_gamma )
		fillWithOnes( grid=react_gamma )
		fillWithOnes( grid=innen0außen1_gamma )
		fillWithOnes( grid=vel_gamma )
	
	if s.timeTotal<200:
		densityInflow( flags=flags, density=density, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=heat, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=fuel, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=react, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )

	processBurn( fuel=fuel, density=density, react=react, heat=heat )

	if False:
		advectSemiLagrange( flags=flags, vel=vel, grid=density, order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=heat,   order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=fuel,   order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=react, order=2 )

		#advectSemiLagrange( flags=flags, vel=vel, grid=innen0außen1, order=2 )
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=innen0außen1, gammaCumulative=innen0außen1_gamma)

		advectSemiLagrange( flags=flags, vel=vel, grid=vel,   order=2 )
		#massMomentumConservingAdvect( flags=flags, vel=vel, grid=vel, gammaCumulative=gamma_cumulative)
	else:
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=density, gammaCumulative=density_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=heat, gammaCumulative=heat_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=fuel, gammaCumulative=fuel_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=react, gammaCumulative=react_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=innen0außen1, gammaCumulative=innen0außen1_gamma)

		massMomentumConservingAdvect( flags=flags, vel=vel, grid=vel, gammaCumulative=vel_gamma)
		#advectSemiLagrange( flags=flags, vel=vel, grid=vel,   order=2 )

	if doOpen:
		resetOutflow( flags=flags, real=density )

	# Apply global and fuel-based flame vorticity
	flame.copyFrom(fuel)
	flame.multConst(vortFlames) # temporarily misuse flame grid
	vorticityConfinement( vel=vel, flags=flags, strength=vortGlobal, strengthCell=flame)

	addBuoyancy( flags=flags, density=density, vel=vel, gravity=(gravity*smokeDensity ) )
	addBuoyancy( flags=flags, density=heat,    vel=vel, gravity=(gravity*smokeTempDiff) )

	setWallBcs( flags=flags, vel=vel )
	solvePressure( flags=flags, vel=vel, pressure=pressure )

	updateFlame( react=react, flame=flame )

	#timings.display()
	stats = calculateMass(grid=innen0außen1).split(",")
	mantaMsg(f"Total \"mass\" inside grid: {stats[0]}, min: {stats[1]}, max: {stats[2]}")

	if s.frame == 0:
		test_start = stats[0]
	test_min = min(test_min, float(stats[0]))
	test_max = max(test_max, float(stats[0]))

	s.step()

mantaMsg(f"Summary of Total Mass:\nStart: {test_start}, Min: {test_min}, Max: {test_max}")
