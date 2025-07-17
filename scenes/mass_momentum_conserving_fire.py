#
# Simulation of a flame with smoke (and with adaptive time-stepping)
#
from manta import *
from data_collection import * 
import json

params = {}
with open("../scenes/test_cases/params_gas_low_cfl.json") as f:
	params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionX"]
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# buoyancy parameters
smokeDensity = -0.001
smokeTempDiff = 0.1

doOpen = False
doObstacle = True
doConserving = False
exportData = False
exportImages = False

# set time step range
s.cfl         = params["maxCFL"]
s.frameLength = params["dt"]     
s.timestep    = s.frameLength 
s.timestepMin = 0.001 
s.timestepMax = 20000

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
curl = s.create(RealGrid)

density_gamma = s.create(RealGrid)
react_gamma = s.create(RealGrid)
fuel_gamma = s.create(RealGrid)
heat_gamma = s.create(RealGrid)
innen0außen1_gamma = s.create(RealGrid)
vel_gamma = s.create(MACGrid)

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

if doObstacle:
	obsPos = vec3(0.5, 0.63, 0)
	obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
	obsSize = 0.11

	obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
	phiObs = obs.computeLevelset()

	setObstacleFlags(flags=flags, phiObs=phiObs, boundaryWidth=4) 

	flags.fillGrid()

	obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

if doOpen:
	setOpenBound( flags, bWidth,'Y',FlagOutflow|FlagEmpty )

if (GUI):
	gui = Gui()
	gui.show(True)
	#gui.pause()

# source: cube in center of domain (x, y), standing on bottom of the domain
boxSize = vec3(res/8, 0.05*res, res/8)
boxCenter = gs*vec3(0.5, 0.15, 0.5)
sourceBox = s.create( Box, center=boxCenter, size=boxSize )

# Data Collection and Export
data_collector = Data_collectior(title="Gas_Low_CFL", params=params, export_data=exportData, export_images=exportImages,
						trackable_grid_names=[["density", density], [], [], [], [], [], ["fixed_volume", innen0außen1], ["curl", curl], [], [], [], [], [], []], 
						tracked_grids_indeces=[0, 6, 7])

data_collector.init()

firstFrame = True
# main loop
while (s.timeTotal < params["max_time"]):
	
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)

	print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}")
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	if firstFrame:
		densityInflow( flags=flags, density=innen0außen1, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		fillWithOnes( grid=density_gamma )
		fillWithOnes( grid=heat_gamma )
		fillWithOnes( grid=fuel_gamma )
		fillWithOnes( grid=react_gamma )
		fillWithOnes( grid=innen0außen1_gamma )
		fillWithOnes( grid=vel_gamma )
		firstFrame = False
	
	if s.timeTotal<200:
		densityInflow( flags=flags, density=density, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=heat, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=fuel, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )
		densityInflow( flags=flags, density=react, noise=noise, shape=sourceBox, scale=1, sigma=0.5 )

	processBurn( fuel=fuel, density=density, react=react, heat=heat )

	if not doConserving:
		advectSemiLagrange( flags=flags, vel=vel, grid=density, order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=heat,   order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=fuel,   order=2 )
		advectSemiLagrange( flags=flags, vel=vel, grid=react, order=2 )

		advectSemiLagrange( flags=flags, vel=vel, grid=innen0außen1, order=2 )
		#massMomentumConservingAdvect( flags=flags, vel=vel, grid=innen0außen1, gammaCumulative=innen0außen1_gamma)

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

	calculateCurl(vel=vel, curl=curl)
	data_collector.step(solver=s, flags=flags, vel=vel, gui=gui)

	s.step()

data_collector.finish()
data_collector.printStats()