#
# Simple example scene for a 2D simulation
# Simulation of a buoyant smoke density plume with open boundaries at top & bottom
#
from manta import *
from data_collection import * 
import json
import os

print(os.getcwd())

params = {}
with open("../scenes/test_cases/gas/params_gas_low_cfl.json") as f:
	params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionX"]
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# scene params
doOpen = False
doObstacle = True
doConserving = True
exportData = True
exportImages = False
exportVideos = False
title = "Gas_Low_CFL_" + ("Conserving" if doConserving else "Traditional")
title = "____gas_with_gas"

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
pressure = s.create(RealGrid)
innen0außen1 = s.create(RealGrid)
curl = s.create(RealGrid)

vel_gamma = s.create(MACGrid)
density_gamma = s.create(RealGrid)
innen0außen1_gamma = s.create(RealGrid)

bWidth=1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

if doOpen:
	setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

if doObstacle:
	obsPos = vec3(0.5, 0.63, 0)
	obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
	obsSize = 0.11

	obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
	phiObs = obs.computeLevelset()

	setObstacleFlags(flags=flags, phiObs=phiObs, boundaryWidth=4) 

	flags.fillGrid()
	obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

if (GUI):
	gui = Gui()
	gui.show( True ) 
	#gui.pause()

source = s.create(Cylinder, center=gs*vec3(0.5,0.1,0.5), radius=res*0.14, z=gs*vec3(0, 0.02, 0))

#Data Colleciton and Export
data_collector = Data_collectior(title=title, params=params, export_data=exportData, export_images=exportImages, export_videos=exportVideos,
						trackable_grid_names=[["density", density], [], ["fixed_volume", innen0außen1], ["curl", curl], [], [], []], 
						tracked_grids_indeces=[0, 2, 3])

data_collector.init()

firstFrame = True
#main loop
while (s.timeTotal < params["max_time"]):
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)

	print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}")
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))


	if firstFrame:
		source.applyToGrid(grid=innen0außen1, value=1)
		fillWithOnes( grid=density_gamma )
		fillWithOnes( grid=innen0außen1_gamma )
		fillWithOnes( grid=vel_gamma )
		firstFrame = False

	if s.timeTotal<300:
		source.applyToGrid(grid=density, value=1)

	# optionally, dissolve smoke
	#dissolveSmoke(flags=flags, density=density, speed=4)

	if not doConserving:
		advectSemiLagrange(flags=flags, vel=vel, grid=density,      order=1) 
		advectSemiLagrange(flags=flags, vel=vel, grid=innen0außen1, order=1)
		advectSemiLagrange(flags=flags, vel=vel, grid=vel,          order=1)
	else:
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=density, gammaCumulative=density_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=innen0außen1, gammaCumulative=innen0außen1_gamma)

		massMomentumConservingAdvect( flags=flags, vel=vel, grid=vel, gammaCumulative=vel_gamma)

	if doOpen:
		resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-3,0), flags=flags)

	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	#timings.display()    
	calculateCurl(vel=vel, curl=curl)
	data_collector.step(solver=s, flags=flags, vel=vel, gui=gui)

	s.step()

data_collector.finish()
data_collector.printStats()