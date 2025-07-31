import sys
from manta import *
from data_collection import * 
import json
import os

params = {}

param_path = "../scenes/test_cases/simple_plume_tests/params_simple_plume_3d_test.json"
#param_path = "../scenes/test_cases/gas/params_gas_low_cfl.json"

if len(sys.argv) > 1:
    param_path = sys.argv[1]

with open(param_path) as f:
	params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionX"]
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# scene params
if False:
	doOpen = False
	doObstacle = False
	doConserving = False
	exportData = True
	exportImages = False
	exportVideos = False
	exportVDBs = False if GUI else True
	exportVDBs = True
else:
	doOpen = params["doOpen"]
	doObstacle = params["doObstacle"]
	doConserving = params["doConserving"]
	exportData = params["exportData"]
	exportImages = params["exportImages"]
	exportVideos = params["exportVideos"]
	exportVDBs = params["exportVDBs"]

if exportVDBs and (exportImages or exportVideos):
	raise Exception("Cannot export both VDBs and Images") 

title = "simple_plume_cfl_20_" + ("conserving" if doConserving else "advectSemiLagrange")
title = "testtestestesets_volume"
title = params["title"]

# set time step range
s.cfl         = params["maxCFL"]
s.timestep    = params["dt"]
s.frameLength = params["dt"]     
s.timestepMin = 0.001 
s.timestepMax = 20000

timings = Timings()

# create grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)
innen0außen1 = s.create(RealGrid)
curl = s.create(RealGrid)

vel_gamma = s.create(MACGrid)
density_gamma = s.create(RealGrid)
innen0außen1_gamma = s.create(RealGrid)

#prepare grids
bWidth=2
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

if doOpen:
	setOpenBound(flags, bWidth,'xXYzZ',FlagOutflow|FlagEmpty) 

if doObstacle:
	obsPos = vec3(0.5, 0.55, 0.5)
	obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
	obsSize = 0.11

	obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
	phiObs = obs.computeLevelset()

	setObstacleFlags(flags=flags, phiObs=phiObs, boundaryWidth=4) 

	flags.fillGrid()
	obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

source = s.create(Cylinder, center=gs*vec3(0.5,0.12,0.5), radius=res*0.14, z=gs*vec3(0, 0.04, 0))

source.applyToGrid(grid=innen0außen1, value=1)
fillWithOnes( grid=density_gamma )
fillWithOnes( grid=innen0außen1_gamma )
fillWithOnes( grid=vel_gamma )

# guid and data collection
gui = None
if (GUI) and not exportVDBs:
	gui = Gui()
	gui.show( True ) 
	gui.pause()

data_collector = Data_collectior(title=title ,base_dir=f"../exports/simple_plume_3d/{title}/", params=params, export_data=exportData, 
								 export_images=exportImages, export_videos=exportVideos, export_vdbs=exportVDBs, 
								 trackable_grid_names=[["density", density], [], ["fixed_volume", innen0außen1], ["curl", curl], [], []], tracked_grids_indeces=[0, 2])

data_collector.init()

firstFrame = True
#main loop
while (s.timeTotal < params["max_time"] and data_collector.current_frame < 1000):
	
	#maxvel = vel.getMax()
	maxvel = compMaxInFluid(flags, vel)

	if firstFrame:
		maxvel = 15     
		firstFrame = False

	s.adaptTimestep(maxvel)

	print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}, maxvel: {maxvel}")
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))



	if s.timeTotal<3000:
		source.applyToGrid(grid=density, value=1)

	# optionally, dissolve smoke
	#dissolveSmoke(flags=flags, density=density, speed=4)

	if not doConserving:
		advectSemiLagrange(flags=flags, vel=vel, grid=density,      order=1) # ziemlich scheiße, hauptsachlich da es explicit Euler verwendet, nicht RK4 wie simpleSLAdvect 
		advectSemiLagrange(flags=flags, vel=vel, grid=innen0außen1, order=1) # ziemlich scheiße, hauptsachlich da es explicit Euler verwendet, nicht RK4 wie simpleSLAdvect
		advectSemiLagrange(flags=flags, vel=vel, grid=vel,          order=1) # ziemlich scheiße, hauptsachlich da es explicit Euler verwendet, nicht RK4 wie simpleSLAdvect

		#simpleSLAdvect(flags=flags, vel=vel, grid=density,           interpolationType=0) # 0 = Trilinear, 1 = Catmull Rom (doesn't work well at high CFL because of negative weights)
		#simpleSLAdvect(flags=flags, vel=vel, grid=innen0außen1,      interpolationType=0) # 0 = Trilinear, 1 = Catmull Rom (doesn't work well at high CFL because of negative weights)
		#simpleSLAdvect(flags=flags, vel=vel, grid=vel,               interpolationType=0) # 0 = Trilinear, 1 = Catmull Rom (doesn't work well at high CFL because of negative weights)

	else:
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=density, gammaCumulative=density_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=innen0außen1, gammaCumulative=innen0außen1_gamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=vel, gammaCumulative=vel_gamma)

	if doOpen:
		resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-3,0), flags=flags)

	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	calculateCurl(vel=vel, curl=curl, flags=flags)

	data_collector.step(solver=s, flags=flags, vel=vel, gui=gui, objects=[density])

	timings.display()    
	s.step()

data_collector.finish()
data_collector.printStats()