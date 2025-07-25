from manta import *
from data_collection import * 
import json

params = {}
with open("../scenes/test_cases/fixed_vel_field/params_fixed_velocity_deformation_flow.json") as f:
	params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionX"]
if dim==2:
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# scenario selection, choose one of the four
#scenario = "diagonal_motion"
#scenario = "zalesak_rotation"
#scenario = "shear_flow"
scenario = "deformation_field"

# scene params
doConserving = True
exportData = True
exportImages = True
exportVideos = True
title = scenario + "_" + ("Conserving" if doConserving else "Traditional")

# set time step range
s.cfl         = params["maxCFL"]
s.timestep    = params["dt"] 
s.frameLength = 10000     
s.timestepMin = 0.001 
s.timestepMax = 1000

timings = Timings()

# create grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
testPhi = s.create(LevelsetGrid)
testField = s.create(RealGrid)
curl = s.create(RealGrid)

testPhiGamma = s.create(RealGrid)
testFieldGamma = s.create(RealGrid)

# prepare grids
bWidth=2
flags.initDomain( boundaryWidth=bWidth )
flags.fillGrid()

if scenario == "diagonal_motion":
	source = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.20)

	source.applyToGrid(grid=testField, value=1)
	testPhi.copyFrom(source.computeLevelset())

if scenario == "zalesak_rotation":
	circle = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.3)
	box = s.create(Box, p0=gs*vec3(0.47,0,0), p1=gs*vec3(.53, .5, 1.))

	testPhi.copyFrom(circle.computeLevelset())

	testPhi.multConst(-1)
	testPhi.join(box.computeLevelset())
	testPhi.multConst(-1)

	fillLevelsetWithOnes( grid=testField, flags=flags, phi=testPhi, level=0 )

if scenario == "shear_flow":
	source = s.create(Sphere, center=gs*vec3(0.5,0.3,0.5), radius=res*0.2)

	source.applyToGrid(grid=testField, value=1)
	testPhi.copyFrom(source.computeLevelset())

if scenario == "deformation_field":
	source = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.2)

	source.applyToGrid(grid=testField, value=1)
	testPhi.copyFrom(source.computeLevelset())


fillWithOnes( grid=testPhiGamma )
fillWithOnes( grid=testFieldGamma )
setVelocityField(vel=vel, flags=flags, functionName=scenario)

#Data Colleciton and Export
if GUI:
	gui = Gui()
	gui.show( True ) 

	data_collector = Data_collectior(title=title, params=params, export_data=exportData, export_images=exportImages, export_videos=exportVideos,
						trackable_grid_names=[["testPhi", testPhi], ["testField", testField], ["curl", curl], [], []], 
						tracked_grids_indeces=[0, 1], fixed_volume="testField")

	data_collector.init()

	gui.pause()

#main loop
while (s.timeTotal < params["max_time"]):
	maxvel = vel.getMax()
	s.adaptTimestep(maxvel)

	print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}, maxVel: {maxvel}")
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	if not doConserving:
		advectSemiLagrange(flags=flags, vel=vel, grid=testPhi, order=1) 
		advectSemiLagrange(flags=flags, vel=vel, grid=testField, order=1) 
	else:
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=testPhi, gammaCumulative=testPhiGamma)
		massMomentumConservingAdvect( flags=flags, vel=vel, grid=testField, gammaCumulative=testFieldGamma)

	#timings.display()    
	calculateCurl(vel=vel, curl=curl, flags=flags)

	if GUI:
		data_collector.step(solver=s, flags=flags, vel=vel, gui=gui)

	# optionally save some of the simulation objects to an OpenVDB file (requires compilation with -DOPENVDB=1)
	if True and s.frame % 1 == 0 and not GUI:
		# note: when saving pdata fields, they must be accompanied by and listed before their parent pp
		objects = [flags, pressure, density, vel, curl]
		save( name='../analysis/test_vdb_export/test_vdb_export_%04d.vdb' % s.frame, objects=objects )

	s.step()

if GUI:
	data_collector.finish()
	data_collector.printStats()