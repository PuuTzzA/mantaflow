from manta import *
from data_collection import * 
import json
import sys

params = {}
param_path = "../scenes/test_cases/test_tests/fixed_velocity_test.json"
EXPORTS_BASE_DIR = "../exportsIgnore/test/"

if len(sys.argv) > 1:
    param_path = sys.argv[1]
    EXPORTS_BASE_DIR = "../exports/1_fixed_vel_shear_flow_low_cfl__/"

with open(param_path) as f:
    params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionX"]
if dim==2:
    gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# scenario selection, choose one of the four: "diagonal_motion","zalesak_rotation","shear_flow","deformation_field"
scenario = params["scenario"]

# scene params
title = params["title"]

doConserving = params["doConserving"]
exportData = params["exportData"]
exportImages = params["exportImages"]
exportVideos = params["exportVideos"]
exportVDBs = params["exportVDBs"]

if exportVDBs and (exportImages or exportVideos):
    raise Exception("Cannot export both VDBs and Images") 

interpolationMethod = params["interpolationMethod"] # only important for tracingMethod == "RK4"
tracingMethod = params["tracingMethod"] # only important for notCoserving (EE1: Explicit Euler Order 1, EE2: Explicit Euler with MAC cormack, RK$, Runge Kutta 4)
redistributeClamped = params["redistributeClamped"]

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
testPhi = s.create(LevelsetGrid)
testField = s.create(RealGrid)
velocity_magnitude = s.create(RealGrid)

testPhiGamma = s.create(RealGrid)
testFieldGamma = s.create(RealGrid)

# prepare grids
bWidth=1
flags.initDomain( boundaryWidth=bWidth )
flags.fillGrid()

#setOpenBound(flags, bWidth,'xXyY',FlagOutflow|FlagEmpty) 

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

    gui.nextRealGrid()
    gui.windowSize(1000, 1000)
    gui.setCamPos(0, 0, -1.2)
    gui.update()
    gui.screenshot(str(Path(EXPORTS_BASE_DIR).expanduser().resolve() / f"{params["scenario"]}_first_frame.png"))
 
    #gui.pause()

data_collector = Data_collectior(title=title, base_dir=EXPORTS_BASE_DIR, params=params, 
                                 export_data=exportData, export_images=exportImages, export_videos=exportVideos, export_vdbs=exportVDBs,
                                 trackable_grid_names=[["testField", testField], ["vel_magnitude", velocity_magnitude], [], [], ["testPhi", testPhi]], 
                                 tracked_grids_indeces=[0, 1], image_grids_indeces=[0], graph_grids=[["vel_magnitude", "max"], ["testField", "sum"]])

data_collector.init()

maxvel = vel.getMax()
s.adaptTimestep(maxvel)

#main loop
while (s.timeTotal < params["max_time"]):

    print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}, maxVel: {maxvel}")
    mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

    if not doConserving:
        if tracingMethod.startswith("EE"):
            order = int(tracingMethod[-1])  # Extracts 1 or 2 from "EE1" or "EE2"
            advectSemiLagrange(flags=flags, vel=vel, grid=testPhi,   order=order) 
            advectSemiLagrange(flags=flags, vel=vel, grid=testField, order=order)

        elif tracingMethod == "RK4": #Runge Kutta 4 
            simpleSLAdvect(flags=flags, vel=vel, grid=testPhi,  interpolationType=interpolationMethod) # 0 = linear, 1 = cubic, 2 = polynomial interpolation, 3 = monotone hermite
            simpleSLAdvect(flags=flags, vel=vel, grid=testField,interpolationType=interpolationMethod) # 0 = linear, 1 = cubic, 2 = polynomial interpolation, 3 = monotone hermite 

    else:
        print(data_collector.current_frame)
        massMomentumConservingAdvect( flags=flags, vel=vel, grid=testPhi, gammaCumulative=testPhiGamma    ,interpolationType=interpolationMethod, redistributeClamped=redistributeClamped)
        massMomentumConservingAdvect( flags=flags, vel=vel, grid=testField, gammaCumulative=testFieldGamma,interpolationType=interpolationMethod, redistributeClamped=redistributeClamped)

    #timings.display()    
    computeVelocityMagnitude(dest=velocity_magnitude, vel=vel)
    maxVel = getMaxVal(grid=velocity_magnitude, flags=flags) # flags param does nothign for now

    data_collector.step(solver=s, flags=flags, maxVel=maxVel, gui=gui, objects=[])

    s.step()

data_collector.finish()
data_collector.printStats()