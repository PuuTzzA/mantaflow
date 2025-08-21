import sys
from manta import *
from data_collection import * 
import json
import os

params = {}
param_path = "../scenes/test_cases/test_tests/mass_momentum_conserving_highres_3d_test.json"
EXPORTS_BASE_DIR = "../exportsIgnore/highres_3d_yay/"

OBSTACLE_MESH_PATH = "../resources/highres_scene_obstacle_triangulated.obj"

if len(sys.argv) > 1:
    param_path = sys.argv[1]
    #EXPORTS_BASE_DIR = "../exports/3d_final/simple_plume_3d_high"
    EXPORTS_BASE_DIR = "../exports/7_highres_3d"

with open(param_path) as f:
    params = json.load(f)

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"], params["resolutionZ"])
res = params["resolutionY"]
if dim==2:
    gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# scene params
doOpen = params["doOpen"]
doObstacle = params["doObstacle"]
doConserving = params["doConserving"]
exportData = params["exportData"]
exportImages = params["exportImages"]
exportVideos = params["exportVideos"]
exportVDBs = params["exportVDBs"]

if exportVDBs and (exportImages or exportVideos):
    raise Exception("Cannot export both VDBs and Images") 

interpolationMethod = params["interpolationMethod"] # only important for tracingMethod == "RK4" or doConserving ==  true
tracingMethod = params["tracingMethod"] # only important for notCoserving (EE1: Explicit Euler Order 1, EE2: Explicit Euler with MAC cormack, RK$, Runge Kutta 4)
tracingFunction = params["tracingFunction"] # only for RK4 or doConserving, 0 = NORMAL, 1 = LOCAL_CFL
redistributeClamped = True

title = params["title"]

# set time step range
s.frameLength = params["dt"]     
s.cfl         = params["maxCFL"]
s.timestep    = s.frameLength
s.timestepMin = 0.1 
s.timestepMax = s.frameLength

timings = Timings()

# create grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)
curl = s.create(RealGrid)
velocity_magnitude = s.create(RealGrid)
phiObs = s.create(LevelsetGrid)

vel_gamma = s.create(MACGrid)
density_gamma = s.create(RealGrid)

if doObstacle:
    mesh = s.create(Mesh)
    mesh.load(OBSTACLE_MESH_PATH)
    #mesh.load("../resources/simpletorus.obj")
    mesh.scale( vec3(res) )
    mesh.offset( gs* (Vec3(0.5, 0.5, 0.5)) ) # center + slight offset

    mesh.computeLevelset(phiObs, 1)

#prepare grids
bWidth=0
flags.initDomain(boundaryWidth=bWidth) 

if doObstacle:
    setObstacleFlags(flags=flags, phiObs=phiObs) #, fractions=fractions)

if doOpen or True:
    setOpenBound(flags, bWidth,'xXYzZ',FlagOutflow|FlagEmpty) 

flags.fillGrid()

source = s.create(Cylinder, center=gs*vec3(0.5,0.125,0.5), radius=res*0.15, z=gs*vec3(0, 0.051, 0))
#source = s.create(Cylinder, center=gs*vec3(0.5,0.12,0.5), radius=res*0.15, z=gs*vec3(0, 0.04, 0))

fillWithOnes( grid=density_gamma )
fillWithOnes( grid=vel_gamma )

# guid and data collection
gui = None
if (GUI) and not exportVDBs:
    gui = Gui()
    gui.show( True ) 
    gui.pause()

data_collector = Data_collectior(title=title ,base_dir=EXPORTS_BASE_DIR, params=params, export_data=exportData, 
                                 export_images=exportImages, export_videos=exportVideos, export_vdbs=exportVDBs, 
                                 trackable_grid_names=[["density", density], [], ["curl", curl], ["vel_magnitude", velocity_magnitude], [], [], []], 
                                 tracked_grids_indeces=[0, 3], image_grids_indeces=[0], graph_grids=[["vel_magnitude", "max"]])

data_collector.init()

num_steps = 0
firstFrame = True
#main loop
while s.timeTotal < params["max_time"] and num_steps < 120:
    num_steps += 1
    
    computeVelocityMagnitude(dest=velocity_magnitude, vel=vel)
    maxvel = getMaxVal(grid=velocity_magnitude, flags=flags) # flags param does nothign for now

    if doConserving == False and tracingMethod == "EE1":
        maxvel = vel.getMax()

    if firstFrame:
        maxvel = 20     
        firstFrame = False

    s.adaptTimestep(maxvel)

    print(f"cfl number?: {maxvel * s.timestep}, timestep: {s.timestep}, maxvel: {maxvel}")
    mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

    if s.timeTotal<3000:
        source.applyToGrid(grid=density, value=1)

    # optionally, dissolve smoke
    # dissolveSmoke(flags=flags, density=density, speed=4)

    #Advection
    if not doConserving:
        if tracingMethod.startswith("EE"):
            order = int(tracingMethod[-1])  # Extracts 1 or 2 from "EE1" or "EE2"
            advectSemiLagrange(flags=flags, vel=vel, grid=density,      order=order) # ziemlich scheiße, hauptsachlich da es explicit Euler verwendet, nicht RK4 wie simpleSLAdvect 
            advectSemiLagrange(flags=flags, vel=vel, grid=vel,          order=order) # ziemlich scheiße, hauptsachlich da es explicit Euler verwendet, nicht RK4 wie simpleSLAdvect

            setOutflowToZero(flags=flags, grid=density)
            setOutflowToZero(flags=flags, grid=vel)

        elif tracingMethod == "RK4":
            simpleSLAdvect(flags=flags, vel=vel, grid=density,     interpolationType=interpolationMethod, tracingMethod=tracingFunction) # 0 = Trilinear, 1 = Cubic, 2= Polynomial Interpolation, 3 = monotonue cubib (hermite)
            simpleSLAdvect(flags=flags, vel=vel, grid=vel,         interpolationType=interpolationMethod, tracingMethod=tracingFunction) # 0 = Trilinear, 1 = Cubic, 2= Polynomial Interpolation, 3 = monotonue cubib (hermite)

    else:
        #simpleSLAdvect(flags=flags, vel=vel, grid=density,      interpolationType=1) # 0 = Trilinear, 1 = Catmull Rom
        massMomentumConservingAdvect( flags=flags, vel=vel, grid=density, gammaCumulative=density_gamma,          interpolationType=interpolationMethod, tracingMethod=tracingFunction, redistributeClamped=redistributeClamped)
        massMomentumConservingAdvect( flags=flags, vel=vel, grid=vel, gammaCumulative=vel_gamma,                  interpolationType=interpolationMethod, tracingMethod=tracingFunction, redistributeClamped=redistributeClamped)

    if doOpen:
        resetOutflow(flags=flags,real=density) 

    # Pressure Solve
    setWallBcs(flags=flags, vel=vel)    
    addBuoyancy(density=density, vel=vel, gravity=vec3(0,-4e-3,0), flags=flags)

    solvePressure(flags=flags, vel=vel, pressure=pressure)
    
    # Output
    calculateCurl(vel=vel, curl=curl, flags=flags)
    computeVelocityMagnitude(dest=velocity_magnitude, vel=vel)
    maxVel = getMaxVal(grid=velocity_magnitude, flags=flags) # flags param does nothign for now

    if doConserving == False and tracingMethod == "EE1":
        maxvel = vel.getMax()
        storeVelocityMagnitude(dest=velocity_magnitude, vel=vel)

    data_collector.step(solver=s, flags=flags, maxVel=maxVel, gui=gui, objects=[density, velocity_magnitude])

    timings.display()    
    s.step()

data_collector.finish()
data_collector.printStats()