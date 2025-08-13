from manta import *
from data_collection import *
import json
import sys

params = {}
param_path = "../scenes/test_cases/test_tests/mmc_liquid_test.json"
EXPORTS_BASE_DIR = "../exports/test/"

if len(sys.argv) > 1:
    param_path = sys.argv[1]
    EXPORTS_BASE_DIR = "../exports/mmc_liquid/"

with open(param_path) as f:
    params = json.load(f)


LEVEL = 0
doOpen = params["doOpen"]
doConserving = params["doConserving"]
doParticleLevelSet = params["doParticleLevelSet"]
interpolationMethod = params["interpolationMethod"]

exportData = params["exportData"]
exportImages = params["exportImages"]
exportVideos = params["exportVideos"]
exportVDBs = params["exportVDBs"]

layout = 0 if params["layout"] == "dam" else 1 # 0=dam break, 1=falling drop

title = params["title"]

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"],params["resolutionZ"])
if (dim==2):
    gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# Adaptive time stepping
s.cfl         = params["maxCFL"]
s.timestep    = params["dt"]
s.frameLength = params["dt"]     
s.timestepMin = 0.001 
s.timestepMax = 20000

# prepare grids and particles
flags_n          = s.create(FlagGrid)
flags_n_plus_one = s.create(FlagGrid)

innen1außen0       = s.create(RealGrid)
innen1außen0_gamma = s.create(RealGrid)
curl 			   = s.create(RealGrid)
pressure           = s.create(RealGrid)
phi_gamma          = s.create(RealGrid)

vel              = s.create(MACGrid)
vel_gamma        = s.create(MACGrid)
vel_extrapolated = s.create(MACGrid)
velOld           = s.create(MACGrid)
tmpVec3          = s.create(VecGrid)

phi_fluid = s.create(LevelsetGrid)

gIdx = s.create(IntGrid)

# FLIP particles
pp       = s.create(BasicParticleSystem) 
pVel     = pp.create(PdataVec3) 
pT       = pp.create(PdataInt)
gIdxSys  = s.create(ParticleIndexSystem)

# Particle Level Set particles
level_set_particles = s.create(BasicParticleSystem)
particle_radii      = level_set_particles.create(PdataReal)
particle_pT         = level_set_particles.create(PdataInt)

# scene setup
bWidth=2
flags_n.initDomain( boundaryWidth=bWidth )
flags_n.fillGrid()
flags_n_plus_one.initDomain( boundaryWidth=bWidth)
flags_n_plus_one.fillGrid()
flags_n_plus_one.initDomain( boundaryWidth=bWidth )

if layout == 0:
    obs = Box( parent=s, p0=gs*vec3(0.5,0,0), p1=gs*vec3(0.6,0.1,1))
    phiObs = obs.computeLevelset()

    setObstacleFlags(flags=flags_n, phiObs=phiObs) 
    setObstacleFlags(flags=flags_n_plus_one, phiObs=phiObs)

flags_n.fillGrid()
flags_n_plus_one.fillGrid()

phiInit = None
fluidbox = None
fb2 = None

if layout == 0:
    fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.25, 0.5, 1)) # breaking dam
    phiInit = fluidbox.computeLevelset()
elif layout == 1:
    fluidbox = Sphere( parent=s , center=gs*Vec3(0.5,0.5,0.5), radius=params["resolutionX"]*0.125) # centered falling sphere
    fluidbox2 = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1, .1, 1)) # centered falling block
    phiInit = fluidbox.computeLevelset()
    phiInit.join(fluidbox2.computeLevelset())
    fb2 = fluidbox2

flags_n.updateFromLevelset(phiInit)
flags_n_plus_one.updateFromLevelset(phiInit)

# FLIP particles
begin = pp.pySize()
sampleShapeWithParticles(shape=fluidbox, flags=flags_n_plus_one, parts=pp, discretization=5, randomness=0, notiming=True)
if layout == 1:
    sampleShapeWithParticles(shape=fb2, flags=flags_n_plus_one, parts=pp, discretization=5, randomness=0, notiming=True)

end = pp.pySize()
pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)

# level set method 
phi_fluid.copyFrom(phiInit)
phi_fluid.reinitMarching( flags=flags_n )

sampleLevelsetBorderWithParticles( phi=phi_fluid, flags=flags_n, particles=level_set_particles, radii=particle_radii)

fluidbox.applyToGrid(grid=flags_n_plus_one, value=FlagFluid, respectFlags=flags_n_plus_one)
if layout == 1:
    fb2.applyToGrid(grid=flags_n_plus_one, value=FlagFluid, respectFlags=flags_n_plus_one)

fillWithOnes( grid=vel_gamma )
fillLevelsetWithOnes( grid=innen1außen0, flags=flags_n, phi=phi_fluid, level=LEVEL )
fillWithOnes( grid=innen1außen0_gamma )
fillWithOnes( grid=phi_gamma )

if (GUI):
    gui = Gui()
    gui.show()
    gui.setCamPos(0, 0, -1.3)
    
    gui.nextRealGrid()
    gui.nextRealGrid()
    gui.nextRealGrid()
    gui.nextRealGrid()
    gui.nextRealGrid()
    
    gui.nextVec3Grid()
    gui.nextVec3Grid()

    if layout == 1:
        pass

    if doConserving and doParticleLevelSet:
        gui.nextParts()
        gui.nextParts()
    
    if doConserving and not doParticleLevelSet:
        gui.nextParts()
        
    gui.pause()

# Data collection and exportation
data_collector = None
if layout == 0:
    data_collector = Data_collectior(title=title ,base_dir=EXPORTS_BASE_DIR, params=params, export_data=exportData, 
                                    export_images=exportImages, export_videos=exportVideos, export_vdbs=exportVDBs, 
                                    trackable_grid_names=[["phi_fluid", phi_fluid], [], [], ["fixed_volume", innen1außen0], [], ["curl", curl], [], []], 
                                    tracked_grids_indeces=[0, 3], image_grids_indeces=[0], graph_grids=[["fixed_volume", "max"]])

if layout == 1:
    data_collector = Data_collectior(title=title ,base_dir=EXPORTS_BASE_DIR, params=params, export_data=exportData, 
                                    export_images=exportImages, export_videos=exportVideos, export_vdbs=exportVDBs, 
                                    trackable_grid_names=[["phi_fluid", phi_fluid], [], ["fixed_volume", innen1außen0], [], ["curl", curl], [], []], 
                                    tracked_grids_indeces=[0, 2], image_grids_indeces=[0], graph_grids=[["fixed_volume", "max"]])

data_collector.init()

firstFrame = True
#main loop
while (s.timeTotal < params["max_time"]):
    computeVelocityMagnitude(dest=curl, vel=vel)
    maxvel = getMaxVal(grid=curl, flags=flags_n) # flags param does nothign for now

    if firstFrame:
        maxvel = 15     
        firstFrame = False

    s.adaptTimestep(maxvel)

    print(f"cfl number?: {vel.getMaxAbs() * s.timestep}, timestep: {s.timestep}")
    mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

    if not doConserving:
        ## set (just for visualization)
        markPhiFromFlagGrid(phi_fluid, flags_n)
        reinitializeLevelset( phi=phi_fluid, flags=flags_n )

        # --- Standard FLIP simulation loop for comparison ---
        pp.advectInGrid(flags=flags_n, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )
        mapPartsToMAC(vel=vel, flags=flags_n, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 )
        extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )
        markFluidCells( parts=pp, flags=flags_n )
        # Standard advection for the test grid
        vel_extrapolated.copyFrom(vel)
        extrapolateMACSimple( flags=flags_n, vel=vel_extrapolated, distance=3, intoObs=True )
        advectSemiLagrange( flags=flags_n, vel=vel_extrapolated, grid=innen1außen0, order=1 )
        
        addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))
        setWallBcs(flags=flags_n, vel=vel)
        solvePressure(flags=flags_n, vel=vel, pressure=pressure)
        setWallBcs(flags=flags_n, vel=vel)
        extrapolateMACSimple( flags=flags_n, vel=vel )
        flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags_n, parts=pp, partVel=pVel, flipRatio=0.97 )

    else:
        vel_extrapolated.copyFrom(vel)
        extrapolateMACSimple( flags=flags_n, vel=vel_extrapolated, distance=10, intoObs=True )

        if not doParticleLevelSet:
            simpleSLAdvect(flags=flags_n, vel=vel_extrapolated, grid=phi_fluid, interpolationType=2, all=True) # 0 = Trilinear, 1 = cubic convolutional, 2 = polynomial, 3 = monotone hermite
            #advectSemiLagrange( flags=flags_n, vel=vel_extrapolated, grid=phi_fluid, order=2 )
            #reinitializeLevelset( phi=phi_fluid, flags=flags_n )

            setFlagsFromParticleLevelset( phi=phi_fluid, flags=flags_n_plus_one, level=LEVEL )

        else:
            """             advectParticlesForward( particles=level_set_particles, vel=vel_extrapolated, flags=flags_n)
                        #level_set_particles.advectInGrid(flags=flags_n, vel=vel_extrapolated, integrationMode=IntRK4, deleteInObstacle=False )

                        level_set_particles.projectOutOfBnd(flags=flags_n, bnd=bWidth + 1, plane='xXyYzZ', ptype=particle_pT)
                        if layout == 0:
                            pushOutofObs(parts=level_set_particles, flags=flags_n, phiObs=phiObs, thresh=0.5*0.5, ptype=particle_pT)
                        
                        #simpleSLAdvect(flags=flags_n, vel=vel_extrapolated, grid=phi_fluid, interpolationType=2, all=True) # 0 = Trilinear, 1 = cubic convolutional, 2 = polynomial, 3 = monotone hermite
                        advectSemiLagrange( flags=flags_n, vel=vel_extrapolated, grid=phi_fluid, order=2 )
                        #massMomentumConservingAdvect( flags=flags_n, vel=vel, grid=phi_fluid, gammaCumulative=phi_gamma)

                        correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )
                        reinitializeLevelset( phi=phi_fluid, flags=flags_n )
                        correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )

                        reinitializeRadii( particles=level_set_particles, radii=particle_radii, phi=phi_fluid )
            """            
            advectParticleLevelSet( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, vel=vel_extrapolated, flags=flags_n )

            if (data_collector.current_frame % 20 == 0):
                pass
                reseedParticles(phi=phi_fluid, flags=flags_n, particles=level_set_particles, radii=particle_radii )

            setFlagsFromParticleLevelset( phi=phi_fluid, flags=flags_n_plus_one, level=LEVEL )

        # Advect grid quantities using the custom conserving advection scheme.
        massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=innen1außen0, gammaCumulative=innen1außen0_gamma, phi=phi_fluid, interpolationType=interpolationMethod)
        massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=vel, gammaCumulative=vel_gamma, phi=phi_fluid                  , interpolationType=interpolationMethod)
        
        flags_n.copyFrom(flags_n_plus_one)
        
        # Apply grid-based forces (gravity)
        addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))

        # Solve for pressure to make the velocity field divergence-free
        setWallBcs(flags=flags_n, vel=vel)    
        solvePressure(flags=flags_n, vel=vel, pressure=pressure)
        #solvePressure(flags=flags_n, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=5e-4, phi=phi_fluid )

        setWallBcs(flags=flags_n, vel=vel)

    if doOpen:
        resetOutflow( flags=flags_n)
        

    # Data Collection
    computeVelocityMagnitude(dest=curl, vel=vel)
    maxVel = getMaxVal(grid=curl, flags=flags_n) # flags param does nothign for now
    calculateCurl(vel=vel, curl=curl, flags=flags_n)
    data_collector.step(solver=s, flags=flags_n, maxVel=maxVel, gui=gui, objects=[flags_n])

    s.step()

data_collector.finish()
data_collector.printStats()