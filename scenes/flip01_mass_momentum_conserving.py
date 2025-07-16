from manta import *
from data_collection import *
import json

params = {}
with open("../scenes/test_cases/params_high_clf.json") as f:
	params = json.load(f)

data_collector = Data_collectior(title="test_multiple_images", params=params, export_data=False, export_images=False,
								 trackable_grid_names=["phi_fluid", "_", "onesInFluid", "_", "curl", "pressure"], tracked_grids_indeces=[0, 2, 4, 5])
data_collector.init()

LEVEL = 0
doOpen = False
doConserving = False
doParticleLevelSetThomas = False
layout = 1 # 0=dam break, 1=falling drop

# solver params
dim = params["dimension"]
gs = vec3(params["resolutionX"], params["resolutionY"],params["resolutionZ"])
if (dim==2):
	gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)

# Adaptive time stepping
s.cfl         = params["maxCFL"]          				
s.frameLength = params["dt"]     
s.timestep    = s.frameLength 
s.timestepMin = 0.001 
s.timestepMax = 20000

# prepare grids and particles
flags_n          = s.create(FlagGrid)
flags_n_plus_one = s.create(FlagGrid)

innen1außen0       = s.create(RealGrid)
innen1außen0_gamma = s.create(RealGrid)
curl 			   = s.create(RealGrid)
pressure           = s.create(RealGrid)

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
reinitializeLevelset( phi=phi_fluid, flags=flags_n )

fluidbox.applyToGrid(grid=flags_n_plus_one, value=FlagFluid, respectFlags=flags_n_plus_one)
if layout == 1:
	fb2.applyToGrid(grid=flags_n_plus_one, value=FlagFluid, respectFlags=flags_n_plus_one)

if (GUI):
	gui = Gui()
	gui.show()
	gui.setCamPos(0, 0, -1.3)
	
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	#gui.nextVec3Grid()

	if layout == 1:
		pass

	if doConserving and doParticleLevelSetThomas:
		gui.nextParts()
		gui.nextParts()

	gui.pause()

firstFrame = True
#main loop
while (s.timeTotal < params["max_time"]):
	if firstFrame:
		fillWithOnes( grid=vel_gamma )
		fillLevelsetWithOnes( grid=innen1außen0, flags=flags_n, phi=phi_fluid, level=LEVEL )
		fillWithOnes( grid=innen1außen0_gamma )
		firstFrame = False

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
		advectSemiLagrange( flags=flags_n, vel=vel_extrapolated, grid=innen1außen0, order=2 )
		
		addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))
		setWallBcs(flags=flags_n, vel=vel)
		solvePressure(flags=flags_n, vel=vel, pressure=pressure)
		setWallBcs(flags=flags_n, vel=vel)
		extrapolateMACSimple( flags=flags_n, vel=vel )
		extrapolateMACSimple( flags=flags_n, vel=vel )
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags_n, parts=pp, partVel=pVel, flipRatio=0.97 )

	else:
		if not doParticleLevelSetThomas:
			extrapolateMACSimple( flags=flags_n, vel=vel_extrapolated, distance=5, intoObs=True )

			pp.advectInGrid(flags=flags_n, vel=vel_extrapolated, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagEmpty)
			eulerStep(parts=pp, vel=pVel, ptype=pT, exclude=FlagFluid)

			pp.projectOutOfBnd(flags=flags_n, bnd=bWidth + 1, plane='xXyYzZ', ptype=pT)
			if layout == 0:
				pushOutofObs(parts=pp, flags=flags_n, phiObs=phiObs, thresh=0.5*0.5, ptype=pT)
	

			markFluidCells(parts=pp, flags=flags_n_plus_one, ptype=pT)
			setPartType(parts=pp, ptype=pT, mark=FlagFluid, stype=FlagEmpty, flags=flags_n_plus_one, cflag=FlagFluid)
			markIsolatedFluidCell(flags=flags_n_plus_one, mark=FlagEmpty)
			setPartType(parts=pp, ptype=pT, mark=FlagEmpty, stype=FlagFluid, flags=flags_n_plus_one, cflag=FlagEmpty)

			gridParticleIndex(parts=pp, flags=flags_n_plus_one, indexSys=gIdxSys, index=gIdx)
			unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=flags_n_plus_one, index=gIdx, phi=phi_fluid, radiusFactor=1.0)
			extrapolateLsSimple(phi=phi_fluid, distance=5, inside=True)
			extrapolateLsSimple(phi=phi_fluid, distance=5, inside=False)

			flags_n_plus_one.copyFrom(flags_n_plus_one)
		
		else:
			vel_extrapolated.copyFrom(vel)
			extrapolateMACSimple( flags=flags_n, vel=vel_extrapolated, distance=3, intoObs=True )

			advectParticlesForward( particles=level_set_particles, vel=vel_extrapolated, flags=flags_n)
			#level_set_particles.advectInGrid(flags=flags_n, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

			level_set_particles.projectOutOfBnd(flags=flags_n, bnd=bWidth + 1, plane='xXyYzZ', ptype=particle_pT)
			if layout == 0:
				pushOutofObs(parts=level_set_particles, flags=flags_n, phiObs=phiObs, thresh=0.5*0.5, ptype=particle_pT)
			
			simpleSLAdvection( flags=flags_n, vel=vel_extrapolated, grid=phi_fluid )
			#advectSemiLagrange( flags=flags_n, vel=vel_extrapolated, grid=phi_fluid, order=2 )
			
			correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )
			reinitializeLevelset( phi=phi_fluid, flags=flags_n )
			correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )

			if (s.frame % 10 == 0):
				reseedParticles(phi=phi_fluid, flags=flags_n, particles=level_set_particles )

			reinitializeRadii( particles=level_set_particles, radii=particle_radii, phi=phi_fluid )

			setFlagsFromParticleLevelset( phi=phi_fluid, flags=flags_n_plus_one, level=LEVEL )

		# Advect grid quantities using the custom conserving advection scheme.
		massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=innen1außen0, gammaCumulative=innen1außen0_gamma, phi=phi_fluid)
		massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=vel_extrapolated, gammaCumulative=vel_gamma, phi=phi_fluid)
		
		vel.copyFrom(vel_extrapolated)
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
		
	calculateCurl(vel=vel, curl=curl)

	data_collector.step(solver=s, tracked_grids=[[innen1außen0, "test_real_grid"], [curl, "curl"]], flags=flags_n, vel=vel, gui=gui)
	s.step()

data_collector.finish()
data_collector.printStats()