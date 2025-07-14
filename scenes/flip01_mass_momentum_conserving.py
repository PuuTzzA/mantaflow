# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

RESOLUTION = 100
TITLE = "FLIP simple high CFL"
FILENAME = f'../analysis/data/{TITLE.replace(" ", "_")}.json'
CFL = 1
MAX_TIME = 600
NUM_FRAMES_RENDERED = 6
EXPORT = False

LEVEL = 0
doOpen = False
doConserving = True

# solver params
dim = 2
particleNumber = 2
res = RESOLUTION
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
	particleNumber = 3 # use more particles in 2d
s = Solver(name='main', gridSize = gs, dim=dim)

# Adaptive time stepping
s.cfl         = CFL          # maximal velocity per cell and timestep, 3 is fairly strict
s.frameLength = 0.8                 # length of one frame (in "world time")
s.timestep    = s.frameLength 
s.timestepMin = 0.001 
s.timestepMax = 20000

# prepare grids and particles
flags_n          = s.create(FlagGrid)
flags_n_plus_one = s.create(FlagGrid)
test_real_grid   = s.create(RealGrid)
test_real_grid_gamma = s.create(RealGrid)
vel              = s.create(MACGrid)
vel_gamma        = s.create(MACGrid)
vel_extrapolated = s.create(MACGrid)
velOld           = s.create(MACGrid)
pressure         = s.create(RealGrid)
tmpVec3          = s.create(VecGrid)

pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 

# scene setup
bWidth=2
flags_n.initDomain( boundaryWidth=bWidth )
flags_n.fillGrid()
flags_n_plus_one.initDomain( boundaryWidth=bWidth)
flags_n_plus_one.fillGrid()

obsPos = vec3(0.5, 0.1, 0)
obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
obsSize = 0.11

obs = Box( parent=s, p0=gs*vec3(0.5,0,0), p1=gs*vec3(0.6,0.1,1))
#obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
phiObs = obs.computeLevelset()

setObstacleFlags(flags=flags_n, phiObs=phiObs) 
setObstacleFlags(flags=flags_n_plus_one, phiObs=phiObs)

flags_n.fillGrid()
flags_n_plus_one.fillGrid()

# enable one of the following
fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.25, 0.5, 1)) # breaking dam
#fluidbox = Box( parent=s, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags_n.updateFromLevelset(phiInit)
flags_n_plus_one.updateFromLevelset(phiInit )

# level set method 
phi_fluid = fluidbox.computeLevelset()
# phi_fluid.initFromFlags( flags=flags_n ) # This is redundant due to updateFromLevelset above
phi_fluid.reinitMarching( flags=flags_n )

level_set_particles = s.create(BasicParticleSystem)
particle_radii = level_set_particles.create(PdataReal)
sampleLevelsetBorderWithParticles( phi=phi_fluid, flags=flags_n, particles=level_set_particles, radii=particle_radii)
reinitializeLevelset( phi=phi_fluid, flags=flags_n )

# note, there's no resampling here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags_n, parts=pp, discretization=particleNumber, randomness=0.2 )

if (GUI):
	gui = Gui()
	gui.show()
	gui.setCamPos(0, 0, -1.3)
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextRealGrid()
	gui.nextVec3Grid()
	gui.nextVec3Grid()
	if doConserving:
		gui.nextParts()
	#gui.pause()

pp.clearFile(FILENAME)

test_start = 0
test_min = 1000000
test_max = -1000000

#main loop
for t in range(MAX_TIME):
	if s.frame < 1:
		fillWithOnes( grid=vel_gamma )
		fillLevelsetWithOnes( grid=test_real_grid, flags=flags_n, phi=phi_fluid, level=LEVEL )
		fillWithOnes( grid=test_real_grid_gamma )

	maxVel = vel.getMax()
	s.adaptTimestep( maxVel )

	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	if not doConserving:
		# --- Standard FLIP simulation loop for comparison ---
		pp.advectInGrid(flags=flags_n, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )
		mapPartsToMAC(vel=vel, flags=flags_n, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 )
		extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )
		markFluidCells( parts=pp, flags=flags_n )
		# Standard advection for the test grid
		advectSemiLagrange( flags=flags_n, vel=vel, grid=test_real_grid, order=2 )
		
		addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))
		setWallBcs(flags=flags_n, vel=vel)
		solvePressure(flags=flags_n, vel=vel, pressure=pressure)
		setWallBcs(flags=flags_n, vel=vel)
		extrapolateMACSimple( flags=flags_n, vel=vel )
		flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags_n, parts=pp, partVel=pVel, flipRatio=0.97 )

	else:
		# Level Set
		vel_extrapolated.copyFrom(vel)
		extrapolateMACSimple( flags=flags_n, vel=vel_extrapolated, distance=5, intoObs=True )
		#extrapolateVelFSM( phi=phi_fluid, flags=flags_n, vel=vel_extrapolated, steps=5 )

		advectParticlesForward( particles=level_set_particles, vel=vel_extrapolated, flags=flags_n)
		simpleSLAdvection( flags=flags_n, vel=vel_extrapolated, grid=phi_fluid )
		correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )
		reinitializeLevelset( phi=phi_fluid, flags=flags_n )
		correctErrorsWithParticles( phi=phi_fluid, particles=level_set_particles, radii=particle_radii, flags=flags_n )

		if (t % 10 == 0):
			reseedParticles(phi=phi_fluid, flags=flags_n, particles=level_set_particles )

		reinitializeRadii( particles=level_set_particles, radii=particle_radii, phi=phi_fluid )

		setFlagsFromParticleLevelset( phi=phi_fluid, flags=flags_n_plus_one, level=LEVEL )

		# Advect grid quantities using the custom conserving advection scheme.
		# This is the core of the comparison.
		massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=test_real_grid, gammaCumulative=test_real_grid_gamma, phi=phi_fluid)
		massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel_extrapolated, grid=vel, gammaCumulative=vel_gamma, phi=phi_fluid)

		# The fluid domain has now officially moved. Update flags_n for the pressure solve.
		flags_n.copyFrom(flags_n_plus_one)
		
		# Apply grid-based forces (gravity)
		addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))

		# Solve for pressure to make the velocity field divergence-free
		setWallBcs(flags=flags_n, vel=vel)    
		solvePressure(flags=flags_n, vel=vel, pressure=pressure)
		setWallBcs(flags=flags_n, vel=vel)

	if doOpen:
		resetOutflow( flags=flags_n)
		
	stats = calculateMass(grid=test_real_grid).split(",")
	mantaMsg(f"Total \"mass\" inside grid: {stats[0]}, min: {stats[1]}, max: {stats[2]}")

	if (GUI and EXPORT and t % (MAX_TIME / NUM_FRAMES_RENDERED) == 0):
		gui.screenshot( f'../analysis/images/{TITLE.replace(" ", "_")}_{str(t).zfill(4)}.png')

	if s.frame == 0:
		test_start = stats[0]
	test_min = min(test_min, float(stats[0]))
	test_max = max(test_max, float(stats[0]))

	s.step()

mantaMsg(f"Summary of Total Mass:\nStart: {test_start}, Min: {test_min}, Max: {test_max}")