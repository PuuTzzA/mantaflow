#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

RESOLUTION = 64
TITLE = "FLIP simple high CFL"
FILENAME = f'../analysis/data/{TITLE.replace(" ", "_")}.json'
CFL = 1
MAX_TIME = 1000
NUM_FRAMES_RENDERED = 6
EXPORT = False

doOpen = False

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
s.cfl         = 1           # maximal velocity per cell and timestep, 3 is fairly strict
s.frameLength = 0.8                 # length of one frame (in "world time")
s.timestep    = s.frameLength 
s.timestepMin = 0.001 
s.timestepMax = 200

# prepare grids and particles
flags_n    = s.create(FlagGrid)
flags_n_plus_one = s.create(FlagGrid)
vel      = s.create(MACGrid)
vel_gamma = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)

pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 

# scene setup
bWidth=1
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

#obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

if doOpen:
	setOpenBound( flags_n, bWidth,'yY',FlagOutflow|FlagEmpty )
	setOpenBound( flags_n_plus_one, bWidth, 'yY', FlagOutflow|FlagEmpty )

# enable one of the following
fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.25, 0.5, 1)) # breaking dam
#fluidbox = Box( parent=s, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags_n.updateFromLevelset(phiInit)
flags_n_plus_one.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags_n, parts=pp, discretization=particleNumber, randomness=0.2 )

if (GUI):
	gui = Gui()
	gui.show()
	gui.windowSize(800, 800)
	gui.setCamPos(0, 0, -1.3)
	gui.nextRealDisplayMode()
	gui.nextRealDisplayMode()
	gui.nextRealDisplayMode()
	#gui.pause()

pp.clearFile(FILENAME)
#main loop
for t in range(MAX_TIME):
	if s.frame < 1:
		fillWithOnes( grid=vel_gamma )

	maxVel = vel.getMax()
	s.adaptTimestep( maxVel )

	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# FLIP 
	if False:
		pp.advectInGrid(flags=flags_n, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) # advect with velocities stored in vel
		mapPartsToMAC(vel=vel, flags=flags_n, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) # maps velocity from particles to grid
		extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) # idk
		markFluidCells( parts=pp, flags=flags_n )

	else:
		#pp.advectInMACGrid(vel=vel)
		pp.advectInGrid(flags=flags_n, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) # advect with velocities stored in vel
		#advectParticlesForward( particles=pp, vel=vel)
		
		markFluidCells( parts=pp, flags=flags_n_plus_one)

		#advectSemiLagrange( flags=flags_n, vel=vel, grid=vel,   order=2 )
		massMomentumConservingAdvectWater( flags_n=flags_n, flags_n_plus_one=flags_n_plus_one, vel=vel, grid=vel, gammaCumulative=vel_gamma)

		flags_n.copyFrom(flags_n_plus_one)
	
	if doOpen:
		resetOutflow( flags=flags_n)


	addGravity(flags=flags_n, vel=vel, gravity=(0,-0.002,0))

	if EXPORT:
		pp.getCurrentData(FILENAME, TITLE, CFL, RESOLUTION, flags=flags_n, vel=vel, lastFrame=t == (MAX_TIME - 1))

	# pressure solve
	setWallBcs(flags=flags_n, vel=vel)    
	solvePressure(flags=flags_n, vel=vel, pressure=pressure)
	setWallBcs(flags=flags_n, vel=vel)

	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	extrapolateMACSimple( flags=flags_n, vel=vel )
	
	# FLIP velocity update
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags_n, parts=pp, partVel=pVel, flipRatio=0.97 ) # add difference to per particle velocity
	
	#gui.screenshot( 'flipt_%04d.png' % t );
	if (EXPORT and t % (MAX_TIME / NUM_FRAMES_RENDERED) == 0):
		gui.screenshot( f'../analysis/images/{TITLE.replace(" ", "_")}_{str(t).zfill(4)}.png')

	s.step()

