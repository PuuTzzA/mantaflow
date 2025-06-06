#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

RESOLUTION = 64
TITLE = "FLIP simple high CFL"
FILENAME = f'../analysis/data/{TITLE.replace(" ", "_")}.json'
CFL = 1
MAX_TIME = 300
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
s.cfl         = 30            # maximal velocity per cell and timestep, 3 is fairly strict
s.frameLength = 0.8                 # length of one frame (in "world time")
s.timestep    = s.frameLength 
s.timestepMin = 0.001 
s.timestepMax = 2

# prepare grids and particles
flags    = s.create(FlagGrid)
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
flags.initDomain( boundaryWidth=bWidth )
flags.fillGrid()

obsPos = vec3(0.5, 0.1, 0)
obsVelVec = vec3(0.6,0.2,0.0) * (1./100.) * float(res) # velocity in grid units for 100 steps
obsSize = 0.11

obs = Box( parent=s, p0=gs*vec3(0.5,0,0), p1=gs*vec3(0.6,0.1,1))
#obs = Sphere( parent=s, center=gs*obsPos, radius=res*obsSize)
phiObs = obs.computeLevelset()

setObstacleFlags(flags=flags, phiObs=phiObs) 

flags.fillGrid()

#obs.applyToGrid(grid=density, value=0.) # clear smoke inside, flags

if doOpen:
	setOpenBound( flags, bWidth,'yY',FlagOutflow|FlagEmpty )

# enable one of the following
fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.25, 0.5, 1)) # breaking dam
#fluidbox = Box( parent=s, p0=gs*vec3(0.4,0.72,0.4), p1=gs*vec3(0.6,0.92,0.6)) # centered falling block
phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 )

	
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
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) # advect with velocities stored in vel

	if False:
		mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) # maps velocity from particles to grid
		extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) # idk

	else:
		#advectSemiLagrange( flags=flags, vel=vel, grid=vel,   order=2 )
		massMomentumConservingAdvectWater( flags=flags, vel=vel, grid=vel, gammaCumulative=vel_gamma)
	
	markFluidCells( parts=pp, flags=flags )

	if doOpen:
		resetOutflow( flags=flags)


	addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))

	if EXPORT:
		pp.getCurrentData(FILENAME, TITLE, CFL, RESOLUTION, flags=flags, vel=vel, lastFrame=t == (MAX_TIME - 1))

	# pressure solve
	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)

	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	extrapolateMACSimple( flags=flags, vel=vel )
	
	# FLIP velocity update
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 ) # add difference to per particle velocity
	
	#gui.screenshot( 'flipt_%04d.png' % t );
	if (EXPORT and t % (MAX_TIME / NUM_FRAMES_RENDERED) == 0):
		gui.screenshot( f'../analysis/images/{TITLE.replace(" ", "_")}_{str(t).zfill(4)}.png')

	s.step()

