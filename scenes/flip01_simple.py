#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

RESOLUTION = 64
FILENAME = "../flip_simple.json"
MAX_TIME = 300

# solver params
dim = 2
particleNumber = 2
res = RESOLUTION
gs = vec3(res,res,res)
if (dim==2):
	gs.z=1
	particleNumber = 3 # use more particles in 2d
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 

# scene setup
flags.initDomain(boundaryWidth=0) 
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
	#gui.pause()

def getFromType(type):
	match type:
		case 0:
			return "TypeNone"
		case 1:
			return "TypeFluid"
		case 2:
			return "TypeObstacle"
		case 4: 
			return "TypeEmpty"
		case 8:
			return "TypeInflow"
		case 16:
			return "TypeOutflow"
		case 32:
			return "TypeOpen"
		case 64:
			return "TypeStick"
	return "InternalType"

pp.clearFile(FILENAME)
#main loop
for t in range(MAX_TIME):
	mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))
	
	# FLIP 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
	markFluidCells( parts=pp, flags=flags )

	addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))

	if 1:
		pp.getCurrentData(FILENAME, RESOLUTION, flags=flags, lastFrame=t == (MAX_TIME - 1))

	# pressure solve
	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)

	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	extrapolateMACSimple( flags=flags, vel=vel )
	
	# FLIP velocity update
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )
	
	#gui.screenshot( 'flipt_%04d.png' % t );
	s.step()

