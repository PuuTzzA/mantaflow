#
# Simple example of a free-surface simulation with level-set
# (Optionally, demos outflow boundaries)
# 
from manta import *

# solver params
dim = 2
res = 64
#res = 128 
gs = Vec3(res,res,res)
if (dim==2):
    gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep  = 0.15 

#s.frameLength = 0.6   # length of one frame (in "world time")
#s.timestepMin = 0.01   # time step range
#s.timestepMax = 2.0
#s.cfl         = 0.1   # maximal velocity per cell
#s.timestep    = (s.timestepMax+s.timestepMin)*0.5

# scene file params
ghostFluid  = True
doOpen      = False
accuracy    = 5e-4
# using fast marching is more accurate, but currently causes asymmetries
useMarching = False

visualizerGrid = s.create(LevelsetGrid)

# prepare grids and particles
phi = s.create(LevelsetGrid)
phi_n_plus_one = s.create(LevelsetGrid)
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

flags_n_plus_one = s.create(FlagGrid)
flags_all_fluid = s.create(FlagGrid)
vel_gamma = s.create(MACGrid)
phi_gamma = s.create(RealGrid)
phi_nutzlos = s.create(LevelsetGrid)

fillWithOnes( grid=vel_gamma )

# scene setup
bWidth=1
if False: 
    flags.initDomain(boundaryWidth=bWidth)
    basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(1,0.2,1))
    #basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(0.2,1,1)) # Y
    drop  = Sphere( parent=s , center=gs*Vec3(0.5,0.5,0.5), radius=res*0.125)
    phi.setConst(1e10)
    phi.join(basin.computeLevelset())
    phi.join(drop.computeLevelset())
    flags.updateFromLevelset(phi)
else:
    flags.initDomain(boundaryWidth=bWidth)

    obs = Box( parent=s, p0=gs*vec3(0.5,0,0), p1=gs*vec3(0.6,0.1,1))
    setObstacleFlags(flags=flags, phiObs=obs.computeLevelset()) 
    flags.fillGrid()

    fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.25, 0.5, 1)) # breaking dam
    fluidbox.computeLevelset()

    phi.setConst(1e10)
    phi.join(fluidbox.computeLevelset())

    flags.updateFromLevelset(phi)


# optionally, enable open boundaries here and below...
if doOpen:
    setOpenBound(flags,bWidth,'xXzZ',FlagOutflow|FlagEmpty) 
        
if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()
    

flags_n_plus_one.copyFrom(flags)
flags_all_fluid.copyFrom(flags_n_plus_one)
setAllEmptyFlagsToLiquid(flags=flags_all_fluid)

phi_n_plus_one.copyFrom(phi)
#main loop
for t in range(1000):
    #maxVel = vel.getMax()
    #s.adaptTimestep( maxVel )

    mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))
    
    # update and advect levelset
    if useMarching:
        phi.reinitMarching(flags=flags, velTransport=vel, ignoreWalls=True) 
    else:
        extrapolateLsSimple(phi=phi_n_plus_one, distance=10, inside=False)
        extrapolateLsSimple(phi=phi_n_plus_one, distance=10, inside=True )
        extrapolateMACSimple( flags=flags, vel=vel, distance=10 )

    advectSemiLagrange(flags=flags, vel=vel, grid=phi_n_plus_one, order=2, clampMode=2) 
    #massMomentumConservingAdvect( flags=flags_all_fluid, vel=vel, grid=phi_n_plus_one, gammaCumulative=phi_gamma,interpolationType=0, tracingMethod=0)
    #simpleSLAdvect(flags=flags_all_fluid, vel=vel, grid=phi, interpolationType=0, tracingMethod=0) # 0 = Trilinear, 1 = Cubic, 2= Polynomial Interpolation, 3 = monotonue cubib (hermite)

    phi.setBound(bWidth, 1.) # enforce outside values at border
    if doOpen:
        resetOutflow(flags=flags,phi=phi_n_plus_one) # open boundaries
    
    flags_n_plus_one.updateFromLevelset(phi_n_plus_one)
    
    # velocity self-advection
    #advectSemiLagrange(flags=flags_n_plus_one, vel=vel, grid=vel, order=2)

    massMomentumConservingAdvectWater( flags_n=flags, flags_n_plus_one=flags_n_plus_one, vel=vel, grid=vel, gammaCumulative=vel_gamma, 
                                   phi=phi, interpolationType=0, phi_n_plus_one=phi_n_plus_one)

    #simpleSLAdvect(flags=flags_n_plus_one, vel=vel, grid=vel, interpolationType=0, tracingMethod=0) # 0 = Trilinear, 1 = Cubic, 2= Polynomial Interpolation, 3 = monotonue cubib (hermite)

    visualizeFlags(flags=flags, grid=visualizerGrid, flags_n_plus_one=flags_n_plus_one)

    flags.copyFrom(flags_n_plus_one)
    phi.copyFrom(phi_n_plus_one)
    
    interestedPos = Vec3(16, 2, 0)
    interestedPos2 = Vec3(16, 1, 0)
    print("right after advection")
    #printAtMACPos(vel=vel, pos=interestedPos)

    addGravity(flags=flags, vel=vel, gravity=Vec3(0,-0.025,0))
    #addGravity(flags=flags, vel=vel, gravity=Vec3(-0.025,0,0)) # Y
    
    print("after adding Gravity")
    #printAtMACPos(vel=vel, pos=interestedPos)
    # pressure solve

    extrapolateMACSimple( flags=flags, vel=vel, distance=10 )
    setWallBcs(flags=flags, vel=vel)

    print("after Extrapolation and Boundry conditions")
    #printAtMACPos(vel=vel, pos=interestedPos)
    #printAtMACPos(vel=vel, pos=interestedPos2)
    #printAtMACPos(vel=vel, pos=interestedPos + Vec3(1, 0, 0))
    #printAtMACPos(vel=vel, pos=interestedPos + Vec3(0, 1, 0))
    #printAtMACPos(vel=vel, pos=interestedPos + Vec3(-1, 0, 0))

    if ghostFluid and False:
        solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy, phi=phi )
    else:
        solvePressure(flags=flags, vel=vel, pressure=pressure, cgMaxIterFac=0.5, cgAccuracy=accuracy)
    
    print("after Pressure Projection")
    #printAtMACPos(vel=vel, pos=interestedPos)
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering (only in 3D for now)
    if (dim==3):
        phi.createMesh(mesh)
        #mesh.save('phi%04d.bobj.gz' % t)
    
    findPosVelY(vel=vel, vis=phi_nutzlos)
    s.step()
    #gui.screenshot( 'freesurface_%04d.png' % t );

