#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            SIMULATION STOP            #
#.......................................#
time.stop_time               =   0.0   # Max (simulated) time to evolve
time.max_step                =   0   # Max number of time steps

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#         TIME STEP COMPUTATION         #
#.......................................#
#time.fixed_dt         =   0.006835938        # Use this constant dt if > 0
time.cfl              =   0.5        # CFL factor

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval  =  100   # Steps between plot files
time.checkpoint_interval =   -1  # Steps between checkpoint files
io.output_default_variables = 0
io.outputs = vorticity
io.derived_outputs = q_criterion q_criterion_nondim mag_vorticity
incflo.post_processing = ke enst
ke.type = KineticEnergy
enst.type = Enstrophy

incflo.do_initial_proj=0
incflo.initial_iterations=0


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
incflo.use_godunov = 1
incflo.godunov_type = "plm"
transport.viscosity = 0.001
transport.laminar_prandtl = 1.0
turbulence.model = Laminar


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        ADAPTIVE MESH REFINEMENT       #
#.......................................#
amr.n_cell              =  256 256 256   # Grid cells at coarsest AMRlevel
#amr.max_grid_size       =  512 512 512 
amr.max_level           =   0           # Max AMR level in hierarchy 
#amr.blocking_factor     =  512 512 512

tagging.labels = vm
tagging.vm.type = VorticityMagRefinement
tagging.vm.nondim = true
tagging.vm.values = 0.1 0.1 0.1 0.1
time.regrid_interval = 10

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              GEOMETRY                 #
#.......................................#
geometry.prob_lo        =   -2. -2. -2.  # Lo corner coordinates
geometry.prob_hi        =   2.  2.  2.  # Hi corner coordinates
geometry.is_periodic    =   0   0   0   # Periodicity x y z (0/1)
xlo.type="pressure_outflow"
ylo.type="pressure_outflow"
zlo.type="pressure_outflow"
xhi.type="pressure_outflow"
yhi.type="pressure_outflow"
zhi.type="pressure_outflow"

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          INITIAL CONDITIONS           #
#.......................................#
incflo.physics = VortexRingCollision
vortexringcollision.R = 1.0
vortexringcollision.Gamma = 1.0
vortexringcollision.delta = 0.1
vortexringcollision.dz = 2.0
vortexringcollision.perturbation_modes = 18
vortexringcollision.perturbation_amplitude = 1e-2
