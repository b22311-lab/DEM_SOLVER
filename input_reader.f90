module input_reader
    implicit none
    integer :: dimension ,nparticles,write_interval , gravity_axis
    integer :: contact_search_method = 0
    logical :: verbose = .false.
    real (kind=8) :: k, gamma , density, rad_init ,end_time,time_step,&
    xr_in ,yr_in=0.0d0 ,zr_in=0.0d0 ,gravity_magnitude=9.81d0
    real (kind=8) :: cell_size_factor = 2.5d0
    real (kind=8) , allocatable :: posit_init(:, :), vel_init(:, :)
    real (kind=8) :: gravity_force(3)
    namelist /domain/ dimension, xr_in ,yr_in ,zr_in
    namelist /material_properties/ k , gamma , density , rad_init ,nparticles
    namelist /particle_initialization/ posit_init, vel_init
    namelist /simulation_control/ time_step, end_time, write_interval, verbose, contact_search_method, cell_size_factor
    namelist /gravity/ gravity_axis , gravity_magnitude
    
contains
    subroutine read_input(filename)
        character(len=*) ,intent(in) :: filename
        integer :: iostat
        open(unit=10, file=trim(filename), status='old', action='read',iostat=iostat)
        if (iostat /= 0) then

            print *, " FATAL ERROR: Could not find '", trim(filename), "'!"
            print *, " Please ensure the file is in the same directory."
          
            error stop "Simulation aborted due to missing input file."
        end if
     
        read(10, nml=domain)
        if (verbose) print*, 'read domain'
        read(10, nml=material_properties)
        if (verbose) print*, 'read material properties'
        allocate(posit_init(dimension, nparticles))
        allocate(vel_init(dimension, nparticles))
        posit_init = 0.0d0
        vel_init = 0.0d0
        read(10, nml=particle_initialization)
        if (verbose) print*, 'read particle initialization'
        read(10, nml=simulation_control)
        if (verbose) print*, 'read simulation control'
        read(10, nml=gravity)
        if (verbose) print*, 'read gravity properties'

        gravity_force(gravity_axis) = -gravity_magnitude
        !change the particle initilization and velocity intilization method 
        
        close(10)

    end subroutine read_input
    
end module input_reader