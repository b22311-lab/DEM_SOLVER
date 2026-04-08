program main_loop
    use particle_init
    use domain_mkr
    use input_reader
    use write_output
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none

    integer :: loop_cont, i
    character(len=256) :: input_file

    real(kind=8) :: t0
    real(kind=8) :: total_start, total_end, total_runtime
    real(kind=8) :: t_init, t_gravity, t_contacts, t_walls, t_update, t_output
    integer(kind=8) :: contact_candidates_total, contact_contacts_total

    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) then
        write(*,*) "Usage: ./demo.exe input_file.in"
        stop
    end if

    call read_input(trim(input_file))

    loop_cont = int(end_time / time_step)

    call radius_init(nparticles)
    call def_material()
    call set_domain()

    t_init = 0.0d0
    t_gravity = 0.0d0
    t_contacts = 0.0d0
    t_walls = 0.0d0
    t_update = 0.0d0
    t_output = 0.0d0

    call reset_contact_search_counters()

    if (verbose) then
        print *, "Running with OpenMP threads =", omp_get_max_threads()
    end if

    total_start = omp_get_wtime()

    do i = 1, loop_cont
        if (verbose) then
            print *, "Solving for time:", i * time_step
        end if

        t0 = omp_get_wtime()
        call initialize(i)
        t_init = t_init + (omp_get_wtime() - t0)

        t0 = omp_get_wtime()
        call apply_body_force()
        t_gravity = t_gravity + (omp_get_wtime() - t0)

        t0 = omp_get_wtime()
        call compute_particle_contacts()
        t_contacts = t_contacts + (omp_get_wtime() - t0)

        t0 = omp_get_wtime()
        call compute_wall_contacts()
        t_walls = t_walls + (omp_get_wtime() - t0)

        t0 = omp_get_wtime()
        call update_all_particles(time_step)
        t_update = t_update + (omp_get_wtime() - t0)

        if (mod(i, write_interval) == 0) then
            t0 = omp_get_wtime()
            call VtKwriter(i)
            t_output = t_output + (omp_get_wtime() - t0)
        end if
    end do

    total_end = omp_get_wtime()
    total_runtime = total_end - total_start
    call get_contact_search_counters(contact_candidates_total, contact_contacts_total)

    print *, "PROFILE total_runtime_s", total_runtime
    print *, "PROFILE initialize_s", t_init
    print *, "PROFILE gravity_s", t_gravity
    print *, "PROFILE particle_contacts_s", t_contacts
    print *, "PROFILE wall_contacts_s", t_walls
    print *, "PROFILE integration_s", t_update
    print *, "PROFILE output_s", t_output
    print *, "PROFILE contact_search_method", contact_search_method
    print *, "PROFILE contact_cell_size_factor", cell_size_factor
    print *, "PROFILE contact_candidates_total", contact_candidates_total
    print *, "PROFILE contacts_detected_total", contact_contacts_total

    if (loop_cont > 0) then
        print *, "PROFILE contact_candidates_avg_per_step", real(contact_candidates_total, kind=8) / real(loop_cont, kind=8)
        print *, "PROFILE contacts_detected_avg_per_step", real(contact_contacts_total, kind=8) / real(loop_cont, kind=8)
    end if

    if (total_runtime > 1.0d-14) then
        print *, "PROFILE pct_initialize", 100.0d0 * t_init / total_runtime
        print *, "PROFILE pct_gravity", 100.0d0 * t_gravity / total_runtime
        print *, "PROFILE pct_particle_contacts", 100.0d0 * t_contacts / total_runtime
        print *, "PROFILE pct_wall_contacts", 100.0d0 * t_walls / total_runtime
        print *, "PROFILE pct_integration", 100.0d0 * t_update / total_runtime
        print *, "PROFILE pct_output", 100.0d0 * t_output / total_runtime
    end if

    call free_arrays()

end program main_loop
