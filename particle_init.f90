module particle_init
    use input_reader
    use domain_mkr
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private

    public density, radius, k, gamma
    public :: apply_body_force, get_velocities, get_forces, get_positions, get_mass, free_arrays, radius_init, &
              def_material, initialize, contact_particle, contact_wall, update_velocities, update_positions, &
              kinetic_energy, particle_height, compute_particle_contacts, compute_wall_contacts, update_all_particles, &
              reset_contact_search_counters, get_contact_search_counters

    real(kind=8), allocatable :: mass(:), radius(:), velocities(:,:), positions(:,:), force(:,:)
    integer(kind=8) :: contact_candidates_total = 0_8
    integer(kind=8) :: contact_contacts_total = 0_8

contains
    pure real(kind=8) function contact_force_magnitude(stiffness, damping, overlap, normal_velocity)
        real(kind=8), intent(in) :: stiffness, damping, overlap, normal_velocity
        contact_force_magnitude = max(0.0d0, stiffness * overlap - damping * normal_velocity)
    end function contact_force_magnitude

    subroutine reset_contact_search_counters()
        contact_candidates_total = 0_8
        contact_contacts_total = 0_8
    end subroutine reset_contact_search_counters

    subroutine get_contact_search_counters(total_candidates, total_contacts)
        integer(kind=8), intent(out) :: total_candidates, total_contacts
        total_candidates = contact_candidates_total
        total_contacts = contact_contacts_total
    end subroutine get_contact_search_counters

!subroutine for initializing the particles with user defined radius and material properties    
    subroutine initialize(tstep)    
        integer , intent(in) :: tstep
        if (tstep ==1) then
            allocate ( velocities(dimension, nparticles), positions(dimension, nparticles), force(dimension, nparticles))
            velocities = 0.0d0
            positions = 0.0d0
            force = 0.0d0

            if (nparticles == 1 .or. any(abs(posit_init) > 1.0d-14)) then
                positions = posit_init
            else
                call initialize_positions_lattice()
            end if

            velocities = vel_init
        end if

        force = 0.0d0

        if (verbose) then
            print *, "Particles initialized with zero forces........."
        end if
    end subroutine initialize

    subroutine initialize_positions_lattice()
        integer :: i, ix, iy, iz
        integer :: nx, ny, nz
        real(kind=8) :: spacing
        real(kind=8) :: x_min, x_max, y_min, y_max, z_min, z_max

        spacing = 2.2d0 * rad_init

        x_min = domain_boundaries(1, 1) + rad_init
        x_max = domain_boundaries(2, 1) - rad_init
        y_min = domain_boundaries(1, 2) + rad_init
        y_max = domain_boundaries(2, 2) - rad_init

        if (dimension == 2) then
            nx = int((x_max - x_min) / spacing) + 1
            ny = int((y_max - y_min) / spacing) + 1

            if (nx * ny < nparticles) then
                error stop "Not enough domain volume to place all particles on lattice in 2D."
            end if

            i = 0
            do iy = 0, ny - 1
                do ix = 0, nx - 1
                    i = i + 1
                    if (i > nparticles) exit
                    positions(1, i) = x_min + ix * spacing
                    positions(2, i) = y_min + iy * spacing
                end do
                if (i >= nparticles) exit
            end do

        else if (dimension == 3) then
            z_min = domain_boundaries(1, 3) + rad_init
            z_max = domain_boundaries(2, 3) - rad_init

            nx = int((x_max - x_min) / spacing) + 1
            ny = int((y_max - y_min) / spacing) + 1
            nz = int((z_max - z_min) / spacing) + 1

            if (nx * ny * nz < nparticles) then
                error stop "Not enough domain volume to place all particles on lattice in 3D."
            end if

            i = 0
            do iz = 0, nz - 1
                do iy = 0, ny - 1
                    do ix = 0, nx - 1
                        i = i + 1
                        if (i > nparticles) exit
                        positions(1, i) = x_min + ix * spacing
                        positions(2, i) = y_min + iy * spacing
                        positions(3, i) = z_min + iz * spacing
                    end do
                    if (i >= nparticles) exit
                end do
                if (i >= nparticles) exit
            end do
        end if
    end subroutine initialize_positions_lattice

    subroutine radius_init(n)
        integer , intent(in) :: n
        allocate(radius(n))
       
        radius = rad_init
    end subroutine radius_init

    subroutine def_material()
        allocate(mass(nparticles))
        if (verbose) then
            print *, "Defining material properties for the particles........"
        end if

        mass = density * (4.0d0/3.0d0) * 3.14159265358979323846 * radius**3
   
    end subroutine def_material

! implementing getter function for positions, velocities, forces, mass

    pure subroutine get_velocities(array)
        real(kind=8), intent(out) :: array(:,:)
        array = velocities
    end subroutine get_velocities

    pure subroutine get_forces(array)
        real(kind=8), intent(out) :: array(:,:)
        array = force
    end subroutine get_forces

    pure subroutine get_positions(array)
        real(kind=8), intent(out) :: array(:,:)
        array = positions
    end subroutine get_positions

    pure subroutine get_mass(array)
        real(kind=8), intent(out) :: array(:)
        array = mass
    end subroutine get_mass
    


    !implemeting body force function for the particles
    subroutine apply_body_force()
        integer :: j

        !$omp parallel do private(j) shared(nparticles, gravity_axis, force, mass, gravity_force) if(nparticles > 32)
        do j = 1, nparticles
            force(gravity_axis, j) = force(gravity_axis, j) + mass(j) * gravity_force(gravity_axis)
        end do
        !$omp end parallel do
    end subroutine apply_body_force

! calculation particle contatct forces 
    subroutine contact_particle(particle_iter)
        integer, intent(in) :: particle_iter
        integer :: i, j
        real(kind=8) :: distance, overlap, normal_velocity, new_mag
        real(kind=8) :: unit_normal(dimension), relative_position(dimension), relative_velocity(dimension)

        do j = particle_iter + 1, nparticles
            do i = 1, dimension
                relative_position(i) = positions(i, j) - positions(i, particle_iter)
            end do

            distance = sqrt(sum(relative_position**2))
            if (distance <= 1.0d-14) cycle

            unit_normal = relative_position / distance
            overlap = radius(particle_iter) + radius(j) - distance

            if (overlap > 0.0d0) then
                relative_velocity = velocities(:, j) - velocities(:, particle_iter)
                normal_velocity = sum(relative_velocity * unit_normal)
                new_mag = contact_force_magnitude(k, gamma, overlap, normal_velocity)
                if (new_mag > 0.0d0) then
                    force(:, particle_iter) = force(:, particle_iter) + new_mag * unit_normal
                    force(:, j) = force(:, j) - new_mag * unit_normal
                end if
            end if
        end do
    end subroutine contact_particle

    subroutine compute_particle_contacts()
        integer(kind=8) :: step_candidates, step_contacts

        step_candidates = 0_8
        step_contacts = 0_8

        select case (contact_search_method)
        case (0)
            call compute_particle_contacts_allpairs(step_candidates, step_contacts)
        case (1)
            call compute_particle_contacts_cell(step_candidates, step_contacts)
        case default
            error stop "Unsupported contact_search_method. Use 0 (all-pairs) or 1 (cell-linked)."
        end select

        contact_candidates_total = contact_candidates_total + step_candidates
        contact_contacts_total = contact_contacts_total + step_contacts
    end subroutine compute_particle_contacts

    subroutine compute_particle_contacts_allpairs(candidates, contacts)
        integer(kind=8), intent(out) :: candidates, contacts
        integer :: i, j, d, tid, nthreads
        real(kind=8) :: distance, overlap, normal_velocity, new_mag
        real(kind=8) :: unit_normal(dimension), relative_position(dimension), relative_velocity(dimension)
        real(kind=8), allocatable :: force_private(:, :, :)

        candidates = 0_8
        contacts = 0_8

        if (nparticles < 2) return

        nthreads = max(1, omp_get_max_threads())
        allocate(force_private(dimension, nparticles, nthreads))
        force_private = 0.0d0

        !$omp parallel &
        !$omp& private(i, j, d, tid, distance, overlap, normal_velocity, new_mag, &
        !$omp& unit_normal, relative_position, relative_velocity) &
        !$omp& shared(nthreads, nparticles, dimension, positions, radius, velocities, force_private, k, gamma) &
        !$omp& reduction(+:candidates, contacts)
        tid = omp_get_thread_num() + 1

        !$omp do schedule(static)
        do i = 1, nparticles - 1
            do j = i + 1, nparticles
                candidates = candidates + 1_8

                do d = 1, dimension
                    relative_position(d) = positions(d, j) - positions(d, i)
                end do

                distance = sqrt(sum(relative_position**2))
                if (distance <= 1.0d-14) cycle

                unit_normal = relative_position / distance
                overlap = radius(i) + radius(j) - distance

                if (overlap > 0.0d0) then
                    contacts = contacts + 1_8
                    relative_velocity = velocities(:, j) - velocities(:, i)
                    normal_velocity = sum(relative_velocity * unit_normal)
                    new_mag = contact_force_magnitude(k, gamma, overlap, normal_velocity)

                    if (new_mag > 0.0d0) then
                        do d = 1, dimension
                            force_private(d, i, tid) = force_private(d, i, tid) + new_mag * unit_normal(d)
                            force_private(d, j, tid) = force_private(d, j, tid) - new_mag * unit_normal(d)
                        end do
                    end if
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel

        do tid = 1, nthreads
            force = force + force_private(:, :, tid)
        end do

        deallocate(force_private)
    end subroutine compute_particle_contacts_allpairs

    subroutine compute_particle_contacts_cell(candidates, contacts)
        integer(kind=8), intent(out) :: candidates, contacts
        integer :: nx, ny, nz, total_cells
        integer :: ix, iy, iz, jx, jy, jz
        integer :: i, j, p, cell_id, neigh_id
        integer, allocatable :: head(:), nextp(:)
        real(kind=8) :: cell_size
        real(kind=8) :: x_min, y_min, z_min
        real(kind=8) :: distance, overlap, normal_velocity, new_mag
        real(kind=8) :: unit_normal(dimension), relative_position(dimension), relative_velocity(dimension)

        candidates = 0_8
        contacts = 0_8

        if (nparticles < 2) return

        cell_size = max(2.0d0 * rad_init, cell_size_factor * rad_init)
        if (cell_size <= 1.0d-14) then
            cell_size = 2.0d0 * rad_init
        end if

        x_min = domain_boundaries(1, 1)
        y_min = domain_boundaries(1, 2)
        nx = max(1, int((domain_boundaries(2, 1) - x_min) / cell_size))
        ny = max(1, int((domain_boundaries(2, 2) - y_min) / cell_size))

        if (dimension == 3) then
            z_min = domain_boundaries(1, 3)
            nz = max(1, int((domain_boundaries(2, 3) - z_min) / cell_size))
        else
            z_min = 0.0d0
            nz = 1
        end if

        total_cells = nx * ny * nz
        allocate(head(total_cells), nextp(nparticles))
        head = 0
        nextp = 0

        if (dimension == 2) then
            do p = 1, nparticles
                ix = int((positions(1, p) - x_min) / cell_size) + 1
                iy = int((positions(2, p) - y_min) / cell_size) + 1

                ix = max(1, min(nx, ix))
                iy = max(1, min(ny, iy))

                cell_id = ix + (iy - 1) * nx
                nextp(p) = head(cell_id)
                head(cell_id) = p
            end do

            do iy = 1, ny
                do ix = 1, nx
                    cell_id = ix + (iy - 1) * nx

                    i = head(cell_id)
                    do while (i > 0)
                        j = nextp(i)
                        do while (j > 0)
                            call evaluate_pair(i, j, candidates, contacts)
                            j = nextp(j)
                        end do
                        i = nextp(i)
                    end do

                    do jy = max(1, iy - 1), min(ny, iy + 1)
                        do jx = max(1, ix - 1), min(nx, ix + 1)
                            neigh_id = jx + (jy - 1) * nx
                            if (neigh_id <= cell_id) cycle

                            i = head(cell_id)
                            do while (i > 0)
                                j = head(neigh_id)
                                do while (j > 0)
                                    call evaluate_pair(i, j, candidates, contacts)
                                    j = nextp(j)
                                end do
                                i = nextp(i)
                            end do
                        end do
                    end do
                end do
            end do

        else if (dimension == 3) then
            do p = 1, nparticles
                ix = int((positions(1, p) - x_min) / cell_size) + 1
                iy = int((positions(2, p) - y_min) / cell_size) + 1
                iz = int((positions(3, p) - z_min) / cell_size) + 1

                ix = max(1, min(nx, ix))
                iy = max(1, min(ny, iy))
                iz = max(1, min(nz, iz))

                cell_id = ix + (iy - 1) * nx + (iz - 1) * nx * ny
                nextp(p) = head(cell_id)
                head(cell_id) = p
            end do

            do iz = 1, nz
                do iy = 1, ny
                    do ix = 1, nx
                        cell_id = ix + (iy - 1) * nx + (iz - 1) * nx * ny

                        i = head(cell_id)
                        do while (i > 0)
                            j = nextp(i)
                            do while (j > 0)
                                call evaluate_pair(i, j, candidates, contacts)
                                j = nextp(j)
                            end do
                            i = nextp(i)
                        end do

                        do jz = max(1, iz - 1), min(nz, iz + 1)
                            do jy = max(1, iy - 1), min(ny, iy + 1)
                                do jx = max(1, ix - 1), min(nx, ix + 1)
                                    neigh_id = jx + (jy - 1) * nx + (jz - 1) * nx * ny
                                    if (neigh_id <= cell_id) cycle

                                    i = head(cell_id)
                                    do while (i > 0)
                                        j = head(neigh_id)
                                        do while (j > 0)
                                            call evaluate_pair(i, j, candidates, contacts)
                                            j = nextp(j)
                                        end do
                                        i = nextp(i)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if

        deallocate(head, nextp)

    contains
        subroutine evaluate_pair(i_idx, j_idx, candidate_count, contact_count)
            integer, intent(in) :: i_idx, j_idx
            integer(kind=8), intent(inout) :: candidate_count, contact_count
            integer :: d

            candidate_count = candidate_count + 1_8

            do d = 1, dimension
                relative_position(d) = positions(d, j_idx) - positions(d, i_idx)
            end do

            distance = sqrt(sum(relative_position**2))
            if (distance <= 1.0d-14) return

            unit_normal = relative_position / distance
            overlap = radius(i_idx) + radius(j_idx) - distance

            if (overlap > 0.0d0) then
                contact_count = contact_count + 1_8
                relative_velocity = velocities(:, j_idx) - velocities(:, i_idx)
                normal_velocity = sum(relative_velocity * unit_normal)
                new_mag = max(0.0d0, k * overlap - gamma * normal_velocity)

                if (new_mag > 0.0d0) then
                    do d = 1, dimension
                        force(d, i_idx) = force(d, i_idx) + new_mag * unit_normal(d)
                        force(d, j_idx) = force(d, j_idx) - new_mag * unit_normal(d)
                    end do
                end if
            end if
        end subroutine evaluate_pair
    end subroutine compute_particle_contacts_cell

    subroutine contact_wall(particle_index, debug)
        integer, intent(in) :: particle_index
        logical, intent(in), optional :: debug
        integer :: j
        real(kind=8) :: overlap, normal_velocity, unit_normal, force_mag
        real(kind=8) :: lower_overlap, upper_overlap

        do j = 1, dimension
            lower_overlap = radius(particle_index) - (positions(j, particle_index) - domain_boundaries(1, j))
            upper_overlap = radius(particle_index) - (domain_boundaries(2, j) - positions(j, particle_index))

            overlap = 0.0d0
            unit_normal = 0.0d0

            if (lower_overlap > 0.0d0) then
                overlap = lower_overlap
                unit_normal = 1.0d0
            else if (upper_overlap > 0.0d0) then
                overlap = upper_overlap
                unit_normal = -1.0d0
            end if

            if (overlap > 0.0d0) then
                normal_velocity = velocities(j, particle_index) * unit_normal
                force_mag = contact_force_magnitude(k, gamma, overlap, normal_velocity)
                force(j, particle_index) = force(j, particle_index) + force_mag * unit_normal
            end if

            if (present(debug)) then
                print *, "wall_overlap=", overlap, " force_component=", force(j, particle_index)
            end if
        end do
    end subroutine contact_wall

    subroutine compute_wall_contacts()
        integer :: i

        !$omp parallel do private(i) shared(nparticles) if(nparticles > 32)
        do i = 1, nparticles
            call contact_wall(i)
        end do
        !$omp end parallel do
    end subroutine compute_wall_contacts

    subroutine update_velocities(particle_index, tstep)
        integer, intent(in) :: particle_index
        real(kind=8), intent(in) :: tstep

        velocities(:, particle_index) = velocities(:, particle_index) + (force(:, particle_index) / mass(particle_index)) * tstep
    end subroutine update_velocities

    subroutine update_positions(particle_index, tstep, debug)
        integer, intent(in) :: particle_index
        real(kind=8), intent(in) :: tstep
        logical, intent(in), optional :: debug
        integer :: j

        positions(:, particle_index) = positions(:, particle_index) + velocities(:, particle_index) * tstep

        do j = 1, dimension
            if (positions(j, particle_index) < domain_boundaries(1, j) .OR. &
                positions(j, particle_index) > domain_boundaries(2, j)) then

                if (present(debug)) then
                    print *, positions(j, particle_index), ".....................", domain_boundaries(:, j)
                end if
                print *, "Warning: Particle ", particle_index, " is out of bounds in dimension ", j
                error stop "Unphysical behavior detected: Particle out of bounds. Simulation terminated."
            end if
        end do
    end subroutine update_positions

    subroutine update_all_particles(tstep)
        real(kind=8), intent(in) :: tstep
        integer :: i

        !$omp parallel do private(i) shared(nparticles, tstep) if(nparticles > 32)
        do i = 1, nparticles
            call update_velocities(i, tstep)
        end do
        !$omp end parallel do

        !$omp parallel do private(i) shared(nparticles, tstep) if(nparticles > 32)
        do i = 1, nparticles
            call update_positions(i, tstep)
        end do
        !$omp end parallel do
    end subroutine update_all_particles

    subroutine kinetic_energy(array)
        real(kind=8), intent(out) :: array(:)
        integer :: i

        !$omp parallel do private(i) shared(nparticles, array, mass, velocities) if(nparticles > 32)
        do i = 1, nparticles
            array(i) = 0.5d0 * mass(i) * sum(velocities(:, i)**2)
        end do
        !$omp end parallel do
    end subroutine kinetic_energy

    pure subroutine particle_height(array)
        real(kind=8), intent(out) :: array(:)
        array(:) = positions(2, :)

    end subroutine particle_height

    ! freeing the allocated arrays
    subroutine free_arrays()
        if (allocated(velocities)) deallocate(velocities)
        if (allocated(positions)) deallocate(positions)
        if (allocated(force)) deallocate(force)
    end subroutine free_arrays

end module particle_init
