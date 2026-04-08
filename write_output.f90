module write_output
    use input_reader
    use particle_init
    implicit none
    
    contains
    !abhi ke liye pyvista pe switch ho rhe kyunki paraview ki visualization mein kuch issue aa rhe the, 
    !guyz skill issue ho gya :(  toh abhi ke liye vtk writer mein pyvista ke liye switch kar diya hai, future mein agar paraview ke liye bhi support add karna ho toh uske liye code mein thoda modification karna padega 
    subroutine VtKwriter(iter)
        integer, intent(in) :: iter
        integer :: i
        character(len=256) :: filename
        real(kind=8), allocatable :: positions(:,:), velocities(:,:), heights(:), KE_per_particle(:)


        ! create filename using integer iter to avoid long formatted reals
        write(filename,'(A,I0,A)') 'output_', iter, '.vtk'
        open(unit=20,file=trim(filename), status='unknown', action='write', position='REWIND')
        write(20,'(A)') '# vtk DataFile Version 2.0'
        write(20,'(A,F12.6)') 'Particle data at time ', iter*time_step
        write(20,'(A)') 'ASCII'
        write(20,'(A)') 'DATASET POLYDATA'

        allocate(positions(dimension,nparticles))
        allocate(velocities(dimension, nparticles))
        allocate(heights(nparticles))
        allocate(KE_per_particle(nparticles))
        call get_positions(positions)
        call get_velocities(velocities)
        call particle_height(heights)
        call kinetic_energy(KE_per_particle)

        if (dimension == 3) then
            write(20,'(A,I10,A)') 'POINTS ', nparticles,' float'
            do i=1, nparticles
                write(20,'(3F12.6)') positions(:, i)
            end do
        else if (dimension == 2) then
            write(20,'(A,I10,A)') 'POINTS ', nparticles ,' float'
            do i=1, nparticles
                write(20,'(3F12.6)') positions(1, i), positions(2, i), 0.0d0
            end do
        end if
        
        i = 0
        write(20,'(A,I10)') 'POINT_DATA ', nparticles

        write(20,'(A)') 'SCALARS radius float'
        write(20,'(A)') 'LOOKUP_TABLE default'
        do i=1, nparticles
            write(20,'(F12.6)') rad_init
            end do

        write(20,'(A)') 'SCALARS height float'
        write(20,'(A)') 'LOOKUP_TABLE default'
        do i=1, nparticles
            write(20,'(F12.6)') heights(i)
        end do

        
        write(20,'(A)') 'SCALARS kinetic_energy float'
        write(20,'(A)') 'LOOKUP_TABLE default'
        do i=1, nparticles
            write(20,'(E15.8)') KE_per_particle(i)
        end do

        write(20,'(A)') 'VECTORS velocity float'
        if (dimension == 3) then
            do i=1, nparticles
                write(20,'(3F12.6)') velocities(:, i)
            end do
        else if (dimension == 2) then
            do i=1, nparticles
                write(20,'(3F12.6)') velocities(:, i) , 0.0d0
            end do
        end if

        deallocate(positions, velocities, heights, KE_per_particle)
        close(20)
    end subroutine VtKwriter
end module write_output