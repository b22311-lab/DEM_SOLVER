module domain_mkr
  use input_reader
  implicit none
  real (kind=8),allocatable :: domain_boundaries(:,:)
  


  contains
    subroutine set_domain()
       
        
        allocate(domain_boundaries(2, dimension))
        domain_boundaries = 0.0d0
        if(dimension == 2) then
            domain_boundaries(2, 1) =  xr_in
            domain_boundaries(2, 2) =  yr_in
        else if (dimension == 3) then
            domain_boundaries(2, 1) =  xr_in
            domain_boundaries(2, 2) =  yr_in
            domain_boundaries(2, 3) =  zr_in
        end if

    end subroutine set_domain


end module domain_mkr