module initialize_module  
  implicit none
  public :: initializeSpins
contains
  subroutine initializeSpins(ndim, spins)
    integer, intent(in) :: ndim
    integer, intent(inout) :: spins(ndim, ndim)
    integer :: ix, iy
    real*8 :: r

    ! initialize spins randomly
    spins = 0
    do ix = 1, ndim
      do iy = 1, ndim
        call random_number(r)
        if ( r > 0.5d0 ) then
          spins(ix, iy) = 1
        else 
          spins(ix, iy) = -1
        end if
      end do
    end do
    
  end subroutine initializeSpins

end module initialize_module  
