module mcLoop_module
  implicit none
  public :: mcLoop, ergi
contains
  subroutine mcLoop(ndim, mcsteps, label, spins,jval,kt,file_num)
    ! variables 
    integer, intent(in) :: ndim, mcsteps, label, file_num
    integer, intent(inout) :: spins(ndim, ndim)
    integer :: nsteps, ix, iy, istep, ndim2
    integer,parameter :: eq_step = 50000
    real*8, intent(in) :: jval,kt
    real*8 :: de, ergi0, ergi1, erg
    real*8 :: r, ptemp, ptrans
    
    ! mc loop
    ndim2 = ndim * ndim
    nsteps = ndim2 * mcsteps
  do istep = 1, nsteps
    ! pick one site randomly
    call random_number(r)
    ix = floor(r*ndim) + 1
    call random_number(r)
    iy = floor(r*ndim) + 1
    ! decide whether flip or not
    ergi0 = ergi(ndim, spins, jval, ix, iy)
    spins(ix, iy) = -spins(ix, iy)
    ergi1 = ergi(ndim, spins, jval, ix, iy)
    de = ergi1 - ergi0
    ptemp = exp(-de/kT)
    ptrans = min(1.0d0, ptemp)
    call random_number(r)
    if ( r < ptrans ) then
      erg = erg + de
    else 
      spins(ix, iy) = -spins(ix, iy)
    end if
    if (istep >= eq_step .and. mod(istep,5000) == 0) then 
      write(file_num, '(F5.2, ",", *(I2, ","), I1)') kt, &
      ((spins(ix, iy), iy=1,ndim), ix=1, ndim), &
      label
    endif
  end do
  end subroutine mcLoop


  function ergi(ndim, spins, jval, ix, iy) 
  integer, intent(in) ::ndim, ix, iy,spins(ndim, ndim)
  real*8, intent(in) :: jval
  real*8 :: ergi
  integer:: jx, jy

  ergi = 0.0d0  
! left
  jx = ix - 1
  if (jx <= 0) jx = ndim
  jy = iy
  ergi = ergi - jval * dble(spins(ix, iy))*dble(spins(jx, jy))
!right
  jx = ix + 1
  if (jx > ndim) jx = 0
  jy = iy
  ergi = ergi - jval * dble(spins(ix, iy))*dble(spins(jx, jy))
!upper
  jy = iy + 1
  if (jy > ndim) jy = 0
  jx = ix
  ergi = ergi - jval * dble(spins(ix, iy))*dble(spins(jx, jy))
!lower
  jy = iy - 1
  if (jy <= 0) jy = ndim
  jx = ix
  ergi = ergi - jval * dble(spins(ix, iy))*dble(spins(jx, jy))

  end function ergi


end module mcLoop_module
