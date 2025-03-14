program ising2d

  implicit none
  ! variables
  integer,parameter :: mcsteps = 1000
  integer,parameter :: ndim = 20
  integer,parameter :: eq_step = 50000
  integer :: ndim2, nsteps,istep
  integer ::ix, iy 
  integer :: label, spins(ndim, ndim)
  real*8 :: kT 
  real*8 ::erg,ergi0, ergi1,de,mag
  real*8 :: r, ptemp, ptrans
  real*8, parameter :: jval = 1.0d0
  real*8, external :: ergi
  
  ndim2 = ndim* ndim

  ! input kt
  write(6, *) "input kt"
  read(5,*) kt

  if(kt > 2.27d0) then 
    label = 0
  else 
    label = 1
  endif

  ! initial state in hot start
  spins = 0
  do ix = 1, ndim
    do iy = 1, ndim
      call random_number(r)
      if ( r > 0.5d0 ) then 
        spins(ix, iy) = 1
      else 
        spins(ix, iy) = -1
      endif
    end do
  end do


  open(10, file="example.csv", status="replace")

  ! Metropolis MC loop
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
      write(10, '(F5.2, ",", *(I2, ","), I1)') kt, &
      ((spins(ix, iy), iy=1,ndim), ix=1, ndim), &
      label
    endif
  end do

  close(10)

end program ising2d

! calculate energy of per site
function ergi(ndim, spins, jval, ix, iy) 
  implicit none
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
