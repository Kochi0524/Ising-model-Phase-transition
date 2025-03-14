program ising2d
  use mcLoop_module
  use initialize_module 

  implicit none
  ! variables
  integer,parameter :: mcsteps = 1000
  integer,parameter :: ndim = 20
  integer :: label, spins(ndim, ndim), itemp, file_num
  integer, parameter :: ntrain = 30
  integer, parameter :: nval = 20
  real*8 :: kT ,r
  real*8, parameter :: jval = 1.0d0
  
  ! make train data
  file_num = 10
  open(file_num, file="train.csv", status="replace")
  do itemp = 1, ntrain
    call random_number(r)
    if (itemp < ntrain / 2) then 
      kt = 2.0d0 * r + 0.1d0
      label = 1
      call initializeSpins(ndim, spins)
      call mcLoop(ndim, mcsteps, label, spins,jval, kt,file_num)
    else 
      kt = 2.5d0 + (10.0d0 - 2.5d0) * r
      label = 0
      call initializeSpins(ndim, spins)
      call mcLoop(ndim, mcsteps, label, spins,jval, kt,file_num)
    endif
  end do
  close(file_num)

  ! make validation data
  file_num = 11
  open(file_num, file="val.csv", status="replace")
  do itemp = 1, nval
    call random_number(r)
    kt = 1.0d0 + (3.3d0 - 1.0d0) * r
    if (kt < 2.27d0) then 
      label = 1
      call initializeSpins(ndim, spins)
      call mcLoop(ndim, mcsteps, label, spins,jval, kt,file_num)
    else 
      label = 0
      call initializeSpins(ndim, spins)
      call mcLoop(ndim, mcsteps, label, spins,jval, kt,file_num)
    endif
  end do
  close(file_num)

end program ising2d
