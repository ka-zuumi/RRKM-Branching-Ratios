program numericallyCalculateBranchingRatio
implicit none

character(100) :: intText
integer :: iostate

integer :: Nmolecules
double precision, allocatable :: rates(:,:)
double precision, allocatable :: x0(:,:), x(:,:), dx(:,:), dxdt(:,:)

double precision :: dt, dxdtThreshold

integer :: i, j, k

! Example usage:
! gfortran numericallyCalculateBranchingRatio.f90 -o numericallyCalculateBranchingRatio.o
! file=rateConstantMatrix.dat; ./numericallyCalculateBranchingRatio.o $(sed 1,2d $file | wc -l) < <(sed 1,2d $file | awk '{$1=""; print $0}')

! The threshold for when deciding the system
! is in "steady state" (so we can stop)
dxdtThreshold = 1.0d-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (iargc() /= 1) then
  write(0,FMT="(A)") "Error. Needs exactly one argument:"
  write(0,FMT="(A)") "    (1) Number of molecules"
  stop
end if

call getarg(1,intText)

read(intText,FMT=*,iostat=iostate) Nmolecules
if (iostate /= 0) then
  write(0,FMT="(A)") "Error. First argument is incorrect:"
  write(0,FMT="(A)") "    (1) Number of molecules"
  stop
end if

! The input has, for each row, the rate constant to
! isomerize from molecule ROW to molecule COLUMN
allocate(rates(Nmolecules,Nmolecules))
do i = 1, Nmolecules
  read(5,FMT=*,iostat=iostate) rates(i,1:Nmolecules)
  if (iostate /= 0) then
    write(0,FMT="(A)") "Error. Input is incorrectly formatted"
    write(0,FMT="(A)") "       Must be rate constant matrix"
    stop
  end if
end do

allocate(x0(Nmolecules,1),x(Nmolecules,1),&
         dx(Nmolecules,1),dxdt(Nmolecules,1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the total rate for each species
do i = 1, Nmolecules
  rates(i,i) = -1.0 * sum(rates(i,1:Nmolecules))
end do

rates(1:Nmolecules,1:Nmolecules) = &
  transpose(rates(1:Nmolecules,1:Nmolecules))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For each molecule, simulate what happens if we have
! 1 mole (or unitless quantity) of that species
! initially
do i = 1, Nmolecules

  ! Decide on the correct time stemp for
  ! numerical integration (initially)
  dt = 1.0d-4 / maxval(rates(1:Nmolecules,1:Nmolecules))

  ! Set the initial concentrations
  x0 = 0.0d0
  x0(i,1) = 1.0d0

  ! Integrate forward in time
  x = x0
  j = 0
  do

    dxdt(1:Nmolecules,1:1) = &
            matmul(rates(1:Nmolecules,1:Nmolecules),&
                   x(1:Nmolecules,1:1))

!   x = x + dxdt*dt
!   x = max(x + dxdt*dt, 0.0d0)
!   dx = max(dxdt*dt-x)
    dx = dxdt*dt
    do k = 1, Nmolecules
      if (dx(k,1) < -x(k,1)) then
        dx = - dx * x(k,1) / dx(k,1)
      end if
    end do
    x = x + dx
    j = j + 1

    ! If the system is in "steady state", exit
    if (maxval(abs(dxdt(1:Nmolecules,1))) < dxdtThreshold) exit

!   print *, i, sum(x(1:Nmolecules,1)), maxval(abs(x(1:Nmolecules,1)))
    if (j > 1000) then
      !dt = max(dt, 1.0d-6 * maxval(abs(x(1:Nmolecules,1)/dxdt(1:Nmolecules,1))))

      dt = max(dt, 1.0d-4 / maxval(abs(dxdt(1:Nmolecules,1))))
      !dt = 1.0d-8 / maxval(abs(dxdt(1:Nmolecules,1)))
!     print *, i, sum(x(1:Nmolecules,1)), maxval(x(1:Nmolecules,1)), maxloc(x(1:Nmolecules,1)),&
!             maxval(abs(dxdt(1:Nmolecules,1))), dt
      j = 0
    end if
  end do

  ! Output the final concentration
  ! or the "branching ratios" if the inital
  ! concentration was "1"
  write(6,FMT="(200(F12.8,1x))") x(1:Nmolecules,1)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(rates)

deallocate(x0,x,dx,dxdt)

return
end program numericallyCalculateBranchingRatio
