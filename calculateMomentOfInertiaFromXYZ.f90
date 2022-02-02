program calculateMomentOfInertiaFromXYZ
implicit none

integer,parameter :: Natoms_max = 200
integer :: Natoms

      ! Atom symbols must be in the correct case
      ! (first letter is capitalized)
      character(2),dimension(26) :: atom_symbols = (/ &
           "H ",                                                 "He", &
           "Li","Be",                   "B ","C ","N ","O ","F ","Ne", &
           "Na","Mg",                   "Al","Si","P ","S ","Cl","Ar", &
           "K ","Ca",                   "Ga","Ge","As","Se","Br","Kr" &
                            /)
      real,dimension(26) :: atomic_masses = (/ &
           1.008,                                               4.003, &
           6.941, 9.012,    10.811,12.011,14.007,15.999,18.998,20.179, &
          22.990,24.305,    26.982,28.086,30.974,32.066,35.453,39.948, &
          39.098,40.078,    69.723,72.590,74.922,78.960,79.904,83.800 &
                            /)

character(2),dimension(Natoms_max) :: symbol
double precision,dimension(3,Natoms_max) :: coords, coordsCM
real,dimension(Natoms_max) :: mass
double precision, dimension(3,3) :: I33, Ir33
double precision :: rii, rij
double precision, dimension(3) :: rCM, l
double precision, dimension(1,1) :: RE1, RE2, RE3, RE

double precision :: ABSTOL, VL, VU
integer, parameter :: LWMAX = 1000
double precision, dimension(LWMAX) :: WORK

double precision, dimension(3) :: eigenvalues
double precision, dimension(3) :: W
double precision, dimension(3,3) :: Z, AF, Z1, Z2, Z3
double precision :: RCOND
double precision, dimension(1) :: FERR, BERR
integer,dimension(3) :: IFAIL
integer,dimension(5*3) :: IWORK
integer :: INFO, LWORK, M

integer,parameter :: filechannel1 = 666
character(200) :: filename
logical :: exists
integer :: iostate

integer :: Nstep

integer :: i, j, k

!
! This program assumes that the XYZs have atoms in the same
! order (e.g. H, Br, H, Cl order throughout the whole file)
!

!
! To compile:
! module load toolchain/intel/2018.5.274
! gfortran -c ls_rmsd_original_v.f90; gfortran ls_rmsd_original_v.o calculateMomentOfInertiaFromXYZ.f90 -o calculateMomentOfInertiaFromXYZ.o -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Nstep = 0
do
  read(5,FMT=*,iostat=iostate) Natoms
  if (iostate /= 0) exit
  read(5,FMT=*)

  ! Read in the coordinates
  do i = 1, Natoms
    read(5,FMT=*) symbol(i), coords(1:3,i)
  end do
  Nstep = Nstep + 1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Assign masses to each of the atoms
  if (Nstep == 1) then
    do i = 1, Natoms
      j = 0
      do k = 1, 26
        if (atom_symbols(k) == symbol(i)) then
          j = k
          exit
        end if
      end do
  
      if (j == 0) then
        write(0,FMT="(A)") "Error. Unkown atom symbol read in XYZ file!"
        write(0,FMT="(A)") "  Hint: Atom symbols are case sensitive"
        stop
      end if

      mass(i) = atomic_masses(j)
    end do
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the center of mass vector
  rcm = 0.0d0
  do i = 1, Natoms
    rcm = rcm + coords(1:3,i) * mass(i)
  end do
  rcm = rcm / sum(mass(1:Natoms))

  do i = 1, Natoms
    coordsCM(1:3,i) = coords(1:3,i) - rcm
  end do

  ! Calculate the moment of inertia tensor
  I33 = 0.0d0
  do i = 1, Natoms

    Ir33 = 0.0d0
    rii = sum(coordsCM(1:3,i)*coordsCM(1:3,i))
    do j = 1, 3
      Ir33(j,j) = Ir33(j,j) + rii
      do k = j+1, 3
        rij = coordsCM(j,i) * coordsCM(k,i)
        Ir33(k,j) = Ir33(k,j) - rij
        Ir33(j,k) = Ir33(j,k) - rij
      end do
    end do

    I33 = I33 + Ir33 * mass(i)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LWORK = 3 * 3

  call DSYSVX( "N", "U", 3, 1, I33(1:3,1:3),&
               3, AF(1:3,1:3), 3, Z(1:3,1:3),&
               reshape(l(1:3),(/3,1/)),&
               3, W, 3,&
               RCOND, FERR, BERR, WORK, LWORK,&
               IWORK, INFO )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !ABSTOL means using the default value
  ABSTOL = -1.0
  
  !Query the optimal workspace.
  LWORK = -1
  
  !Solve eigenproblem.
  CALL DSYEVX('Vectors', 'Indices', 'Upper', 3,&
              I33(1:3,1:3),&
              3, VL, VU, 1,&
              3, ABSTOL, M, eigenvalues(1:3), &
              Z(1:3,1:3), &
              3, WORK, LWORK, &
              IWORK(1:5*3), IFAIL(1:3), &
              INFO )
  
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! write(6,FMT="(1000(ES16.6,1x))") reshape(I33(1:3,1:3),(/9/))
  !Solve eigenproblem.
  CALL DSYEVX('Vectors', 'Indices', 'Upper', 3, &
              I33(1:3,1:3),&
              3, VL, VU, 1,&
              3, ABSTOL, M, eigenvalues(1:3), &
              Z(1:3,1:3), &
              3, WORK, LWORK, &
              IWORK(1:5*3), IFAIL(1:3), &
              INFO )
  ! By choosing 'upper' the subroutine preserves the
  ! upper triangular portion of I33 but messes up the
  ! rest; thus, we correct it here:
  do i = 1, 2
    do j = i+1, 3
      I33(i,j) = I33(j,i)
    end do
  end do
! write(6,FMT="(1000(ES16.6,1x))") reshape(I33(1:3,1:3),(/9/))
  
  !Check for convergence.
  IF( INFO.GT.0 ) THEN
    WRITE(6,*)'The algorithm failed to compute eigenvalues.'
    STOP
  END IF
  
  !Print the number of eigenvalues found.
  !WRITE(6,'(/A,I2)')' The total number of eigenvalues found: ', M
  
  !Print the eigenvalues:
  !WRITE(6,'(A,100(F12.6,1x))') ' Eigenvalues: ', eigenvalues(1:3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  write(6,FMT="(3(ES16.6,1x))")&
          eigenvalues(1:3)

end do

return
end program calculateMomentOfInertiaFromXYZ
