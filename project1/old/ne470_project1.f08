!ne470_project1.f08

!Title:			NE 470 Project #1
!Engineer:		Tom Looby
!Date:			10/03/2017
!Description:	Solves 1D 1-group Neutron Diffusion Equation, given
!               user input parameters, and outputs results to screen
!               and CSV file.

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program ne470_project1
implicit none

!=======================================================================
!            Variables, Literals
!=======================================================================

integer :: i, j, ierror, meshN, geom, init_a
integer, parameter :: dp = selected_real_kind(15, 307)
real(dp) :: del_x, S, S_mult, D, sigma_tr, sigma_a, c, S2, sqr


!=== Initialize arrays based upon user input
real(dp), ALLOCATABLE :: mat_A(:,:), A_inv(:,:)
real(dp), ALLOCATABLE :: mat_B(:), AinvB(:)
real(dp), ALLOCATABLE :: F(:), F_an(:), F_num(:), x(:)
logical :: exists
character(len=200) :: outfile, dummy


!=======================================================================
!           Print Stuff / Instructions
!=======================================================================

print *, '======================Project 1 ============================='
print *, 'Engineers: Group 1'
print *, 'Date: 10/07/2017'
print *, ' '
print *, 'This program accepts user defined geometry, cross sections, '
print *, 'source multiplier, and mesh resolution.  It solves a 2nd '
print *, 'order ODE of the form -D*f" + sigma_a*f = S.  S = S when '
print *, 'x = 0, and S = 0 if x > 0.  This equation is commonly '
print *, 'referred to as the 1D 1-Group Neutron Diffusion Equation. '
print *, ' '
print *, 'The output is in two formats.  1) Direct output to the'
print *, 'screen, and 2) to a CSV file in the local directory.'
print *, 'The user declares CSV name in the input.  '
print *, ' '
print *, 'Please provide necessary program input below...'
print *, '=============================================================' 
print *, ' '



!=======================================================================
!           User Input with Format Error Checking  
!=======================================================================

!===User defined geometry with error checking===
do
  write(*,*) "Enter Geometry from Following Options (integer):"
  write(*,*) "1 Cartesian"
  write(*,*) "2 Cylindrical"
  write(*,*) "3 Spherical"
  
  read(*,'(i10)',iostat=ierror) geom
  if ( ierror == 0 .AND. geom <= 3 .AND. geom >= 1)  then
    exit
  endif
  write(*,*) 'Please enter an integer between 1 and 3.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Geometry Integer: ', geom
write(*,*) ' '
ierror = 0

!===User defined Mesh Node Number with error checking===
do
  write(*,*) "Enter Number of desired mesh points (integer):"

  read(*,'(i10)',iostat=ierror) meshN
  if ( ierror == 0 .AND. meshN >= 2)  then
    exit
  endif
  write(*,*) 'Please enter an integer greater than 2.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Mesh Points: ', meshN
write(*,*) ' '
ierror = 0

!===User defined Source Multiplier with error checking===
do
  write(*,*) "Enter Desired Secondary Source Multiplier:"

  read(*,'(F10.0)',iostat=ierror) S_mult
  if ( ierror == 0 .AND. S_mult >= 0)  then
    exit
  endif
  write(*,*) 'Please enter a numeric greater than 0.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Source Multiplier: ', S_mult
write(*,*) ' '
ierror = 0


!===Absorber / Moderator parameters with error checking===
do
  write(*,*) "Enter Macroscopic Absorbtion Cross Section (cm^-1):"

  read(*,'(F10.0)',iostat=ierror) sigma_a
  if ( ierror == 0 )  then
    exit
  endif
  write(*,*) 'Please enter a numeric.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Absorbtion Macroscopic Cross Section: ', sigma_a
write(*,*) ' '
ierror = 0

do
  write(*,*) "Enter Macroscopic Transport Cross Section (cm^-1):"

  read(*,'(F10.0)',iostat=ierror) sigma_tr
  if ( ierror == 0 )  then
    exit
  endif
  write(*,*) 'Please enter a numeric.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Macroscopic Transport Cross Section: ', sigma_tr
write(*,*) ' '
ierror = 0



!===CSV output with file overwrite error checking===
do
  write(*,*) "Enter CSV Filename (ie test.csv):"
  read(*,*,iostat=ierror) outfile
  if ( ierror == 0 )  then
    exit
  endif
  write(*,*) 'Please enter a string (ie test.csv).  Try Again...'
  write(*,*) ' '
enddo

INQUIRE(File=outfile, Exist=exists)
IF (exists) THEN
   print *, ' ' 
   print *, 'FILE EXISTS!  CSV data will be overwritten if you proceed!'
   print *, 'Press any key then ENTER to Proceed'
   read(*,*,iostat=ierror) dummy
END IF
write(*,*) ' '

!===Initialize dynamic arrays based upon received user input
init_a = meshN - 2
ALLOCATE(mat_A(init_a,init_a), A_inv(init_a,init_a))
ALLOCATE(mat_B(init_a), AinvB(init_a))
ALLOCATE(F(meshN), F_an(meshN), F_num(meshN), x(meshN))

!=======================================================================
!            Literals and Initial Conditions
!=======================================================================
!===Create c value from geometry selection
c = geom - 1
S = (10**8)
S2 = (10**8)*S_mult
!===Diffusion Coefficient Definition
D = 1/(3*sigma_tr)

!===Initial Conditions
F(1) = S/D
F(meshN) = S2/D

!===Step Size (Normalized Against Width ie 1, 2, 3, etc.)
do i=1,(meshN)
   x(i) = i-1
end do




!=======================================================================
!            Analytical Solution         
!=======================================================================
!This section can be populated with analytical solution to 
!create 1 CSV file with both data sets

sqr = SQRT(sigma_a/D)
do i=1,meshN

   F_an(i) = S*(exp(sqr*x(i)) - exp(sqr*(2*meshN - x(i))))/((1-exp(2*meshN*sqr)))

end do



!=======================================================================
!           Numerical Solution             
!=======================================================================

do i=1, (meshN - 2)
   ! Account for X vector index offset
   del_x = x(i+1) - x(i)
   do j=1,(meshN - 2)

      IF (j == i) THEN
         mat_A(i,j) = (((2*D)/(del_x**2)) + sigma_a)

      ELSE IF (j == i+1) THEN
         mat_A(i,j) = (((-D)/(del_x**2))*(1 + c/(2*i - 1)))
                  
      ELSE IF (j == i-1) THEN
         mat_A(i,j) = (((-D)/(del_x**2))*(1 - c/(2*i - 1)))
         ! Test Stuff
         !print *,'i:', i
         !print *,'i - 1 value:', mat_A(i,j)
         !print *,'C Part:', c/(2*i - 1)
         !print *,'c:', c
                  
      ELSE
         mat_A(i,j) = 0

      END IF

   end do
end do

!  Initialize B
do i=2, (meshN - 3)
   mat_B(i) = 0
end do
mat_B(1) = F(1)
mat_B(meshN - 2) = F(meshN)
!print *,'B:', mat_B


!print *,'A Array OUTSIDE:', mat_A
!print *, 'SHAPE: ', SHAPE(mat_A)
A_inv = inv(mat_A)
AinvB = matmul(A_inv, mat_B)
!print *,'A INV B:', AinvB
do i=2,(meshN - 1)
   F(i) = AinvB(i-1)
end do

F_num = F



!=======================================================================
!            Results
!=======================================================================
!===Print Results
print *,' '
print *,' '
print *,'=============================================================='
print *,'                      Results'
print *,'=============================================================='
print *,' '

print *,'Diffusion Coef.:', D
!print *,'Source:', S
!print *,'c:', c
!print *,'A Array:', mat_A
!print *,'B Array:', mat_B
!print *,'A_inv Array:', A_inv
!print *,'AinvB Array:', AinvB
print *,'SQRT:', sqr

!print *, 'Analytical Solution:'
!print *,'F_an Array:', F_an
print *,'Numerical Solution:'
print *,'F_num Array:', F_num
print *, ' '


!=======================================================================
!            Write to CSV
!=======================================================================

OPEN(UNIT=1,FILE=outfile,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
do i=1,meshN
!   write(1,*) x(i), ',', F(i)
   write(1,*) x(i), ',', F(i), ',', F_an(i)
end do
CLOSE(UNIT=1)
print *, 'File Written Successfully'   
print *, ' '
print *, 'Program Exexuted Successfully'


Contains
!=======================================================================
!            Functions
!=======================================================================

! Returns the inverse of a matrix calculated by 
! finding the LU decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  !real, dimension(3,3), intent(in) :: A
  !real, dimension(3,3) :: Ainv

  ! work array for LAPACK
  real(dp), dimension(size(A,1)) :: work  
  ! pivot indices
  integer, dimension(size(A,1)) :: ipiv   
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being 
  ! overwritten by LAPACK
  Ainv = A
  n = size(A,1)
!  print *,'A Array INSIDE:', Ainv

!  ! DGETRF computes an LU factorization of a general
!  ! M-by-N matrix A using partial pivoting with row
!  ! interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)
!  print *, 'info value: ', info
  if (info /= 0) stop 'Matrix is numerically singular!'

!  ! DGETRI computes the inverse of a matrix using the
!  ! LU factorization computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)
  if (info /= 0) stop 'Solution of the linear system failed!'



end function inv


end program ne470_project1
