program chvala_lu
! LAPACK linear system solver example, Ondrej Chvala <ochvala@utk.edu> 2014-10-31
!
! Uses LAPACK's DGETRF to compute an LU factorization of a matrix A
! using partial pivoting with row interchanges, and then DGETRS to
! solve a system of linear equations A*x=b by the LU factorization
! already computed by DGETRF.
!
use iso_fortran_env
implicit none

real(real64), allocatable :: A(:,:)  ! matrix to LU decompose
real(real64), allocatable :: b(:)    ! b vector
real(real64), allocatable :: ipiv(:) ! pivoting array for LAPACK
integer                   :: n       ! problem size
integer                   :: info    ! success flag for LAPACK
character(len=200)        :: fmt     ! formatting string for pretty 
printing
integer                   :: i, j, ioerr
character(*), parameter :: outfile = 'results.txt'
real(real64), parameter :: diagA =   2.0 ! diagonal element
real(real64), parameter :: nexdA =  -1.0 ! next to diagonal element
real(real64), parameter :: S1    =  10.0 ! "source"
real(real64), parameter :: Width = 100.0 ! width [cm]
real(real64) :: x
! External procedures defined in LAPACK
external DGETRF, DGETRS
! Input data
print *,"Size of the problem n? "
read *,n
allocate (A(n,n),b(n),ipiv(n))
!print *,"Elements of the A matrix ?"
!read *,A
!print *,"Elements of the b vector ?"
!read *,b
! Allocate matrix A and vector b
A = 0.0
b = 0.0
b(1) = S1
do i = 1, n
do j = 1, n
if(i .eq. j) A(i,j) = diagA
if((i.eq. j+1).or.(i.eq. j-1))  A(i,j) = nexdA
enddo
enddo
A(1,1) = diagA/2.0

call DGETRF(n, n, A, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'

call DGETRS('N', n, 1, A, n, ipiv, b, n, info)
if (info /= 0) stop 'Solution of the linear system failed!'
b = b/real(N,real64) ! normalization


! Print solution
!!write (fmt,"(a,i3,a)") "(",n,"(f10.4,1x))"
!print fmt, b
!! Write results to a file
!open(10, file=outfile, status='replace', action='write', iostat=ioerr)
!if(ioerr.ne.0) stop 'File open error'
!write(10,*) "# i   x   flux(x)   cos(x)"
!do i = 1, n
!x = Width * real(i-1,real64) / real(N-1,real64)
!write (10,999) i, x, b(i), S1*cos(2.0_real64 * atan(1.0) * real(i-
!1,real64) / real(N-1,real64) )
!enddo
!close(10)
!999 format(i3, f10.3, f10.3, f10.4)
end program chvala_lu
