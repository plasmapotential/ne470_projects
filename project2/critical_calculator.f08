!ne470_projecttester.f08

!Title:			NE 470 Project Test Script
!Engineer:		Tom Looby
!Date:			10/28/2017
!Description:	Solves for k

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program tester
use iso_fortran_env
implicit none

!=======================================================================
!            Variables
!=======================================================================

integer :: i, j, ind, counter, newton, n, init_a, info
real(real64) :: del_x, D, sigma_tr, sigma_a, W, test
real(real64) :: sigma_f, k, k_old, tol1, tol2, w_new, w_old, k1, k2
real(real64) :: react_new, react_old, flx_dif

integer, ALLOCATABLE :: ipiv(:)

!=== Initialize arrays based upon user input
real(real64), ALLOCATABLE :: mat_A(:,:), sigma_mat(:,:), flux(:), S(:), Ainv(:,:)
real(real64), ALLOCATABLE :: mat_B(:), S_old(:), S_new(:), flux_old(:), work(:)
!real(real64), ALLOCATABLE :: ipiv(:)

external DGETRF, DGETRS

print *,"Size of the problem n? "
read *,n



!===Initialize dynamic arrays based upon received user input
init_a = n
ALLOCATE(mat_A(init_a,init_a), sigma_mat(init_a,init_a))
ALLOCATE(mat_B(init_a), flux(init_a), flux_old(init_a))
ALLOCATE(ipiv(init_a), work(init_a), Ainv(init_a, init_a))
ALLOCATE(S_old(init_a), S_new(init_a), S(init_a))
!=======================================================================
!            Literals and Initial Conditions
!=======================================================================
!Testing Stuff
sigma_tr = 0.0362
sigma_a = 0.1532
sigma_f = 0.1570


!===Diffusion Coefficient Definition
D = 1/(3*sigma_tr)

!Tolerances For Loops
tol1 = 0.00001
tol2 = 0.00001

!===Counters and the like
counter = 0
ind = 1
newton = 1
j=1
i=1


!Guess Analytical Critical Width
!W = (3.1415)*sqrt(D/(sigma_f-sigma_a))
W=100

del_x = W/n


!===Sources and Fission Matrix
!S = (10**8)

do i=1,n
   do j=1,n
      IF (i==j) THEN
         sigma_mat(i,j) = sigma_f
      ELSE
         sigma_mat(i,j) = 0
      END IF
   end do
end do
sigma_mat(1,1) = sigma_f/2

!Deconstruction Matrix A
mat_A=0   
do i=1,n
   do j=1,n
      if (i==j) mat_A(i,j) = ((2*D)/(del_x**2) + sigma_a)
      if (i==j+1) mat_A(i,j) = -D/(del_x**2) 
      if (i==j-1) mat_A(i,j) = -D/(del_x**2)

   end do
end do
mat_A(1,1) = mat_A(1,1)/2

!Initial Guesses
flux = 1
k = 1
S = 1/k*(matmul(sigma_mat,flux))
test = 1
flx_dif = 1


Ainv = 0
Ainv = mat_A
call DGETRF(n, n, Ainv, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'
call DGETRI(n, Ainv, n, ipiv, work, n, info)  
if (info /= 0) stop 'Solution of the linear system failed!'
!print *, 'AINV:', Ainv

!=======================================================================
!            Iterations
!=======================================================================

do while(test >= tol1)
      flux_old = flux
      k_old = k
      S_old = S   
      flux = matmul(Ainv,S_old)

      S = (1/k_old)*(matmul(sigma_mat,flux))
      k = k_old*(sum(S)/sum(S_old))
!      print *, 'Current K:', k

      test = abs((k-k_old)/k)
      counter = counter + 1
end do

print *, 'D', D
print *, 'TEST', test
print *, 'Current K:', k
print *, 'W', W
print *, 'delx', del_x
end program tester
  
   
   



