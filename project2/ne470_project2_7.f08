!ne470_project2_7.f08

!Title:			NE 470 Project Test Script
!Engineer:		Tom Looby
!Date:			10/28/2017
!Description:	Solves for k

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program ne470_project2_7
use iso_fortran_env
implicit none

!=======================================================================
!            Variables
!=======================================================================

integer :: i, j, ind, counter, newton, n, init_a, info
real(real64) :: del_x, D, sigma_tr, sigma_a, W, test
real(real64) :: sigma_f, k, k_old, tol1, tol2, w_new, w_old, k1, k2
real(real64) :: react_new, react_old, k_orig, w_orig

integer, ALLOCATABLE :: ipiv(:)

!=== Initialize arrays based upon user input
real(real64), ALLOCATABLE :: mat_A(:,:), sigma_mat(:,:), flux(:), S(:), Ainv(:,:)
real(real64), ALLOCATABLE :: mat_B(:), S_old(:), S_new(:), flux_old(:), work(:)
!real(real64), ALLOCATABLE :: ipiv(:)

external DGETRF, DGETRS

print *,"How many mesh Nodes? "
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
w_old = 0
w_new = 0
W=100
w_orig = W
k=0
k1 = 0
k2 = 0

!===Sources and Fission Matrix
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

!=======================================================================
!            Iterations
!=======================================================================
!OUTER LOOP: Find crtitical W
do while(1 .NE. 0)
del_x = W/n

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

!Calculate Inverse A Matrix Outside of Loop Once (faster)
Ainv = 0
Ainv = mat_A
call DGETRF(n, n, Ainv, n, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'
call DGETRI(n, Ainv, n, ipiv, work, n, info)  
if (info /= 0) stop 'Solution of the linear system failed!'
!print *, 'AINV:', Ainv


!Initial Guesses
flux = 1
k = 1
S = 1/k*(matmul(sigma_mat,flux))
test = 1

   !INNER LOOP: Find k for W
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

   if (k>=1-tol2 .AND. k<=1+tol2) then
      print *, 'FOUND SOLUTION:', k
      print *, ' '   
      exit
   end if



   if (newton == 1) then
      w_new = W
      k2 = k
      k_orig = k
      W = w_new*1.01
      newton = 2
   
   else
      w_old = w_new
      w_new = W
      k1 = k2
      k2 = k
      react_new = (k2-1)/k2
      react_old = (k1-1)/k1
      W = (w_old*react_new - w_new*react_old)/(react_new - react_old)
   end if
   
   ind = ind+1
end do

!=======================================================================
!            Results
!=======================================================================
!===Print Results
print *,' '
print *,' '
print *,'=============================================================='
print *,'             Project 2 Final Results'
print *,'=============================================================='
print *,' '

print *,'Diffusion Coef.:', D
print *, ' '
write(*,*) 'Inner Loop Tolerance = ', tol1
write(*,*) 'Outer Loop Tolerance = ', tol2
print *, ' '
write(*,*) 'Outer Loop Iterations = ', ind
write(*,*) 'Inner Loop Iterations = ', counter
print *, ' '
write(*,*) 'Original Width Guess [cm] = ', w_orig
write(*,*) 'Original Width k Value= ', k_orig
print *, ' '
write(*,*) 'Final Multiplication Factor (k) = ', k
write(*,*) 'Critical Width [cm] = ', W*2
print *, ' '

   
   


end program ne470_project2_7
