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

integer :: i, j, ind, counter, newton, n, init_a, info, nmod, temp
real(real64) :: del_x, D_mult, sigma_tr, sigma_a, W, test
real(real64) :: sigma_f, k, k_old, tol1, tol2, w_new, w_old, k1, k2
real(real64) :: react_new, react_old, k_orig, w_orig
real(real64) :: W_m, w_old_m, w_new_m, w_orig_m, del_x_m
real(real64) :: cross1, cross2, D_mod, sigma_a_m, sigma_f_m, sigma_tr_m

integer, ALLOCATABLE :: ipiv(:)

!=== Initialize arrays based upon user input
real(real64), ALLOCATABLE :: mat_A(:,:), sigma_mat(:,:), flux(:), S(:), Ainv(:,:)
real(real64), ALLOCATABLE :: mat_B(:), S_old(:), S_new(:), flux_old(:), work(:)
real(real64), ALLOCATABLE :: x_step(:)
character(len=200) :: outfile

external DGETRF, DGETRS

print *,"How many Multiplying Nodes? "
read *,n

print *,"How many Moderator Nodes? "
read *,nmod


!===Initialize dynamic arrays based upon received user input
init_a = n + nmod
ALLOCATE(mat_A(init_a,init_a), sigma_mat(init_a,init_a))
ALLOCATE(mat_B(init_a), flux(init_a), flux_old(init_a))
ALLOCATE(ipiv(init_a), work(init_a), Ainv(init_a, init_a))
ALLOCATE(S_old(init_a), S_new(init_a), S(init_a), x_step(init_a))
!=======================================================================
!            Literals and Initial Conditions
!=======================================================================
!Multiplier
sigma_tr = 0.0362
sigma_a = 0.1532
sigma_f = 0.157

!Moderator
sigma_a_m = 0.00808
sigma_f_m = 0
sigma_tr_m = 0.0179


!===Diffusion Coefficient Definition
D_mult = 1/(3*sigma_tr)
D_mod = 1/(3*sigma_tr_m)

!Tolerances For Loops
tol1 = 0.000001
tol2 = 0.000001

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

w_old_m = 0
w_new_m = 0
W_m = 10
w_orig_m = W_m 


k=0
k1 = 0
k2 = 0

outfile = 'test.csv'





!=======================================================================
!            Iterations
!=======================================================================
!OUTER LOOP: Find crtitical W
do while(1 .NE. 0)
del_x = W/n
del_x_m = W_m/nmod

!===Deconstruction Matrix A
!A Matrix - Multiplier
mat_A=0   
do i=1,n
   do j=1,n
      if (i==j) mat_A(i,j) = ((2*D_mult)/(del_x**2) + sigma_a)
      if (i==j+1) mat_A(i,j) = -D_mult/(del_x**2) 
      if (i==j-1) mat_A(i,j) = -D_mult/(del_x**2)

   end do
end do
mat_A(1,1) = mat_A(1,1)/2


!A Matrix - Moderator
do i=n,init_a
   do j=n, init_a
      if (i==j) mat_A(i,j) = ((2*D_mod)/(del_x_m**2) + sigma_a_m)
      if (i==j+1) mat_A(i,j) = -D_mod/(del_x_m**2) 
      if (i==j-1) mat_A(i,j) = -D_mod/(del_x_m**2)

   end do
end do

!Multiplier - Moderator Interface Term
cross1 = D_mult/(del_x**2) + D_mod/(del_x_m**2)
cross2 = sigma_a + sigma_a_m
mat_A(n,n) =  cross1 + cross2/2

!do i=1,init_a
!      print *, 'A:', mat_A(i,:)
!end do
!read *, temp

sigma_mat = 0
!===Sources and Fission Matrix
!Multiplier => fission
do i=1,n
   do j=1,n
      IF (i==j) THEN
         sigma_mat(i,j) = sigma_f
      ELSE
         sigma_mat(i,j) = 0
      END IF
   end do
end do
sigma_mat(1,1) = sigma_mat(1,1)/2


!Calculate Inverse A Matrix Outside of Loop Once (faster)
Ainv = 0
Ainv = mat_A
call DGETRF(init_a, init_a, Ainv, init_a, ipiv, info)
if (info /= 0) stop 'Matrix is numerically singular!'
call DGETRI(init_a, Ainv, init_a, ipiv, work, init_a, info)  
if (info /= 0) stop 'Solution of the linear system failed!'
!print *, 'AINV:', Ainv


!===Initial Guesses
flux = 1
k = 1
S = 1/k*(matmul(sigma_mat, flux))

test = 1

   !INNER LOOP: Find k for W
   do while(test >= tol1)
         flux_old = flux
         k_old = k
         S_old = S   
         flux = matmul(Ainv,S_old)
   
         S = (1/k_old)*(matmul(sigma_mat,flux))
         k = k_old*(sum(S)/sum(S_old))

         !print *, 'Current K:', k
         !read *,sigma_tr
!         print *, 'flux:', flux
!         read *,sigma_tr   
   
   
         test = abs((k-k_old)/k)
         counter = counter + 1
   end do

   if (k>=1-tol2 .AND. k<=1+tol2) then
      print *, 'FOUND SOLUTION:', k
      print *, ' '   
      exit
   end if
   
!   print *, '=====BREAK====='
!   print *, 'Current K:', k
!   print *, 'W:', W
!   print *, ' '
!   read *,sigma_tr


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
!            Write to CSV
!=======================================================================

do i=1,n
   x_step(i) = (i-1)*W/(n-1)
end do

do i=1,nmod
   x_step(i+n) = (i)*W_m/(nmod-1) + x_step(n)
end do

OPEN(UNIT=1,FILE=outfile,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
do i=1,init_a
!   write(1,*) x(i), ',', F(i)
   write(1,*) x_step(i), ',', flux(i)
end do
CLOSE(UNIT=1)
print *, 'File Written Successfully'   
print *, ' '
print *, 'Program Exexuted Successfully'



!=======================================================================
!            Results
!=======================================================================
!===Print Results
print *,' '
print *,' '
print *,'=============================================================='
print *,'             1Group w/ Moderator Final Results'
print *,'=============================================================='
print *,' '

print *,'Diffusion Coef.:', D_mult
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
write(*,*) 'Moderator Width [cm] = ', W_m
write(*,*) 'Half Critical Width [cm] = ', W
print *, ' '

   
   


end program ne470_project2_7
