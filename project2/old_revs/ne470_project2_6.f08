!ne470_project2.f08

!Title:			NE 470 Project #2
!Engineer:		Tom Looby
!Date:			10/28/2017
!Description:	Solves for critical reactor size

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program ne470_project2
use iso_fortran_env
implicit none

!=======================================================================
!            Variables, Literals
!=======================================================================

integer :: i, j, ind, counter, newton, n, init_a, info
real(real64) :: del_x, D, sigma_tr, sigma_a, W, test
real(real64) :: sigma_f, k, k_old, tol1, tol2, w_new, w_old, k1, k2
real(real64) :: react_new, react_old, flx_dif

!integer, ALLOCATABLE :: ipiv(:)

!=== Initialize arrays based upon user input
real(real64), ALLOCATABLE :: mat_A(:,:), sigma_mat(:,:), flux(:), S(:) 
real(real64), ALLOCATABLE :: mat_B(:), S_old(:), S_new(:), flux_old(:)
real(real64), ALLOCATABLE :: ipiv(:)

external DGETRF, DGETRS

!=======================================================================
!           Print Stuff / Instructions
!=======================================================================

print *, '======================Project 1 ============================='
print *, 'Engineers: Group 1'
print *, 'Date: 10/07/2017'
print *, ' '
print *, 'This program calculates critical reactor size, given various '
print *, 'input parameters '
print *, '=============================================================' 
print *, ' '

print *,"How Many Nodes? "
read *,n


!===Initialize dynamic arrays based upon received user input

init_a = n
ALLOCATE(mat_A(init_a,init_a), sigma_mat(init_a,init_a))
ALLOCATE(mat_B(init_a), flux(init_a), flux_old(init_a))
ALLOCATE(ipiv(init_a))
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
test = 1
!test = .064855739573128068

!===Counters and the like
counter = 0
ind = 1
newton = 1
j=1
i=1


!===Sources
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

W = 100
react_new = 0
react_old = 0





!=======================================================================
!           Numerical Solution             
!=======================================================================
do while (1 .NE. 0)
   del_x = W/(n-1)
   flux = 1
   !flux(1) = sigma_f/2
   k = 1
   k_old = 1
   S = (matmul(sigma_mat,flux))
   mat_A=0
   !flux(1) = 0
   !flux(n) = 0
      
   do i=1,n
      do j=1,n
         if (i==j) mat_A(i,j) = ((2*D)/(del_x**2) + sigma_a)
         if (i==j+1) mat_A(i,j) = -D/(del_x**2) 
         if (i==j-1) mat_A(i,j) = -D/(del_x**2)

      end do
   end do
   mat_A(1,1) = mat_A(1,1)/2
   !print *, 'A:', mat_A
   !read *,n
   
   do while(1 .NE. 0)
      k_old = k
      S_old = S   
      mat_B = (1/k)*S
      flux_old = flux

      !mat_B(n) = 0

     
      call DGETRF(n, n, mat_A, n, ipiv, info)
      if (info /= 0) stop 'Matrix is numerically singular!'
      
      call DGETRS('N', n, 1, mat_A, n, ipiv, mat_B, n, info)  
      if (info /= 0) stop 'Solution of the linear system failed!'

      flux = mat_B*n
      S = (matmul(sigma_mat,flux))

      k = k_old*(sum(S)/sum(S_old))
      !print *, 'Current K:', k
  
      
      if (abs((k-k_old)/k) <= tol1) then 
         !print *, 'BREAK'
         exit
      end if
      
      counter = counter + 1
   end do
   
   
   if (k>=test-tol2 .AND. k<=test+tol2) then
      print *, 'FOUND SOLUTION:', k
      print *, ' '   
      exit
   end if
   
   if (newton == 1) then
      w_new = W
      k2 = k
      W = w_new+0.5
      newton = 2
   
   else
      w_old = w_new
      w_new = W
      k1 = k2
      k2 = k
      !react_new = log(k2)
      !react_old = log(k1)
      react_new = (k2-1)/k2
      react_old = (k1-1)/k1
      !W = (w_old*react_new - w_new*react_old)/(react_new - react_old)
      W = w_new - (react_new*(w_new - w_old))/(react_new - react_old)
      !read *,n
   end if
   
!   print *, 'Upcoming W:', W
!   print *, 'Last W:', w_new
!   print *, 'Current React:', react_new
!   print *, 'Last React:', react_old
!   print *, 'K2:', k2
!   print *, 'K1:', k1
!   print *, 'Current K:', k
!   print *, ' '
   !read *,n

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
write(*,*) 'Final Multiplication Factor (k) = ', k
write(*,*) 'Critical Width [m] = ', W
print *, ' '


!=======================================================================
!            Write to CSV
!=======================================================================

!OPEN(UNIT=1,FILE=outfile,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!do i=1,meshN
!!   write(1,*) x(i), ',', F(i)
!   write(1,*) x(i), ',', F(i), ',', F_an(i)
!end do
!CLOSE(UNIT=1)
!print *, 'File Written Successfully'   
!print *, ' '
!print *, 'Program Exexuted Successfully'

end program ne470_project2
