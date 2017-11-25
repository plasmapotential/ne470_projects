!ne470_project3_1.f08

!Title:			NE 470 Project 3 Tester
!Engineer:		Tom Looby
!Date:			11/12/2017
!Description:	Solves Multigroup Diffusion Equations w/ no moderator

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program ne470_project3_3
use iso_fortran_env
implicit none


!=======================================================================
!            Variables
!=======================================================================

integer :: i, j, g, ind, counter, newton, n, init_a, info, gN, nmod
real(real64) :: del_x, W, test, temp, cross1, cross2, del_x_m
real(real64) :: k, k_old, tol1, tol2, w_new, w_old, k1, k2, Wmod
real(real64) :: react_new, react_old, k_orig, w_orig, ratio

integer, ALLOCATABLE :: ipiv(:)

!=== Initialize arrays based upon user input
real(real64), ALLOCATABLE :: mat_A(:,:,:), sigma_mat(:,:,:), flux(:,:)
real(real64), ALLOCATABLE :: S(:), Ainv(:,:,:), work(:), Atemp(:,:)
real(real64), ALLOCATABLE :: mat_B(:), S_old(:), S_new(:), flux_old(:,:)
!real(real64), ALLOCATABLE :: ipiv(:)

real(real64), ALLOCATABLE :: vsigf(:), sigf(:), siga(:), D(:), sigR(:) 
real(real64), ALLOCATABLE :: D_h2o(:), sigtr_h2o(:), siga_h2o(:)
real(real64), ALLOCATABLE :: sigs_h2o(:,:), sigR_h2o(:), RHS(:), X(:)
real(real64), ALLOCATABLE :: sigs(:,:), scat(:,:), x_step(:)
real(real64), ALLOCATABLE :: fis_mat(:,:), s_mat(:,:)

character(len=200) :: outfile


external DGETRF, DGETRS

print *,"How many Multiplier Nodes? "
read *,n

nmod = 0
print *,"How many Groups? "
read *,gN

!=======================================================================
!            Literals and Initial Conditions
!=======================================================================
!Allocate constants for number of groups
if (gN /= 1 .AND. gN /= 2 .AND. gN /= 4) then
   gN=1
end if

ALLOCATE(vsigf(gN), sigf(gN), siga(gN), D(gN), sigR(gN), D_h2o(gN))
ALLOCATE(sigtr_h2o(gN), siga_h2o(gN), sigs_h2o(gN,gN), sigR_h2o(gN))



!===Initialize dynamic arrays based upon received user input
init_a = n + nmod
ALLOCATE(mat_A(init_a,init_a, gN), sigma_mat(init_a,init_a, gN))
ALLOCATE(mat_B(init_a), flux(init_a, gN), flux_old(init_a,gN))
ALLOCATE(ipiv(init_a), work(init_a), Ainv(init_a, init_a, gN))
ALLOCATE(S_old(init_a), S_new(init_a), S(init_a), Atemp(init_a,init_a))
ALLOCATE(RHS(init_a), X(gN), sigs(gN, gN), scat(init_a,gN), x_step(init_a))
ALLOCATE(fis_mat(init_a, init_a), s_mat(init_a,gN))

if (gN == 1) then
   !Multiplier
   sigR(1) = 0.1532
   siga(1) = 0.1532
   vsigf(1) = 0.157
   
   !Moderator
   siga_h2o(1) = 0.00808
   sigR_h2o(1) = 0.00808
   
   
   !===Diffusion Coefficient Definition
   D = 1/(3*0.0362)
   D_h2o = 1/(3*0.0179)


elseif (gN == 4) then
   !+++Four Group Constants
   vsigf(1) = 0.0009572
   vsigf(2) = 0.001193
   vsigf(3) = 0.01768
   vsigf(4) = 0.18514
   
   sigf(1) = 0.003378
   sigf(2) = 0.0004850
   sigf(3) = 0.006970
   sigf(4) = 0.07527
   
   siga(1) = 0.004946
   siga(2) = 0.002840
   siga(3) = 0.03053
   siga(4) = 0.1210
   
   D(1) = 2.1623
   D(2) = 1.0867
   D(3) = 0.6318
   D(4) = 0.3543
   
   sigR(1) = 0.08795
   sigR(2) = 0.06124
   sigR(3) = 0.09506
   sigR(4) = 0.1210
   
   sigs = 0
   sigs(1,2) = sigR(1)
   sigs(2,3) = sigR(2)
   sigs(3,4) = sigR(3)
   
   
   !H20 cross sections
   sigtr_h2o(1) = 0.20608
   sigtr_h2o(2) =0.60215
   sigtr_h2o(3) =0.56830
   sigtr_h2o(4) =1.2111
   
   siga_h2o(1) = 0.00051
   siga_h2o(2) = 0.00354
   siga_h2o(3) = 0.01581
   siga_h2o(4) = 0.04637
   
   sigs_h2o = 0
   !Original Cross Sections
!   sigs_h2o(1,1) = 0.37045
!   sigs_h2o(1,2) = 0.04152
!   sigs_h2o(1,3) = 0.00001
!   sigs_h2o(2,2) = 0.98285
!   sigs_h2o(2,3) = 0.07459
!   sigs_h2o(2,4) = 0.01371
!   sigs_h2o(3,3) = 0.76110
!   sigs_h2o(3,4) = 0.31856
!   sigs_h2o(4,3) = 0.00085
!   sigs_h2o(4,4) = 1.96607

   !Assuming Direct Coupling
   sigs_h2o(1,1) = 0
   sigs_h2o(1,2) = 0.04152
   sigs_h2o(1,3) = 0
   sigs_h2o(2,2) = 0
   sigs_h2o(2,3) = 0.07459
   sigs_h2o(2,4) = 0
   sigs_h2o(3,3) = 0
   sigs_h2o(3,4) = 0.31856
   sigs_h2o(4,3) = 0
   sigs_h2o(4,4) = 0
   
   do i=1,gN
      D_h2o(i) = 1/(3*sigtr_h2o(i))
   end do
   
  


elseif (gN == 2) then
   !+++Two Group Constants
   vsigf(1) = 0.008476
   vsigf(2) = 0.18514
   
   sigf(1) = 0.003320
   sigf(2) = 0.07537
   
   siga(1) = 0.01207
   siga(2) = 0.1210
   
   D(1) = 1.2627
   D(2) = 0.3543
   
   sigR(1) = 0.02619
   sigR(2) = 0.1210

   
   !sigs(1,1) = sigR(1) - siga(1)
   sigs(1,1) = 0
   sigs(1,2) = sigR(1)
   sigs(2,1) = 0
   sigs(2,2) = 0
   
  
   !H20 Cross Sections 
   D_h2o(1) = 1.13
   D_h2o(2) = 0.16
   
   sigtr_h2o(1) = 0.20608
   sigtr_h2o(2) =0.60215

   siga_h2o(1) = 0.0004
   siga_h2o(2) = 0.0197   
   
   sigR_h2o(1) = 0.0494
   sigR_h2o(2) = siga_h2o(2)
   
   sigs_h2o = 0
   sigs_h2o(1,2) = sigR_h2o(1) - siga_h2o(1) 

   
end if

!Tolerances For Loops
tol1 = 0.000001
tol2 = 0.000001

!===Counters and the like
counter = 0
ind = 1
newton = 1
j=1
i=1

!Width Stuff
w_old = 0
w_new = 0
W=100
w_orig = W
Wmod = 20

!K_effective
k=0
k1 = 0
k2 = 0
!Patition Function - Fission Neutrons only born in Fast Group
X=0
X(1)=1


outfile = 'test.csv'

!=======================================================================
!            Solution
!=======================================================================

!OUTER LOOP: Find crtitical W
do while(1 .NE. 0)
del_x = W/n
del_x_m = Wmod/nmod
ratio=del_x_m/del_x
sigma_mat=0

   !====Deconstruction Matrix A (3-D) with Removal and Moderator
   !Multiplier===
   mat_A=0   
   do g=1,gN
      do i=1,n
         do j=1,n
            if (i==j .AND. g/=gN) then
               mat_A(i,j,g) = ((2*D(g))/(del_x**2) + sigR(g))
            else if (i==j .AND. g==gN)  then
               mat_A(i,j,g) = ((2*D(g))/(del_x**2) + siga(g))
            end if
            if (i==j+1) mat_A(i,j,g) = -D(g)/(del_x**2) 
            if (i==j-1) mat_A(i,j,g) = -D(g)/(del_x**2)
      
         end do
      end do

      mat_A(1,1,g) = mat_A(1,1,g)/2
      
!      do i=1,init_a
!         print *, 'A:', mat_A(i,:,g)
!      end do
!      read *, temp
    
  
      !Calculate Inverse A Matrix (3-D)
      Atemp = 0
      Atemp = mat_A(:,:,g)
      call DGETRF(init_a, init_a, Atemp, init_a, ipiv, info)
      if (info /= 0) stop 'Matrix is numerically singular!'
      call DGETRI(init_a, Atemp, init_a, ipiv, work, init_a, info)  
      if (info /= 0) stop 'Solution of the linear system failed!'
!      print *, 'Atemp:', Atemp
      Ainv(:,:,g) = Atemp

      !Fission Matrix
      do i=1,n
         do j=1,n
            IF (i==j) THEN
               sigma_mat(i,j,g) = vsigf(g)
            ELSE
            sigma_mat(i,j,g) = 0
               
            END IF
         end do
      end do
      sigma_mat(1,1,g) = sigma_mat(1,1,g)/2
    end do


   !Initial Guesses
   flux = 1
   k = 1

   !Build Source Matrix
   do g=1,gN
      s_mat(:,g) = matmul(sigma_mat(:,:,g),flux(:,g))
   end do

   S=sum(s_mat, DIM = 2)
   RHS = 1/k*X(1)*S
   test = 1
   scat=0
   
   !INNER LOOP: Find k for W
   do while(test >= tol1)
      k_old = k
      S_old = S
      scat=0
      
      if (gN/=1) then
         do g=1,gN-1
            flux(:,g) = matmul(Ainv(:,:,g), RHS)

            !Build Scattering Matrix
            scat(1:n,g+1) = (flux(1:n,g)*sigs(g,g+1))
            !scat(n:init_a,g+1) = (flux(n:init_a,g)*sigs_h2o(g,g+1))
            
            RHS = 1/k*X(g+1)*S + scat(:,g+1)
            RHS = scat(:,g+1)

         end do
      end if
      flux(:,gN) = matmul(Ainv(:,:,gN), RHS)
      
      !New Fission Matrix
      do g=1,gN
         s_mat(:,g) = matmul(sigma_mat(:,:,g), flux(:,g))
      end do
      S = sum(s_mat, DIM=2)
      
      k = k_old*sum(S)/sum(S_old)
      RHS = 1/k*X(1)*S


      test = abs((k-k_old)/k)
      counter = counter + 1
!      print *, 'Current S:', S
!      print *, ' '
!      print *, 'Current K:', k
!      print *, ' '
!      print *, 'Fast Flux', flux(:,1)
!      print *, ' '
!      print *, 'Thermal Flux', flux(:,2)
!      print *, ' '
!      read *, temp
   end do


   if (k>=1-tol2 .AND. k<=1+tol2) then
      print *, 'FOUND SOLUTION:', k
      print *, ' '   
      exit
   end if

   print *, '====================BREAK======================'
   print *, 'k = ', k
   print *, 'Fast Flux = ', flux(:,1)
   print *, 'Thermal flux = ', flux(:,2)
   print *, 'Width = ', W
   print *, ' '
   print *, ' '
   read *, temp
   




   if (newton == 1) then
      w_new = W
      k2 = k
      k_orig = k
      W = w_new*1.1
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
   x_step(i+n) = (i)*Wmod/(nmod-1) + x_step(n)
end do

OPEN(UNIT=1,FILE=outfile,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
do i=1,init_a
!   write(1,*) x(i), ',', F(i)
   write(1,*) x_step(i), ',', flux(i, 1), ',', flux(i, 2)!, ',', flux(i, 3), ',', flux(i, 4)
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
print *,'             Project 2 Final Results'
print *,'=============================================================='
print *,' '

print *,'Groups:', gN
print *, ' '
print *,'Diffusion Coef.s:', D
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
write(*,*) 'Moderator Width [cm] = ', Wmod
write(*,*) 'Half Critical Width [cm] = ', W
print *, ' '
   




end program ne470_project3_3
