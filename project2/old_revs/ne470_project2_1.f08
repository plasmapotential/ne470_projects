!ne470_project2.f08

!Title:			NE 470 Project #2
!Engineer:		Tom Looby
!Date:			10/28/2017
!Description:	Solves for critical reactor size

!=======================================================================
!            MAIN PROGRAM
!=======================================================================
program ne470_project2
implicit none

!=======================================================================
!            Variables, Literals
!=======================================================================

integer :: i, j, ierror, meshN, geom, init_a, n, info, ind, newton
integer, parameter :: dp = selected_real_kind(15, 307)
real(dp) :: del_x, S, S_mult, D, sigma_tr, sigma_a, c, S2, sqr, W
real(dp) :: sigma_f, k, k_old, tol1, tol2, w_new, w_old, w_dot, k1, k2
real(dp) :: react_new, react_old, counter, d_ro

!=== Initialize arrays based upon user input
real(dp), ALLOCATABLE :: mat_A(:,:), A_inv(:,:)
real(dp), ALLOCATABLE :: mat_B(:), AinvB(:), S_old(:), S_new(:)
real(dp), ALLOCATABLE :: F(:), F_an(:), F_num(:), x(:), work(:)


integer, ALLOCATABLE :: ipiv(:)
logical :: exists
character(len=200) :: outfile, dummy

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

!===User defined Width with error checking===
do
  write(*,*) "Enter Width Guess (cm):"

  read(*,'(F10.0)',iostat=ierror) W
  if ( ierror == 0 .AND. W >= 0)  then
    exit
  endif
  write(*,*) 'Please enter an numeric greater than 0.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Width Guess: ', W
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

do
  write(*,*) "Enter Macroscopic Fission Cross Section (cm^-1):"

  read(*,'(F10.0)',iostat=ierror) sigma_f
  if ( ierror == 0 )  then
    exit
  endif
  write(*,*) 'Please enter a numeric.  Try Again...'
  write(*,*) ' '

enddo
write(*,*) 'Macroscopic Fission Cross Section: ', sigma_f
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
ALLOCATE(work(init_a), ipiv(init_a))
ALLOCATE(S_old(init_a), S_new(init_a))
!=======================================================================
!            Literals and Initial Conditions
!=======================================================================
!===Create c value from geometry selection
c = geom - 1
S = (10**8)
S2 = (10**8)*S_mult
!===Diffusion Coefficient Definition
D = 1/(3*sigma_tr)

!Tolerances For Loops
tol1 = 0.01
tol2 = 0.01


counter = 0
k = 1
k1 = 1
k_old = 0
ind = 1
w_old = W


!=======================================================================
!           Numerical Solution             
!=======================================================================

do while (1 .NE. 0)
   do newton = 1, 2 
      
      !===Step Size
      do i=1,(meshN)
         x(i) = (i-1)*W/(meshN-1)
      end do
      del_x = x(2) - x(1)
      
      !===Initial Conditions
      F(1) = S+(-(F(1)/(2*D*del_x)))
      F(meshN) = S2+((-F(meshN)/(2*D*del_x)))
      
      
      !=====Destruction Operator (A Matrix)
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
      write(*,*) 'A: ', mat_A
      write(*,*) ' '
      !==========Initial Guesses
      
        
      !Source Matrix
      mat_B(1) = S*(((D)/(del_x**2)) + sigma_a)*1/(c+1)
      mat_B(meshN - 2) =  S2*(((D)/(del_x**2)) + sigma_a + (D/(del_x**2))*(c/(2*meshN - 1)))
      do j=2, (meshN - 3)
            mat_b(j) = sigma_f
      end do
      S_old = mat_B
      k = 1
      k_old = 0
      write(*,*) 'B: ', mat_B
   
      
      
      !+++Inner Loop.  Finding k
         do while (1 .NE. 0)
            !At beginning of loop, mat_B = S(n)
            !After LAPACK, mat_B = Flux
            !Then mat_B is multiplied by sigma_f and k^-1
            !This Yields S(n+1)
            !Doing it this way conserves memory
            
            !print *,'B:', mat_B
           
            !====LAPACK STUFF
            n = size(mat_A,1)
            ! DGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            ! See: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
            call DGETRF(n, n, mat_A, n, ipiv, info)
            if (info /= 0) stop 'Matrix is numerically singular!'
            
            ! DGETRS - solve a system of linear equations A * X = B or 
            ! A' * X = B with a general N-by-N matrix A using the LU 
            ! factorization computed by DGETRF. 
            ! See: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
            call DGETRS('N', n, 1, mat_A, n, ipiv, mat_B, n, info)  
            if (info /= 0) stop 'Solution of the linear system failed!'
         
            !mat_B = Flux here, so fill up the flux matrix
   !         do i=2,(meshN - 1)
   !            F(i) = mat_B(i-1)
   !         end do
         
            !now we convert mat_B to S(n+1)
            S_new = (mat_B*sigma_f)
            !write(*,*) 'S_NEW: ', S_new
            !write(*,*) ' '      
            !Multiplication Factor / Reactivity
            k_old = k
            k = (k_old*sum(S_new))/(sum(S_old))
            
            if (abs((k - k_old)/k) <= tol1) then  
            !if ((k - k_old)/k >= tol1 .AND. (S_new(meshN - 3) - (S_old(meshN-3))/S_new(meshN-2)) >= tol1) then
            !   write(*,*) 'FOUND A K!', k
               exit
            end if
            
            
            !Save this S value
            S_old = S_new
         
            !Create New RHS Source Term: S(n+1)/k(n+1)
            mat_B = S_new/k
         
            !counter = counter + 1
            !write(*,*) 'INSIDE k = ', k
         end do   
      !write(*,*) 'OUT OF INNER LOOP!'
 
      
   
      if (newton == 1) then
         w_old = W
         k1 = k
         W = W + 0.1
      else if (newton == 2) then
         w_new = W
         k2 = k
      end if
      
   end do
   
   
   IF (k1 >= (1-tol2) .AND. k1 <= (1+tol2)) THEN
      write(*,*) 'FOUND SOLUTION!'
      write(*,*) 'ind = ', ind
      !write(*,*) 'W = ', W
      !write(*,*) 'k = ', k
      write(*,*) 'Counter = ', counter
      exit
   END IF   

   
   
   react_new = (k2-1)/k2
   react_old = (k1-1)/k1
   d_ro = (react_new-react_old)/(w_new - w_old)
   W = w_new - (react_new/d_ro)
!   react_new = k2-1
!   react_old = k1-1
!   d_ro = (react_new-react_old)/(w_new - w_old)
!   W = w_new - (react_new/d_ro)
   
   write(*,*) 'K1 = ', k1
   write(*,*) 'K2 = ', k2
   write(*,*) 'react_old = ', react_old 
   write(*,*) 'react_new = ', react_new
   write(*,*) 'd_ro = !', d_ro
   write(*,*) 'w_old = ', w_old
   write(*,*) 'w_new = ', w_new


   write(*,*) 'FOUND A K!', k1
   write(*,*) 'W = ', W   
   write(*,*) 'del_x = ', del_x 
   write(*,*) ' '
  
   ind = ind + 1
         
end do
   

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
write(*,*) 'Multiplication Factor (k) = ', k
write(*,*) 'Critical Width [cm] = ', W
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
