## This is a README for Project 1

To run the Project 1 modules follow these steps:
1)  Save all files to common directory

  fortran_config.txt
  ne470_project1_6.f08
  Makefile
  singleplot
  runall

2)  Make sure .f08 file is compiled by running the following command:

  make

3)  Once an executable file (.exe) is in the directory compilation is
    complete

4)  Edit fortran_config.txt, where each line corresponds to a user input:

  Geometry
  Mesh Points
  Width
  Source #2 Multiplier
  Abs Cross Section
  Tr Cross Section
  CSV filename
  any key to continue

    The fortran_config.txt file makes it easy to run the entire module
    without entering data each time.

5)  Run runall with the following command:

  ./runall

    Plot will be saved as "singleplot.png"
    CSV data will be saved as declared

6)  (Alternative Method) Run the following commands one by one:

  ./ne470_project1_6.exe
  ./singleplot

  If you use this alternative method, you will need to enter parameters for
  fortran file one by one, manually.

Additionally, some sample plots are included in the "plots" directory.  They
contains a myriad of geometries, mesh sizes, etc., in .png format
