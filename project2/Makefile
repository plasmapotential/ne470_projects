EXECUTABLE=ne470_project2_5.exe
SOURCES=ne470_project2_5.f08

OBJECTS=$(SOURCES:.f90=.o)
LDFLAGS=-llapack

FC=gfortran
#CFLAGS=-g -fcheck=all -ftrapv -ffpe-trap=overflow,invalid -Wall -Wextra
CFLAGS=-O3 -Wall -Wextra

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(EXECUTABLE)
