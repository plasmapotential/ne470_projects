EXECUTABLE=ne470_project2_7.exe
SOURCES=ne470_project2_7.f08

#EXECUTABLE=tester_lapack.exe
#SOURCES=tester_lapack.f08

OBJECTS=$(SOURCES:.f08=.o)
LDFLAGS=-llapack -lblas

FC=gfortran
#CFLAGS=-g -fcheck=all -ftrapv -ffpe-trap=overflow,invalid -Wall -Wextra
CFLAGS=-O3 -Wall -Wextra

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@

%.o: %.f08
	$(FC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(EXECUTABLE)
