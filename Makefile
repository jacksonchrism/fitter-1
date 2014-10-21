CC= g++

LDFLAGS=

#CFLAGS=-ggdb

SOURCES=wwfitter.cpp MyMiniSim.cc 

#INCLUDES=-I. -I/home/kate -I${RATROOT}/include/RAT -I${RATROOT}/src/geo -I${RATROOT}/src/stlplus -I${RATROOT}/src/physics/PhotonThinning.cc -I${RATROOT}/include -I${RATROOT}/src/core -I${RATROOT}/src/util -I${RATROOT}/src/physics -L${RATROOT}/lib `geant4-config --cflags` `geant4-config --libs` -lRATEvent_Linux-g++ -lrat_Linux-g++ `clhep-config --libs` `geant4-config --cflags` `geant4-config --libs` -lMinuit `root-config --libs --cflags`
INCLUDES=-I. -I/home/kate -I${RATROOT}/include/RAT -I${RATROOT}/src/geo -I${RATROOT}/src/stlplus -I${RATROOT}/include -I${RATROOT}/src/core -I${RATROOT}/src/util -I${RATROOT}/src/physics -L${RATROOT}/lib `geant4-config --cflags` `geant4-config --libs` -lRATEvent_Linux -lrat_Linux `geant4-config --cflags` `geant4-config --libs` -lMinuit `root-config --libs --cflags`

OBJECTS=$(SOURCES: .c*=.o)

EXECUTABLE=constaging_test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) $(INCLUDES) -o $@

.cpp.o:
	$(CC) $< -o $@
