C_COMPILER     = g++ -std=c++11
C_FLAGS        = -O2 -fopenmp -mavx -g -pg
OUTPUT         = ../executables/dusty_tails.exe
LIBS           = -lm

HEADERS = constants.h functions.h butcher.h particle.h spline.h opacities.h

OBJ = solver.o init.o maths.o forces.o errors.o microphysics.o kvalues.o particles.o ray_tracer.o opacities.o lightcurve.o

dynamics: $(OBJ)
	$(C_COMPILER) $(C_FLAGS)  $(OBJ) $(LIBS)  -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
