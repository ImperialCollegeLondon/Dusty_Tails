C_COMPILER     = g++ -std=c++11
C_FLAGS        = -O3 -mavx -g 
OUTPUT         = ./dynamics_040_kic1255b_sph_tau_2orb.exe
LIBS           = -lm

HEADERS = constants.h functions.h butcher.h particle.h spline.h

OBJ = solver.o init.o maths.o forces.o errors.o microphysics.o kvalues.o particles.o ray_tracer.o tau_spline.cpp

dynamics: $(OBJ)
	$(C_COMPILER) $(C_FLAGS)  $(OBJ) $(LIBS)  -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
