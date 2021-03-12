C_COMPILER     = g++ -std=c++11
C_FLAGS        = -O3 -mavx -g -pg
OUTPUT         = ./dynamics.exe
LIBS           = -lm

HEADERS = constants.h functions.h butcher.h particle.h

OBJ = solver_new_err.o RK_init_update.o maths.o forces.o errors_update.o microphysics.o kvalues.o add_rm.o ray_tracer.o

dynamics: $(OBJ)
	$(C_COMPILER) $(C_FLAGS)  $(OBJ) $(LIBS)  -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
