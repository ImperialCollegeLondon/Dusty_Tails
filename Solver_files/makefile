C_COMPILER     = g++
C_FLAGS        = -O3
OUTPUT         = ../dynamics.exe
LIBS           = -lm

HEADERS = constants.h functions.h butcher.h particle.h

OBJ = solver_update.o RK_init_update.o maths.o forces.o errors.o microphysics.o kvalues.o add_rm.cpp

dynamics: $(OBJ)
	$(C_COMPILER) $(OBJ) $(LIBS) -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
