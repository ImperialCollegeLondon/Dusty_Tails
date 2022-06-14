C_COMPILER     = g++ -std=c++11
C_FLAGS        = -O3 -mavx -pg -g 
OUTPUT         = ./lightcurve
LIBS           = -lm

HEADERS =  opacities.h

OBJ = light_curve.o opacities.o

dynamics: $(OBJ)
	$(C_COMPILER) $(C_FLAGS)  $(OBJ) $(LIBS)  -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
