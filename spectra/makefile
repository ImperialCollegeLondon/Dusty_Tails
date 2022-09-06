C_COMPILER     = g++ -std=c++11
C_FLAGS        = -O3 -mavx -pg -g 
OUTPUT         = ./spectra
LIBS           = -lm

HEADERS = opacities_s.h 

OBJ = spectra.o

dynamics: $(OBJ)
	$(C_COMPILER) $(C_FLAGS)  $(OBJ) $(LIBS)  -o $(OUTPUT)

.cpp.o:
	$(C_COMPILER) $(C_FLAGS) -c $<


clean:
	@rm -f *.o

$(OBJ) : $(HEADERS)
