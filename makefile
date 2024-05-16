# Select compiler
CMPR = g++

OBJDIR = ./obj/

CPPFLAGS=$(ARG)
CPPFLAGS +=  

CPP = 
CFLAGS = -std=c++17 -Wall

OBJ = main.o mesh.o misc.o initconds.o timestep.o setbcs.o fluxSolver.o flags.o output.o  
EXEC = main.x

OBJS = $(addprefix $(OBJDIR), $(OBJ))

$(addprefix ./obj/, %.o):%.cpp 
	$(CMPR) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

# Link into an executable
main: $(OBJS)
	$(CMPR) $(OBJS) -o $(EXEC) 


clean:
	rm -f ./obj/*.o 
	rm -f ./*.x
cleandata:	
	rm -f ./data/*.dat
