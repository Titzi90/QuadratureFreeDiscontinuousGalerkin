########## COMPILER and FLAGS #################################
CC  = gcc
CPP = g++

CFLAGS = -O3 -g -Wall -ansi -std=c++11
RM = rm -f

########## VARIABLES ##########################################
BIN = polyTest
OFILES_poly_test = test.o polynomial.o
OFILES_assembly = assembly.o polynomial.o VTKwriter.o

########## TARGETS ############################################
.PHONY: all test clean
all: test assembly

clean:
	$(RM) $(BIN) $(OFILES_poly_test) $(OFILES_assembly)

test: polyTest
	./polyTest

########## LINKING RULEs ######################################
polyTest: $(OFILES_poly_test)
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS)

assembly: $(OFILES_assembly)
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS)

##########  GENERAL RULEs #####################################
.cpp.o:
	$(CPP) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<


########## DEPENDENCIES ######################################
test.o: test.cpp polynomial.hpp monmomials_and_basefunctions.hpp Grid.hpp DataTypes.hpp assembly.hpp
polynomial.o: polynomial.cpp polynomial.hpp monmomials_and_basefunctions.hpp
assembly.o: assembly.hpp assembly.cpp polynomial.hpp monmomials_and_basefunctions.hpp Grid.hpp DataTypes.hpp
VTKwriter.o: VTKwriter.cpp VTKwriter.hpp Grid.hpp DataTypes.hpp
