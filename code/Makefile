########## COMPILER and FLAGS #################################
CC  = gcc
CPP = g++

CFLAGS = -O3 -g -Wall -ansi -std=c++11
RM = rm -f

########## VARIABLES ##########################################
BIN = test assembly analyticalTest
OFILES_test = test.o polynomial.o Grid.o
OFILES_assembly = assembly.o polynomial.o VTKwriter.o Grid.o
OFILES_analyticalTest = analyticalTest.o polynomial.o VTKwriter.o Grid.o
OFILES_analyticalConvTest = analyticalConvergenceTest.o polynomial.o VTKwriter.o Grid.o

########## TARGETS ############################################
.PHONY: all test clean
all: ${BIN}

clean:
	$(RM) $(BIN) $(OFILES_test) $(OFILES_assembly) ${OFILES_analyticalTest} ${OFILES_analyticalConvTest}

test: test_compile
	./test

########## LINKING RULEs ######################################
test_compile: $(OFILES_test)
	$(CPP) $(CFLAGS) -o test $^ $(LDFLAGS)

assembly: $(OFILES_assembly)
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS)

analyticalTest: $(OFILES_analyticalTest)
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS)

analyticalConvergenceTest: $(OFILES_analyticalConvTest)
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
analyticalTest.o: analyticalTest.cpp stepper.hpp Grid.hpp VTKwriter.hpp assembly.hpp
analyticalConvergenceTest.o: analyticalConvergenceTest.cpp stepper.hpp Grid.hpp VTKwriter.hpp assembly.hpp
Grid.o: Grid.hpp DataTypes.hpp
