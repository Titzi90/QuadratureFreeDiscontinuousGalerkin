########## COMPILER and FLAGS #################################
CC  = icc
CPP = icpc

CFLAGS = -O3 -g -Wall -ansi -std=c++11 -xhost -fopenmp
RM = rm -f

# If you want to make use of the likwid profiler
# uncommend the following two lines
# the evvirement variables LIKWID_LIB and LIKWID_INC has to be set
# on the rrze cluster this variables are already set
LIKWID_FLAGS  = -DUSE_LIKWID $(LIKWID_INC) -DLIKWID_PERFMON
LIKWIDLD_FLAGS = $(LIKWID_LIB) -llikwid -lm

CFLAGS += $(LIKWID_FLAGS)
LDFLAGS = $(LIKWIDLD_FLAGS)


GITVERSION = $(shell git rev-parse --short HEAD)
CFLAGS += -D_GITVERSION=\"$(GITVERSION)\"

########## VARIABLES ##########################################
BIN = test assembly analyticalTest analyticalConvergenceTest solidBodyRotation
OFILES_test = test.o polynomial.o Grid.o
OFILES_assembly = assembly.o polynomial.o VTKwriter.o Grid.o
OFILES_analyticalTest = analyticalTest.o polynomial.o VTKwriter.o Grid.o
OFILES_analyticalConvTest = analyticalConvergenceTest.o polynomial.o VTKwriter.o Grid.o
OFILES_solidBodyRotation =  solidBodyRotation.o polynomial.o VTKwriter.o Grid.o

########## TARGETS ############################################
.PHONY: all test clean
all: ${BIN}

clean:
	$(RM) $(BIN) $(OFILES_test) $(OFILES_assembly) ${OFILES_analyticalTest} \
				${OFILES_analyticalConvTest} ${OFILES_solidBodyRotation}

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

solidBodyRotation: $(OFILES_solidBodyRotation)
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
analyticalTest.o: analyticalTest.cpp stepper.hpp Grid.hpp VTKwriter.hpp assembly.hpp DataTypes.hpp
analyticalConvergenceTest.o: analyticalConvergenceTest.cpp stepper.hpp Grid.hpp VTKwriter.hpp assembly.hpp DataTypes.hpp monmomials_and_basefunctions.hpp
solidBodyRotation.o: solidBodyRotation.cpp stepper.hpp Grid.hpp VTKwriter.hpp assembly.hpp DataTypes.hpp
Grid.o: Grid.hpp DataTypes.hpp monmomials_and_basefunctions.hpp
