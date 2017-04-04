########## COMPILER and FLAGS #################################
CC  = gcc
CPP = g++

CFLAGS = -O3 -g -Wall -ansi -std=c++11
RM = rm -f

########## VARIABLES ##########################################
BIN = polyTest
OFILES_poly_test = test.o polynomial.o

########## TARGETS ############################################
.PHONY: all test clean
all: test

clean:
	$(RM) $(BIN) $(OFILES_poly_test)

test: polyTest
	./polyTest

########## LINKING RULEs ######################################
polyTest: $(OFILES_poly_test)
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

##########  GENERAL RULEs #####################################
.cpp.o:
	$(CPP) $(CFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<


########## DEPENDENCIES ######################################
test.o: test.cpp polynomial.hpp monmomials_and_basefunctions.hpp
polynomial.o: polynomial.cpp polynomial.hpp
