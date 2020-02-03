CC = g++
CXX = g++
CXXFLAGS = -Ofast
LDLIBS += -lntl -lgmp -lm -lpthread

PROGS = main Samplers
//PROGS = test_lwe_oracle test_discrete_gaussian bkw

.PHONY: all clean

all: $(PROGS)

//test_discrete_gaussian: test_discrete_gaussian.o discrete_gaussian.o
//test_discrete_gaussian.o: test_discrete_gaussian.cpp discrete_gaussian.hpp
//discrete_gaussian.o: discrete_gaussian.cpp discrete_gaussian.hpp

main: main.o Samplers.o
main.o: main.cpp Samplers.h
Samplers.o:Samplers.cpp Samplers.h

clean:
	rm -rf *.o $(PROGS) a.out

