CPLEXFLAGS=-O3 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -I/opt/cplex/cplex/include -I/opt/cplex/concert/include  -L/opt/cplex/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L/opt/cplex/concert/lib/x86-64_linux/static_pic -lconcert -lm -pthread -std=c++0x -ldl

CFLAGS=-std=c++11 -static-libstdc++ -static-libgcc -Wall

all:maximumBalancedBicliqueProblem.o
	g++  maximumBalancedBicliqueProblem.o -o maximumBalancedBicliqueProblem.run $(CPLEXFLAGS) $(CFLAGS)


%.o: %.cpp %.hpp
	g++ -c $< -o $@ $(CFLAGS)

maximumBalancedBicliqueProblem.o: maximumBalancedBicliqueProblem.cpp
	g++ -c -o maximumBalancedBicliqueProblem.o maximumBalancedBicliqueProblem.cpp $(CPLEXFLAGS) $(CFLAGS)

clean:
	rm -f *.o

