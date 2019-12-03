SRC = ./src
OBJ = ./obj
TEST = ./test

all: main.o Graph.o Solution.o
	g++ main.o Graph.o Solution.o -o main

main.o: $(SRC)/main.cpp
	g++ -c $(SRC)/main.cpp

Graph.o: $(SRC)/Graph.cpp $(SRC)/Graph.hpp
	g++ -c $(SRC)/Graph.cpp

Solution.o: $(SRC)/Solution.cpp $(SRC)/Solution.hpp
	g++ -c $(SRC)/Solution.cpp

exec: 
	./main <$(TEST)/in

clean:
	rm main
	rm ./*.o