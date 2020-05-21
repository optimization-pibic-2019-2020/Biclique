SRC = ./src
TEST = ./test
DEBUG = ./debug

all: main.o Graph.o Solution.o
	g++ -O3 main.o Graph.o Solution.o -o main
	cp main \debug

main.o: $(SRC)/main.cpp
	g++ -c -O3 $(SRC)/main.cpp 

Graph.o: $(SRC)/Graph.cpp $(SRC)/Graph.hpp
	g++ -c -O3 $(SRC)/Graph.cpp

Solution.o: $(SRC)/Solution.cpp $(SRC)/Solution.hpp
	g++ -c -O3 $(SRC)/Solution.cpp

exec: 
	./main <$(TEST)/in

clear:
	rm main
	rm ./*.o

copy:
	cp main \debug

generatorW:
	python3 $(DEBUG)/generatorWeightedVertices.py

generatorU:
	python3 $(DEBUG)/generatorUniformVertices.py

tester:
	python3 $(DEBUG)/testerDimacs.py
	#python3 $(DEBUG)/testerBhoslib.py
	#python3 $(DEBUG)/testerBarabasiAlbert.py
	#python3 $(DEBUG)/testerErdosRenyi.py

	