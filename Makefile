#TODO: look at each tester/generator and update the path
SRC = ./src
TEST = ./test
DEBUG = ./debug
OUT = ./out
 
all: main.o Graph.o Solution.o NonBipartiteSolution.o BipartiteSolution.o
	g++ -O3 -Wall main.o Graph.o Solution.o NonBipartiteSolution.o BipartiteSolution.o -o main
	cp ./*.o $(OUT)
	rm ./*.o

main.o: $(SRC)/main.cpp
	g++ -c -O3 -Wall $(SRC)/main.cpp 

Graph.o: $(SRC)/Graph.cpp $(SRC)/Graph.hpp
	g++ -c -O3 -Wall $(SRC)/Graph.cpp 

Solution.o: $(SRC)/Solution.cpp $(SRC)/Solution.hpp
	g++ -c -O3 -Wall $(SRC)/Solution.cpp

NonBipartiteSolution.o: $(SRC)/NonBipartiteSolution.cpp $(SRC)/NonBipartiteSolution.hpp
	g++ -c -O3 -Wall $(SRC)/NonBipartiteSolution.cpp

BipartiteSolution.o: $(SRC)/BipartiteSolution.cpp $(SRC)/BipartiteSolution.hpp
	g++ -c -O3 -Wall $(SRC)/BipartiteSolution.cpp

exec: 
	./main <$(TEST)/in

clear:
	rm main
	rm ./*.o

generatorW:
	python3 $(DEBUG)/generatorWeightedVertices.py

generatorU:
	python3 $(DEBUG)/generatorUniformVertices.py

tester:
	python3 $(DEBUG)/testerKonect.py
	python3 $(DEBUG)/testerGNP.py
	#python3 $(DEBUG)/testerGNPTest.py
	#python3 $(DEBUG)/testerGNPBTest.py
	#python3 $(DEBUG)/testerDimacs.py
	#python3 $(DEBUG)/testerBhoslib.py
	#python3 $(DEBUG)/testerBarabasiAlbert.py
	#python3 $(DEBUG)/testerErdosRenyi.py

	