#TODO: look at each tester/generator and update the path
SRC = ./src
TEST = ./test
GENERATOR = ./debug/scripts/generators
TESTER = ./debug/scripts/testers
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
	python3 $(GENERATORS)/generatorWeightedVertices.py

generatorU:
	python3 $(GENERATORS)/generatorUniformVertices.py

tester:
	$(TESTER)/konect_script.sh
	#python3 $(TESTER)/testerKonect.py
	#python3 $(TESTER)/testerGNP.py
	#python3 $(TESTER)/testerGNPTest.py
	#python3 $(TESTER)/testerGNPBTest.py
	#python3 $(TESTER)/testerDimacs.py
	#python3 $(TESTER)/testerBhoslib.py
	#python3 $(TESTER)/testerBarabasiAlbert.py
	#python3 $(TESTER)/testerErdosRenyi.py

	