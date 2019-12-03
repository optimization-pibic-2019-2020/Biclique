#include <bits/stdc++.h>
#include "Solution.hpp"

using namespace std;

int main() {
	try {
		int v, e;
		cin >> v >> e;

		if(v == 0 || e == 0) { // The algorithm can´t run if the number of vertices or edges is equal to zero
			cout << "The number of vertices or edges is equal to 0!" << endl;
			return 0;
		}

		Graph *graph = new Graph(v, e);
		graph->readEdges(); // reads all the edges and put into the adjList
		graph->readWeight(); // read all the weight and put into the weight vector
		graph->sort(); 
		//graph->showGraphInformations(); // It is used to see all the informations of the input Graph
		Solution s(graph);
		s.generateRandomSolution();
		cout << "Initial Random Solution:" << endl;
		s.printSolution();
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}