#include <bits/stdc++.h>
#include "Solution.hpp"

using namespace std;

int main() {
	try {
		int v, e;
		cin >> v >> e;

		if(v == 0 || e == 0) { // The algorithm cannot run if the number of vertices or edges is equal to zero
			cout << "The number of vertices or edges is equal to 0!" << endl;
			return 0;
		}

		Graph *graph = new Graph(v, e);
		graph->readWeight(); // read all the weight and put into the weight vector
		graph->readEdges(); // read all the edges and put into the adjList
		graph->sort(); 
		//graph->showGraphInformations(); // show all the informations of the input Graph
		Solution s(graph); // initialize all the variables and structures for solution
		for(int idx = 0; idx < 1000; idx++) {
			s.generateRandomSolution();
			cout << "Initial Random Solution:" << endl;
			s.printSolution();
			cout << endl;
			s.restartSolution();
		}

		/*// print the random solution
		cout << "Initial Random Solution:" << endl;
		s.printSolution();*/
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}