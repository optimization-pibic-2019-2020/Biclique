#include <bits/stdc++.h>
#include "Solution.hpp"
#define NDEBUG
#include <assert.h>
using namespace std;

int z = 4; // constant related to the number of vertices that will be removed in the function shake()
int k = 20000; // limit of iterations that dont improve the biclique
int total_iterations = 200000; // number of ils iterations

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

		s.generateRandomSolution();
		cout << "Initial Random Solution:" << endl;
		s.printSolution();

		int x = k, best_weight = s.getTotalWeight(), local_weight;
		Solution next_s(s);
		for(int iter = 0; iter < total_iterations; iter++) { // run ILS iterations
			assert(next_s.checkIntegrity());
			next_s.VND();
			local_weight = next_s.getTotalWeight();
			if(local_weight > best_weight) {
				if(next_s.checkBicliqueSize()) {
					s = next_s;
					best_weight = local_weight;
					x = k;
				}
				else {
					next_s.balanceBiclique();
					if(local_weight > best_weight) {
						s = next_s;
						best_weight = local_weight;
						x = k;
					}
					else x--;
				}
			}
			else x--;

			if(x == 0) next_s.generateRandomSolution();
			//else next_s.shake(z);
		}

		if(s.checkBicliqueSize() == false) s.balanceBiclique();

		cout << "Final Solution:" << endl;
		s.printSolution();
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}