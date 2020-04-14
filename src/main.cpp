#include <bits/stdc++.h>
#include "Solution.hpp"
#define NDEBUG
#include <assert.h>
#include <chrono>

using namespace std;

double z = 0.8; // constant related to the number of vertices that will be removed in the function shake()
int k = 5; // limit of iterations that dont improve the biclique
int total_iterations = 1000; // number of ils iterations
double total_time = 0, time_to_best = 0;

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
		Solution best_s(graph);

		// initialize C++ timer
        std::chrono::time_point<std::chrono::system_clock> start, end;

        int loop = 10, best_solution = 0, avarage_solution = 0, K = 3;
        // start timing
        start = std::chrono::system_clock::now();
        while(loop--) {
	        // start the first random solution
			s.generateRandomSolution();
			if(s.checkBicliqueSize() == false) s.balanceBiclique();

			int x = k, best_weight = s.getTotalWeight(), local_weight;
			Solution next_s(s);
			for(int iter = 0; iter < total_iterations; iter++) { // run ILS iterations
				next_s.VND(K);
				local_weight = next_s.getTotalWeight();
				if(local_weight > best_weight) {
					s = next_s;
					best_weight = local_weight;
					x = k;
				}	
				else x--;

				if(x == 0) {
					next_s.restartSolution();
					next_s.generateRandomSolution();
				}
				else next_s.shake(z);

				assert(next_s.checkIntegrity());
				assert(next_s.checkMu());
			}

			if(s.checkBicliqueSize() == false) {
				s.balanceBiclique();
				best_weight = s.getTotalWeight();
			}

			avarage_solution += best_weight;
			if(best_solution <= best_weight) {
				best_s = s;
				best_solution = best_weight;
				std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start; 
				time_to_best = elapsed_seconds.count();
			}

			assert(next_s.checkIntegrity());
			assert(next_s.checkMu());

			s.restartSolution();
		}

		// print final solution
		//s.printSolution();
		
		// end timing
        end = std::chrono::system_clock::now();

        // get difference
        std::chrono::duration<double> elapsed_seconds = end - start;

        //5 digits precision is enough
        total_time += elapsed_seconds.count();
        cout << best_solution << "\t" << avarage_solution / 10 << "\t" << std::setprecision(4) << total_time / 10.0 << "\t" << time_to_best << "\n";
        //best_s.printSolution(); 
        // execution time
	    //cout << "Execution time: " << std::setprecision(4) << total_time << " seconds." << endl;
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}