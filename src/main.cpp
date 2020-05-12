#include <bits/stdc++.h>
#include "Solution.hpp"
#define NDEBUG
#include <assert.h>
#include <chrono>

using namespace std;

int k; // limit of iterations that dont improve the biclique (beta)
double p = 0.0; // variable related to the quality based construction in the RCL list (alpha)
int total_iterations; // number of ils iterations
double total_time = 0, time_to_best = 0;


int main() {
	try {
		int v, e;
		cin >> v >> e;

		total_iterations = e / 100; // 1% of the edges
		if(total_iterations > 20000) total_iterations = 20000; 
		else if(total_iterations < 1000) total_iterations = 3000;

		k = total_iterations / 2;

		if(v == 0 || e == 0) { // The algorithm cannot run if the number of vertices or edges is equal to zero
			cout << "The number of vertices or edges is equal to 0!" << endl;
			return 0;
		}

		Graph *graph = new Graph(v, e);
		graph->readWeight(); // read all the weight and put into the weight vector
		graph->readEdges(); // read all the edges and put into the adjList
		graph->sort(); 
		//graph->showGraphInformations(); // show all the informations of the input Graph

		Solution best_s(graph);

		// initialize C++ timer
        std::chrono::time_point<std::chrono::system_clock> start, end;

        int loop = 10, best_solution = -1, avarage_solution = 0, K = 3;
        // start timing
        start = std::chrono::system_clock::now();
        while(loop--) {
			Solution s(graph); // initialize all the variables and structures for solution
	        // start the first greedy random solution
			s.greedyRandomizedConstructive(p);
			if(s.checkBicliqueSize() == false) s.balanceBiclique();

			int x = k, best_weight = s.getTotalWeight(), local_weight;
			Solution next_s(graph);
			for(int iter = 0; iter < total_iterations; iter++) { // run ILS iterations
				next_s.greedyRandomizedConstructive(p); // starts a new optimal solution
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();

				next_s.VND(K);
				local_weight = next_s.getTotalWeight();
				if(local_weight > best_weight) {
					s = next_s;
					best_weight = local_weight;
					x = k;
				}	
				else if(local_weight == best_weight) x--;

				if(x == 0) break;
				
				p += 0.01;
				p = fmod(p, 0.11);

				next_s.restartSolution();
				assert(next_s.checkIntegrity());
				assert(next_s.checkMu());
			}

			if(s.checkBicliqueSize() == false) {
				s.balanceBiclique();
				best_weight = s.getTotalWeight();
			}

			avarage_solution += best_weight;
			if(best_solution < best_weight) {
				best_s = s;
				best_solution = best_weight;
				std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start; 
				time_to_best = elapsed_seconds.count();
			}

			assert(next_s.checkIntegrity());
			assert(next_s.checkMu());
		}
		
		// end timing
        end = std::chrono::system_clock::now();

        // get difference
        std::chrono::duration<double> elapsed_seconds = end - start;

        //5 digits precision is enough
        total_time += elapsed_seconds.count();
        cout << best_solution << "\t" << avarage_solution / 10 << "\t" << std::setprecision(4) << total_time / 10 << "\t" << time_to_best << "\t";
		best_s.printSolution();
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}