#include <bits/stdc++.h>
#include "Solution.hpp"
#define NDEBUG
#include <assert.h>
#include <chrono>
#include <random>

using namespace std;

int total_iterations = 5000; // number of ils iterations
int beta = 1000; // limit of iterations that dont improve the biclique
int alpha_calibration = 50; // variable related to the reactive grasp adjustment
vector<double> alphas{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; // all alpha values
vector<double> alphaProbability; // all alpha probabilities
vector<int> alphaTimesChosen;  // stores how many each alpha was chosen
vector<int> alphaSolutions;  // stores the sum of each solution for each alpha
double total_time = 0, time_to_best = 0, execution_time = 0;
int alpha_chosen; // variable to store the chosen alpha for the current iteration 


int main(int argc, char* argv[]) {
	try {
		// setting parameters
		if(argc > 1 && strcmp(argv[1], "-help") == 0) {
			cout << "-  " << "I + number of iterations for total iterations of GRASP+VND" << endl;
			cout << "-  " << "B + beta for the parameter beta" << endl;
			cout << "-  " << "C + alpha_calibration for the reactive grasp calibration parameter" << endl;
			cout << "-  " << "If a instruction is not set then the non-used paremeters will have the default value" << endl;
			return 0;
		}

		for(int iter = 1; iter < argc; iter++) {
			if(strcmp(argv[iter], "-I") == 0) total_iterations = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-B") == 0) beta = atoi(argv[++iter]); 
			else if(strcmp(argv[iter], "-C") == 0) alpha_calibration = atoi(argv[++iter]);
			else {
				cout << "Parameters Failed." << endl << "Type -help to use the parameters correctly" << endl;
				return 0;
			}
		}
	
		cout << "Total iterations = " << total_iterations << endl;
		cout << "Beta = " << beta << endl;
		cout << "Alpha Calibration = " << alpha_calibration << endl;
		
		// end of setting parameters

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

		Solution best_s(graph);

		// initialize C++ timer
        std::chrono::time_point<std::chrono::system_clock> start, execution_start, end;

        int loop = 10, best_solution = -1, avarage_solution = 0, K = 3;
        // start timing
        start = std::chrono::system_clock::now();

		cout << "Starting the algorithm...\n" << endl;
        while(loop--) {
			cout << abs(10 - loop) << " execution" << endl;

			Solution s(graph); // initialize all the variables and structures for solution

	        // start the first greedy random solution
			s.greedyRandomizedConstructive(0.0);
			if(s.checkBicliqueSize() == false) s.balanceBiclique();

			int x = beta, y = alpha_calibration, best_weight = s.getTotalWeight(), local_weight;
			double solutionAvarage, totalProbability = 0;
			cout << "Initial Solution: " << s.getTotalWeight() << endl;
			Solution next_s(graph);

			// starting random generator for discrete distribution
			random_device device;
			mt19937 generator(device());

			// initializing alpha vectors
			for(int i = 0; i < 11; i++) {
				alphaProbability.push_back(1.0); 
				alphaTimesChosen.push_back(0);  
				alphaSolutions.push_back(0);
			}

			// start execution timing
			execution_start = std::chrono::system_clock::now();

			for(int iter = 0; iter < total_iterations; iter++) { // run ILS iterations
				// choose the alpha based on reactive grasp
				discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
				alpha_chosen = distribution(generator);
				alphaTimesChosen[alpha_chosen]++;

				next_s.greedyRandomizedConstructive(alphas[alpha_chosen]); // starts a new optimal solution
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();

				next_s.VND(K);
				local_weight = next_s.getTotalWeight();
				alphaSolutions[alpha_chosen] += local_weight;

				if(local_weight > best_weight) {
					s = next_s;
					best_weight = local_weight;
					x = beta;
					std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count();
					cout << "New best found: " << best_weight << endl;
				}	
				else if(local_weight == best_weight) x--;

				if(x == 0) break; // beta parameter

				y--; // parameter related to reactive grasp calibration 

				if(y == 0) { // reactive grasp (alpha calibration)
					y = alpha_calibration; // reset the parameter

					for(int i = 0; i < 11; i++) {
						if(alphaTimesChosen[i] > 0) {
							solutionAvarage = ((double) alphaSolutions[i]) / ((double) alphaTimesChosen[i]); // avarage of alpha solutions
							alphaProbability[i] =  ((double) best_weight) / solutionAvarage; // best weight divided by the avarage of the alpha solutions so far
						} 
					}
				}

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
			}

			assert(next_s.checkIntegrity());
			assert(next_s.checkMu());

			// execution informations
			discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
			std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();
			
			cout << "\nAlpha probabilities (Reactive grasp)" << endl;
			for(int i = 0; i < 11; i++) {
				cout << "Alpha: " << ((double) i) / 10.0 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
			}

			cout << "\nResults:\n";
			cout << "Best: " << best_weight << endl;
			cout << "Total time of execution: " << execution_time << endl;
			cout << "Time to best: " << time_to_best << "s\n" << endl;

			// clearing alpha vectors
			alphaProbability.clear(); 
			alphaTimesChosen.clear();  
			alphaSolutions.clear();
		}
		
		cout << "End of the algorithm" << endl;

		// end timing
        end = std::chrono::system_clock::now();

        // get difference
        std::chrono::duration<double> elapsed_seconds = end - start;

        //5 digits precision is enough
        total_time += elapsed_seconds.count();

		cout << "Best Solution: " << best_solution << endl;
		cout << "Avarage Solution: " << avarage_solution / 10 << endl;
		cout << "Avarage Time: " << std::setprecision(4) << total_time / 10 << endl;
		cout << "Solution: ";
		best_s.printSolution();
		cout << "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" << endl;		
	}
	catch (std::exception &e) {
		cerr << e.what();
	}

	return 0;
}