#include <iostream>
#include <vector>
#include <string.h> 
#include <string> 
#include <iomanip>

#ifndef SOLUTION_HPP
#define SOLUTION_HPP
    
#include "Solution.hpp"

#endif

#include "NonBipartiteSolution.hpp"
#include "BipartiteSolution.hpp"
//#define NDEBUG
#include <assert.h>
#include <chrono>
#include <random>
#include <algorithm> 

using namespace std::chrono;
using namespace std;

// GRAPH Variables 
int input_type, partitionA_size, partitionB_size, v, e;
int verticesRemoved = 0, edgesRemoved = 0;
int maxVerticesRemoved = 0, maxEdgesRemoved = 0;
vector<bool> vertexInGraph; 

// GRASP-VND PARAMETERS
int total_iterations = -1; //5000; // number of grasp iterations 
int loop = 10; 
double execution_time_limit = 60; // time of each execution
int beta_var = 1000000; // TODO --- 2347; // limit of iterations that dont improve the biclique
int alpha_calibration = 231; // variable related to reactive grasp adjustment
int target = -1; // variable related to biclique target
int K = 2; // 3; // number of neighborhood structures

// REACTIVE GRASP VARIABLES AND VECTORS (to help each alpha probability calculation)
vector<double> alphas{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; // all alpha values
vector<double> alphaProbability; // all alpha probabilities
vector<int> alphaTimesChosen;  // stores how many each alpha was chosen
vector<int> alphaSolutions;  // stores the sum of each solution for each alpha
int alpha_chosen; // variable to store the chosen alpha for the current iteration 

// EXECUTION TIME VARIABLES 
double total_time = 0, time_to_best = 0, execution_time = 0, timeLimit, time_to_reduction;

// HYPERPARAMETERS FOR ADAPTIVE MEMORY
int amBetaParamater = 11586;
int amEpsilonParameter = 58768;

// constant related of how much time an individual edge should waste in reduction
// to calculate timeLimit (variable related to time limit in reduction) you need to multiply the quantity of edge (e) per edgeTimeConst
double edgeTimeConst = 0.2614115754; 

void BipartiteReactiveGrasp() {
	try {
		Graph original_graph(v, e);
		original_graph.readWeight(); // read all the weight and put into the weight vector
		original_graph.readEdges(); // read all the edges and put into the adjList
		original_graph.sort(); 
		original_graph.initializeAm();
		original_graph.initializeHIndex();
		original_graph.initializeAccumulatedSum();

		Graph graph = original_graph; // copy original_graph to graph
		//graph.showGraphInformations(); // show all the informations of the input Graph
	
		BipartiteSolution best_s(&graph, partitionA_size, partitionB_size);
		
		// initialize C++ timer
        time_point<system_clock> start, execution_start, end;
		duration<double> elapsed_seconds;

        int best_solution = -1, avarage_solution = 0;

        // start timing
        start = system_clock::now();
		
		cout << "Starting the algorithm...\n" << endl;

        while(loop--) {
			cout << abs(10 - loop) << " execution" << endl;

			BipartiteSolution s(&graph, partitionA_size, partitionB_size); // initialize all the variables and structures for solution

	        // start the first greedy random solution
			s.greedyRandomizedConstructive(0.0);
			if(s.checkBicliqueSize() == false) s.balanceBiclique();
			
			int x = beta_var, y = alpha_calibration, best_weight = s.getTotalWeight(), local_weight, iter = 0;
			double solutionAvarage;
			cout << "Initial Solution: " << s.getTotalWeight() << endl;
			
			BipartiteSolution next_s(&graph, partitionA_size, partitionB_size);

			// starting random generator for discrete distribution
			random_device device;
			mt19937 generator(device());

			// initializing alpha vectors
			for(int i = 0; i < 11; i++) {
				alphaProbability.push_back(1.0); 
				alphaTimesChosen.push_back(0);  
				alphaSolutions.push_back(0);
			}

			// initializing vertex removed vector
			vertexInGraph.resize(v);
			fill(vertexInGraph.begin(), vertexInGraph.end(), true); // initialize the vector with all values set to true meaning that all the vertex are in the graph

			// start execution timing
			execution_start = system_clock::now();
			elapsed_seconds = system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			while((iter < total_iterations) || (execution_time <= execution_time_limit)) { // run ILS iterations
				// choose the alpha based on reactive grasp
				discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
				alpha_chosen = distribution(generator);
				alphaTimesChosen[alpha_chosen]++;

				next_s.greedyRandomizedConstructive(alphas[alpha_chosen]); // starts a new optimal solution
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();
				
				next_s.VND(K);
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();
				local_weight = next_s.getTotalWeight();
				assert(next_s.checkIntegrity());

				// avoiding a problem where many vertices are removed then local weight always results in 0
				// then all probabilites goes to 0 or -nan
				if(local_weight != 0) { alphaSolutions[alpha_chosen] += local_weight; }
				else { alphaTimesChosen[alpha_chosen]--; }

				if(local_weight > best_weight) {	
					elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count(); // update time_to_best 
					cout << "New best found: " << best_weight << endl;

					next_s.restartAm(amBetaParamater);
					amBetaParamater +=  amBetaParamater * amEpsilonParameter;

					// when local_weight surpass target then go to the next iteration
					if(target != -1 && target <= local_weight) { 
						next_s.restartSolution(vertexInGraph);
						assert(next_s.checkIntegrity());
		
						break; 
					} 
					
					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();

					time_to_reduction = (timeLimit < execution_time_limit - execution_time) ? timeLimit : execution_time_limit - execution_time;
					next_s.reduceGraph(vertexInGraph, best_weight, -1, time_to_reduction); // start graph reduction

					s = next_s; // update best local solution
					best_weight = local_weight; // update best_weight
					x = beta_var; // reinitialize parameter x
				} else { 
					next_s.updateAm();
				}	
				
				if(local_weight == best_weight) x--;

				if(x == 0) { break; }// beta_var parameter

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

				if(total_iterations <= 0) {
					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();
				}
				else {
					iter++;
				}

				next_s.restartSolution(vertexInGraph);
				assert(next_s.checkIntegrity());
			}

			// checking if the algorithm is trapped in a loop 
			assert(execution_time < 65);

			// restarting graph
			graph = original_graph;

			if(s.checkBicliqueSize() == false) {
				s.balanceBiclique();
				best_weight = s.getTotalWeight();
			}

			avarage_solution += best_weight;
			if(best_solution < best_weight) {
				best_s = s;
				best_solution = best_weight;
			}

			assert(best_s.checkIntegrity());

			// execution informations
			discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
			elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			cout << "\nAlpha probabilities (Reactive grasp)" << endl;
			for(int i = 0; i < 11; i++) {
				cout << "Alpha: " << ((double) i) / 10.0 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
			}

			verticesRemoved = s.getRemovedVertices();
			edgesRemoved = s.getRemovedEdges();

			maxVerticesRemoved = max(maxVerticesRemoved, verticesRemoved);
			maxEdgesRemoved = max(maxEdgesRemoved, edgesRemoved);

			cout << "\nGraph Reduce results:" << endl;
			cout << "Vertices removed: " << verticesRemoved << " of " << v << endl;
			cout << "Edges removed: " << edgesRemoved << " of " << e << endl;
			cout << ((verticesRemoved * 1.0) / (v * 1.0)) * 100.0 << "% of vertices removed and " << ((edgesRemoved * 1.0) / (e * 1.0))* 100.0 << "% of edges removed"<< endl;

			cout << "\nResults:\n";
			cout << "Best: " << best_weight << endl;
			cout << "Total time of execution: " << execution_time << endl;
			cout << "Time to best: " << time_to_best << "s\n" << endl;
			//cout << best_weight << "," << execution_time << endl;
			//cout << best_weight << "," << time_to_best << endl;

			// clearing alpha vectors
			alphaProbability.clear(); 
			alphaTimesChosen.clear();  
			alphaSolutions.clear();

			// restarting removed vertices stats
			verticesRemoved = 0;
			edgesRemoved = 0;
		}
		
		cout << "End of the algorithm" << endl;

		// end timing
        end = std::chrono::system_clock::now();

        // get difference
        elapsed_seconds = end - start;

        //5 digits precision is enough
        total_time += elapsed_seconds.count();

		cout << "Best Solution: " << best_solution << endl;
		cout << "Avarage Solution: " << avarage_solution / 10 << endl;
		cout << "Avarage Time: " << std::setprecision(4) << total_time / 10 << endl;
		cout << "Solution: ";
		best_s.printSolution();
		cout << "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" << endl;		
		//cout << fixed << setprecision(3) << double(((maxVerticesRemoved * 1.0) / (v * 1.0)) * 100.0) << endl; // Using this print to get the results to generate the graphic
	}
	catch (std::exception &e) {
		cerr << e.what();
	}
}


void NonBipartiteReactiveGrasp() {
	try {
		Graph original_graph(v, e);
		original_graph.readWeight(); // read all the weight and put into the weight vector
		original_graph.readEdges(); // read all the edges and put into the adjList
		original_graph.sort(); 
		original_graph.initializeHIndex();
		original_graph.initializeAccumulatedSum();

		Graph graph = original_graph; // copy original_graph to graph
		//graph.showGraphInformations(); // show all the informations of the input Graph
	
		NonBipartiteSolution best_s(&graph, partitionA_size, partitionB_size);
		
		// initialize C++ timer
        time_point<system_clock> start, execution_start, end;
		duration<double> elapsed_seconds;

        int best_solution = -1, avarage_solution = 0;

        // start timing
        start = system_clock::now();
		
		cout << "Starting the algorithm...\n" << endl;

        while(loop--) {
			cout << abs(10 - loop) << " execution" << endl;

			NonBipartiteSolution s(&graph, partitionA_size, partitionB_size); // initialize all the variables and structures for solution

	        // start the first greedy random solution
			s.greedyRandomizedConstructive(0.0);
			if(s.checkBicliqueSize() == false) s.balanceBiclique();
			
			int x = beta_var, y = alpha_calibration, best_weight = s.getTotalWeight(), local_weight, iter = 0;
			double solutionAvarage;
			cout << "Initial Solution: " << s.getTotalWeight() << endl;
			
			NonBipartiteSolution next_s(&graph, partitionA_size, partitionB_size);

			// starting random generator for discrete distribution
			random_device device;
			mt19937 generator(device());

			// initializing alpha vectors
			for(int i = 0; i < 11; i++) {
				alphaProbability.push_back(1.0); 
				alphaTimesChosen.push_back(0);  
				alphaSolutions.push_back(0);
			}

			// initializing vertex removed vector
			vertexInGraph.resize(v);
			fill(vertexInGraph.begin(), vertexInGraph.end(), true); // initialize the vector with all values set to true meaning that all the vertex are in the graph

			// start execution timing
			execution_start = system_clock::now();
			elapsed_seconds = system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			while((iter < total_iterations) || (execution_time <= execution_time_limit)) { // run ILS iterations
				// choose the alpha based on reactive grasp
				discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
				alpha_chosen = distribution(generator);
				alphaTimesChosen[alpha_chosen]++;

				next_s.greedyRandomizedConstructive(alphas[alpha_chosen]); // starts a new optimal solution
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();
				
				next_s.VND(K);
				local_weight = next_s.getTotalWeight();

				// avoiding a problem where many vertices are removed then local weight always results in 0
				// then all probabilites goes to 0 or -nan
				if(local_weight != 0) { alphaSolutions[alpha_chosen] += local_weight; }
				else { alphaTimesChosen[alpha_chosen]--; } 

				if(local_weight > best_weight) {
					s = next_s; // update best local solution
					best_weight = local_weight; // update best_weight
					x = beta_var; // reinitialize parameter x

					elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count(); // update time_to_best 
					cout << "New best found: " << best_weight << endl;
					// TODO create time limit mechanism (just made in the bipartite solution)
					// TODO add target conditional

					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();
					//next_s.reduceGraph(vertexInGraph, best_weight, -1, timeLimit - execution_time); // start graph reduction
				}	
				else if(local_weight == best_weight) x--;

				if(x == 0) break; // beta_var parameter

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

				if(total_iterations == 0) {
					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();
				}
				else {
					iter++;
				}

				next_s.restartSolution(vertexInGraph);
				assert(next_s.checkIntegrity());
				assert(next_s.checkMu());
			}

			// restarting graph
			graph = original_graph;

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
			elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			cout << "\nAlpha probabilities (Reactive grasp)" << endl;
			for(int i = 0; i < 11; i++) {
				cout << "Alpha: " << ((double) i) / 10.0 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
			}

			verticesRemoved = s.getRemovedVertices();
			edgesRemoved = s.getRemovedEdges();

			maxVerticesRemoved = max(maxVerticesRemoved, verticesRemoved);
			maxEdgesRemoved = max(maxEdgesRemoved, edgesRemoved);

			cout << "\nGraph Reduce results:" << endl;
			cout << "Vertices removed: " << verticesRemoved << " of " << v << endl;
			cout << "Edges removed: " << edgesRemoved << " of " << e << endl;
			cout << ((verticesRemoved * 1.0) / (v * 1.0)) * 100.0 << "% of vertices removed and " << ((edgesRemoved * 1.0) / (e * 1.0))* 100.0 << "% of edges removed"<< endl;

			cout << "\nResults:\n";
			cout << "Best: " << best_weight << endl;
			cout << "Total time of execution: " << execution_time << endl;
			cout << "Time to best: " << time_to_best << "s\n" << endl;

			// clearing alpha vectors
			alphaProbability.clear(); 
			alphaTimesChosen.clear();  
			alphaSolutions.clear();

			// restarting removed vertices stats
			verticesRemoved = 0;
			edgesRemoved = 0;
		}
		
		cout << "End of the algorithm" << endl;

		// end timing
        end = std::chrono::system_clock::now();

        // get difference
        elapsed_seconds = end - start;

        //5 digits precision is enough
        total_time += elapsed_seconds.count();

		cout << "Best Solution: " << best_solution << endl;
		cout << "Avarage Solution: " << avarage_solution / 10 << endl;
		cout << "Avarage Time: " << std::setprecision(4) << total_time / 10 << endl;
		cout << "Solution: ";
		best_s.printSolution();
		cout << "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" << endl;		
			
		//cout << fixed << setprecision(3) << double(((maxVerticesRemoved * 1.0) / (v * 1.0)) * 100.0) << endl; // Using this print to get the results to generate the graphic
	}
	catch (std::exception &e) {
		cerr << e.what();
	}
}

int main(int argc, char* argv[]) {
	try {
		// setting parameters
		if(argc > 1 && strcmp(argv[1], "-help") == 0) {
			cout << "-  " << "You can use only one command (I or T) for execution iteration" << endl;
			cout << "-  " << "-I + number of iterations for total iterations of GRASP+VND" << endl;
			cout << "-  " << "-T + time in seconds to define time per execution of GRASP+VND" << endl;
			cout << "-  " << "-G + an integer to define a target for biclique" << endl;
			cout << "-  " << "-B + beta_var for parameter beta_var" << endl;
			cout << "-  " << "-C + alpha_calibration for reactive grasp calibration parameter" << endl;
			cout << "-  " << "-amB + amBetaParamater for adaptive memory parameter" << endl;
			cout << "-  " << "-amE + amEpsilonParameter for adaptive memory parameter" << endl;
			cout << "-  " << "-eTc + edgeTimeConst for reduction optimization parameter" << endl;
			cout << "-  " << "If a instruction is not set then the non-used paremeters will have the default value" << endl;
			return 0;
		}

		for(int iter = 1; iter < argc; iter++) {
			if(strcmp(argv[iter], "-I") == 0) {
				total_iterations = atoi(argv[++iter]);
				execution_time_limit = 0.0;
			}
			else if(strcmp(argv[iter], "-T") == 0) {
				execution_time_limit = stod(argv[++iter]);
				total_iterations = 0;
			}
			else if(strcmp(argv[iter], "-B") == 0) beta_var = atoi(argv[++iter]); 
			else if(strcmp(argv[iter], "-G") == 0) target = atoi(argv[++iter]); 
			else if(strcmp(argv[iter], "-C") == 0) alpha_calibration = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-amB") == 0) amBetaParamater = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-amE") == 0) amEpsilonParameter = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-eTc") == 0) edgeTimeConst = stod(argv[++iter]);
			else {
				cout << "Parameters Failed." << endl << "Type -help to use the parameters correctly" << endl;
				return 0;
			}
		}
	
		//cout << "Total iterations = " << total_iterations << endl;
		//cout << "Execution Time Limit = " << execution_time_limit << endl;
		//cout << "Beta = " << beta_var << endl;
		//cout << "Alpha Calibration = " << alpha_calibration << endl;
		
		// end of setting parameters

		cin >> input_type;

		if(input_type == 1) { // not bipartitioned instance
			cin >> v >> e;
			partitionA_size = v;
			partitionB_size = v;
		}
		else { // bipartitioned instance
			cin >> partitionA_size >> partitionB_size >> e;
			v = partitionA_size + partitionB_size;
		}
		

		if(v == 0 || e == 0) { // The algorithm cannot run if the number of vertices or edges is equal to zero
			cout << "The number of vertices or edges is equal to 0!" << endl;
			return 0;
		}

		timeLimit = ceil(edgeTimeConst * e);
		(input_type == 1) ? NonBipartiteReactiveGrasp() : BipartiteReactiveGrasp();
	}	
	catch (std::exception &e) {
		cerr << e.what();
	}	

	return 0;
}