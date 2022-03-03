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
#define NDEBUG
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

// TABU PARAMETERS
int total_iterations = 0;//5000; // number of tabu iterations
int loop = 10; // variable related to the total number of executions
double execution_time_limit = 60; // time of each execution 
int target = -1; // variable related to biclique target

// Tabu VARIABLES
int iterWithoutImprovement = 143;
int tabuLimitIteration = 0;

// EXECUTION TIME VARIABLES 
double total_time = 0, time_to_best = 0, execution_time = 0, timeLimit;

// constant related of how much time an individual edge should waste in reduction
// to calculate timeLimit (variable related to time limit in reduction) you need to multiply the quantity of edge (e) per edgeTimeConst
double edgeTimeConst = 0.000005586; 

void BipartiteReactiveGrasp() {
	try {
		Graph original_graph(v, e);
		original_graph.readWeight(); // read all the weight and put into the weight vector
		original_graph.readEdges(); // read all the edges and put into the adjList
		original_graph.sort(); 
		original_graph.initializeHIndex();
		original_graph.initializeAccumulatedSum();

		Graph graph = original_graph; // copy original_graph to graph
		//graph.showGraphInformations(); // show all the informations of the input Graph
	
		BipartiteSolution best_s(&graph, partitionA_size, partitionB_size);
		
		// initialize C++ timer
        time_point<system_clock> start, execution_start, end;
		duration<double> elapsed_seconds;

		// tabu execution variables
		int local_search_iter, iter;
        int best_weight, local_weight, local_search_best_weight, best_solution = -1, avarage_solution = 0;
		int partitionCode; 

		srand (time(NULL));

        // start timing
        start = system_clock::now();
		
		cout << "Starting the algorithm...\n" << endl;

        while(loop--) {
			cout << abs(10 - loop) << " execution" << endl;

			BipartiteSolution s(&graph, partitionA_size, partitionB_size); // initialize all the variables and structures for solution

	        // start the first greedy random solution
			s.randomConstructive();
			if(s.checkBicliqueSize() == false) s.balanceBiclique();
			
			iter = 0;
			best_weight = s.getTotalWeight();

			cout << "Initial Solution: " << s.getTotalWeight() << endl;
			
			BipartiteSolution next_s(&graph, partitionA_size, partitionB_size);
			if(tabuLimitIteration != 0) {
				next_s.setTabuLimitIteration(tabuLimitIteration); 
			}

			// starting random generator for discrete distribution
			random_device device;
			mt19937 generator(device());

			// initializing vertex removed vector
			vertexInGraph.resize(v);
			fill(vertexInGraph.begin(), vertexInGraph.end(), true); // initialize the vector with all values set to true meaning that all the vertex are in the graph

			// start execution timing
			execution_start = system_clock::now();
			elapsed_seconds = system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			next_s.randomConstructive(); // starts a new optimal solution
			if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();
			local_weight = next_s.getTotalWeight();

			while((iter < total_iterations) || (execution_time <= execution_time_limit)) { // run ILS iterations	
				// local search 

				// first neighbor structure
				while(next_s.addPairOfVertices(iter));

				// second neighbor structure
				// it uses tabu list
				local_search_iter = iterWithoutImprovement;
				local_search_best_weight = local_weight;

				while(local_search_iter--) {
					partitionCode = random() % 2;
					next_s.swap1_1(partitionCode, iter);

					if(local_search_best_weight < next_s.getTotalWeight()) {
						local_search_best_weight = next_s.getTotalWeight(); 
						local_search_iter = iterWithoutImprovement;
					}
				}

				// third neighbor structure
				// it uses tabu list
				local_search_iter = iterWithoutImprovement;
				
				while(local_search_iter--) {
					partitionCode = random() % 2;
					next_s.swap1_1(partitionCode, iter);

					if(local_search_best_weight < next_s.getTotalWeight()) {
						local_search_best_weight = next_s.getTotalWeight(); 
						local_search_iter = iterWithoutImprovement;
					}
				}
				// end of local search

				if(next_s.checkBicliqueSize() == false) { next_s.balanceBiclique(); }
				local_weight = next_s.getTotalWeight();

				if(local_weight > best_weight) {
					s = next_s; // update best local solution
					best_weight = local_weight; // update best_weight
						
					elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count(); // update time_to_best 
					cout << "New best found: " << best_weight << endl;
					
					// TODO OLHAR LINHA ABAIXO - REDUCAO COM PROBLEMAS
					//next_s.reduceGraph(vertexInGraph, best_weight, -1, timeLimit); // start graph reduction
				}

				next_s.perturb();
				
				elapsed_seconds = system_clock::now() - execution_start; 
				execution_time = elapsed_seconds.count();
				iter++;

				assert(next_s.checkIntegrity());
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

			// execution informations
			elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

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

	        // start the random solution
			s.randomConstructive();
			if(s.checkBicliqueSize() == false) s.balanceBiclique();
			
			int best_weight = s.getTotalWeight(), local_weight, iter = 0;
			cout << "Initial Solution: " << s.getTotalWeight() << endl;
			
			NonBipartiteSolution next_s(&graph, partitionA_size, partitionB_size);

			// starting random generator for discrete distribution
			random_device device;
			mt19937 generator(device());

			// initializing vertex removed vector
			vertexInGraph.resize(v);
			fill(vertexInGraph.begin(), vertexInGraph.end(), true); // initialize the vector with all values set to true meaning that all the vertex are in the graph

			// start execution timing
			execution_start = system_clock::now();
			elapsed_seconds = system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			while((iter < total_iterations) || (execution_time <= execution_time_limit)) { // run ILS iterations
				next_s.randomConstructive(); // starts a new optimal solution
				if(next_s.checkBicliqueSize() == false) next_s.balanceBiclique();
				
				// TODO VER ABAIXO
				//next_s.VND(neighborCode);
				local_weight = next_s.getTotalWeight();

				if(local_weight > best_weight) {
					s = next_s; // update best local solution
					best_weight = local_weight; // update best_weight

					elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count(); // update time_to_best 
					cout << "New best found: " << best_weight << endl;
					// TODO create time limit mechanism (just made in the bipartite solution)
					// TODO add target conditional
					//next_s.reduceGraph(vertexInGraph, best_weight, -1, timeLimit); // start graph reduction
				}

				if(total_iterations == 0) {
					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();
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
			elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

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
			cout << "-  " << "-iWI + iterWithoutImprovement for local search stop criteria" << endl;
			cout << "-  " << "-tLI + tabuLimitIteration to set quantity of iterations for a tabu movement" << endl;
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
			else if(strcmp(argv[iter], "-G") == 0) target = atoi(argv[++iter]); 
			else if(strcmp(argv[iter], "-iWI") == 0) iterWithoutImprovement = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-tLI") == 0) tabuLimitIteration = atoi(argv[++iter]);
			else if(strcmp(argv[iter], "-eTc") == 0) edgeTimeConst = stod(argv[++iter]);
			else {
				cout << "Parameters Failed." << endl << "Type -help to use the parameters correctly" << endl;
				return 0;
			}
		}
	
		//cout << "Total iterations = " << total_iterations << endl;
		//cout << "Execution Time Limit = " << execution_time_limit << endl;
		
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