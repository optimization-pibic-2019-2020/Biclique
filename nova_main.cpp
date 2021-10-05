void ReactiveGrasp(int input_type) {
	try {
		K = 2; //TODO AQUI ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		Graph original_graph(v, e);
		original_graph.readWeight(); // read all the weight and put into the weight vector
		original_graph.readEdges(); // read all the edges and put into the adjList
		original_graph.sort(); 
		original_graph.initializeHIndex();
		original_graph.initializeAccumulatedSum();

		Graph graph = original_graph; // copy original_graph to graph
		//graph.showGraphInformations(); // show all the informations of the input Graph


		// Initializing all the solution objects

		NonBipartiteSolution nonBipartiteBestSolution(&graph, partitionA_size, partitionB_size);	
		NonBipartiteSolution nonBipartiteSolution(&graph, partitionA_size, partitionB_size);	
		NonBipartiteSolution nonBipartiteNextSolution(&graph, partitionA_size, partitionB_size);	

		BipartiteSolution bipartiteBestSolution(&graph, partitionA_size, partitionB_size);
		BipartiteSolution bipartiteSolution(&graph, partitionA_size, partitionB_size);
		BipartiteSolution bipartiteNextSolution(&graph, partitionA_size, partitionB_size);
		
		// initialize C++ timer
        time_point<system_clock> start, execution_start, end;
		duration<double> elapsed_seconds;

		// starting some variables related to the parameters and solution
        int best_solution = -1, avarage_solution = 0, best_weight, local_weight;
		int x, y, iter;
		double solutionAvarage;

		// starting random generator for discrete distribution
		random_device device;
		mt19937 generator(device());

        // start timing
        start = system_clock::now();
		
		cout << "Starting the algorithm...\n" << endl;

        while(loop--) {
			cout << abs(10 - loop) << " execution" << endl;

			// start the first greedy random solution
			if(input_type == 1) {
				NonBipartiteSolution nonBipartiteSolution(&graph, partitionA_size, partitionB_size);
				best_weight = nonBipartiteSolution.getTotalWeight();
				
				nonBipartiteSolution.greedyRandomizedConstructive(0.0);
				if(nonBipartiteSolution.checkBicliqueSize() == false) nonBipartiteSolution.balanceBiclique();
				
				cout << "Initial Solution: " << nonBipartiteSolution.getTotalWeight() << endl;
				
				NonBipartiteSolution nonBipartiteNextSolution(&graph, partitionA_size, partitionB_size);
			}
			else {	
				BipartiteSolution bipartiteSolution(&graph, partitionA_size, partitionB_size);
				best_weight = bipartiteSolution.getTotalWeight();
				
				bipartiteSolution.greedyRandomizedConstructive(0.0);
				if(bipartiteSolution.checkBicliqueSize() == false) bipartiteSolution.balanceBiclique();
				
				cout << "Initial Solution: " << bipartiteSolution.getTotalWeight() << endl;
				
				BipartiteSolution bipartiteNextSolution(&graph, partitionA_size, partitionB_size);
				verticesRemoved = bipartiteNextSolution.getRemovedVertices();
				edgesRemoved = bipartiteNextSolution.getRemovedEdges();
				cout << verticesRemoved << " " << edgesRemoved << endl;
			}

			x = beta, y = alpha_calibration, iter = 0;

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

				if(input_type == 1) {
					nonBipartiteNextSolution.greedyRandomizedConstructive(alphas[alpha_chosen]); // starts a new optimal solution
					if(nonBipartiteNextSolution.checkBicliqueSize() == false) nonBipartiteNextSolution.balanceBiclique();
					
					nonBipartiteNextSolution.VND(K);
					local_weight = nonBipartiteNextSolution.getTotalWeight();
				}
				else {
					bipartiteNextSolution.greedyRandomizedConstructive(alphas[alpha_chosen]); // starts a new optimal solution
					if(bipartiteNextSolution.checkBicliqueSize() == false) bipartiteNextSolution.balanceBiclique();
					
					bipartiteNextSolution.VND(K);
					local_weight = bipartiteNextSolution.getTotalWeight();
				}
				
				alphaSolutions[alpha_chosen] += local_weight;

				if(local_weight > best_weight) {
					// update best local solution
					if(input_type == 1) {
						nonBipartiteSolution = nonBipartiteNextSolution;
						nonBipartiteNextSolution.reduceGraph(vertexInGraph, best_weight); // start graph reduction
					}
					else {	
						bipartiteSolution = bipartiteNextSolution;
						bipartiteNextSolution.reduceGraph(vertexInGraph, best_weight); // start graph reduction
					}
					
					best_weight = local_weight; // update best_weight
					x = beta; // reinitialize parameter x

					elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
					time_to_best = elapsed_seconds.count(); // update time_to_best 
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

				if(total_iterations == 0) {
					elapsed_seconds = system_clock::now() - execution_start; 
					execution_time = elapsed_seconds.count();
				}
				else {
					iter++;
				}

				if(input_type == 1) {
					nonBipartiteNextSolution.restartSolution(vertexInGraph);
					assert(nonBipartiteNextSolution.checkIntegrity());
					assert(nonBipartiteNextSolution.checkMu());
				}
				else {
					bipartiteNextSolution.restartSolution(vertexInGraph);
					assert(bipartiteNextSolution.checkIntegrity());
					assert(bipartiteNextSolution.checkMu());
				}
				
			}

			// restarting graph
			graph = original_graph;

			if(input_type == 1) {
				if(nonBipartiteSolution.checkBicliqueSize() == false) {
					nonBipartiteSolution.balanceBiclique();
					best_weight = nonBipartiteSolution.getTotalWeight();
				}

				avarage_solution += best_weight;
				if(best_solution < best_weight) {
					nonBipartiteBestSolution = nonBipartiteSolution;
					best_solution = best_weight;
				}

				verticesRemoved = nonBipartiteSolution.getRemovedVertices();
				edgesRemoved = nonBipartiteSolution.getRemovedEdges();

				assert(nonBipartiteBestSolution.checkIntegrity());
				assert(nonBipartiteBestSolution.checkMu());
			}
			else {
				if(bipartiteSolution.checkBicliqueSize() == false) {
					bipartiteSolution.balanceBiclique();
					best_weight = bipartiteSolution.getTotalWeight();
				}

				avarage_solution += best_weight;
				if(best_solution < best_weight) {
					bipartiteBestSolution = bipartiteSolution;
					best_solution = best_weight;
				}

				verticesRemoved = bipartiteSolution.getRemovedVertices();
				edgesRemoved = bipartiteSolution.getRemovedEdges();

				assert(bipartiteBestSolution.checkIntegrity());
				assert(bipartiteBestSolution.checkMu());
			}

			// execution informations
			discrete_distribution<int> distribution(alphaProbability.begin(), alphaProbability.end());
			elapsed_seconds = std::chrono::system_clock::now() - execution_start; 
			execution_time = elapsed_seconds.count();

			cout << "\nAlpha probabilities (Reactive grasp)" << endl;
			for(int i = 0; i < 11; i++) {
				cout << "Alpha: " << ((double) i) / 10.0 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
			}

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
		//best_s.printSolution(); TODO PRINTAR SOLUÇAÕ 1!___________________________________________________________________________________________________________________________
		cout << "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" << endl;		
	}
	catch (std::exception &e) {
		cerr << e.what();
	}
}