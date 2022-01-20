#include <iostream>
#include <vector>
#include <queue>
#include <algorithm> 
#include <random>
#include <ctime>
#include <iomanip>
#include "BipartiteSolution.hpp"
#define NDEBUG
#include <assert.h>
#include <chrono>

using namespace std::chrono;
using namespace std;

random_device BipartiteDevice;
mt19937 BipartiteGenerator(BipartiteDevice());

BipartiteSolution::BipartiteSolution(Graph *graph, int partitionA_size, int partitionB_size) { // initialize all the variables and vectors
	this->graph = graph;
	int V = graph->getV();

    partition_size_A = partitionA_size;
    partition_size_B = partitionB_size;

	solution_A.resize(V);
	solution_B.resize(V);

	free_size_A = partitionA_size;
	free_size_B = partitionB_size;

	position_A.resize(V);
	position_B.resize(V);

	solution_size_A = 0;
	solution_size_B = 0;

	tightness_A.resize(V);
	tightness_B.resize(V);

	total_weight = 0;	
	removed_edges = 0;

    // initializing bipartitioned instance
    removed_size_A = free_size_B; // move all the vertices from partition B to the removed part of partition A  
    removed_size_B = free_size_A; // move all the vertices from partition A to the removed part of partition B
    
    // partition A
    for(int idx = 0; idx < free_size_A; idx++) {
        tightness_A[idx] = 0;
        tightness_B[idx] = 0;

        position_A[idx] = idx;
        position_B[idx] = idx + free_size_B; // move a vertex from partition A to the removed vertices of solution B

        solution_A[idx] = idx;
        solution_B[idx + free_size_B] = idx;
    }
    
    // partition B
    for(int idx = free_size_A; idx < V; idx++) {
  	    tightness_A[idx] = 0;
        tightness_B[idx] = 0;

        position_A[idx] = idx; // move a vertex from partition B to the removed vertices of solution A
        position_B[idx] = idx - free_size_A;

        solution_A[idx] = idx;
        solution_B[idx - free_size_A] = idx;
    }
}

// return how many vertices were removed
int BipartiteSolution::getRemovedVertices() {
	return (removed_size_A - partition_size_B) + (removed_size_B - partition_size_A);
}

// update the free partiton of solution A and solution B
void BipartiteSolution::checkFreePartition() { 
	int vertex, u;

	for(vertex = solution_size_A; vertex < solution_size_A + free_size_A; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_B[u] != solution_size_B) {
			moveFreeToNonFreePartition(u, 0);  
			vertex--;
		}
	}
	for(vertex = solution_size_B; vertex < solution_size_B + free_size_B; vertex++) { // verify the free partition of solution B
		u = solution_B[vertex];
		
		if(tightness_A[u] != solution_size_A) {
			moveFreeToNonFreePartition(u, 1);  
			vertex--;
		}
	}
	
}

// update the non-free partiton of solution A and solution B
void BipartiteSolution::checkNonFreePartition() { 
	int vertex, u, V = graph->getV();
	
	for(vertex = solution_size_A + free_size_A; vertex < V - removed_size_A; vertex++) { // verify the non-free partition of solution A
		u = solution_A[vertex];
		if(tightness_B[u] == solution_size_B) {
			moveNonFreeToFreePartition(u, 0);
			if(vertex >= solution_size_A + free_size_A) vertex--;  
		}
	}
	for(vertex = solution_size_B + free_size_B; vertex < V - removed_size_B; vertex++) { // verify the non-free partition of solution B
		u = solution_B[vertex];
		if(tightness_A[u] == solution_size_A) {
			moveNonFreeToFreePartition(u, 1);  
			if(vertex >= solution_size_B + free_size_B) vertex--; 
		}
	}
}

// remove a vertex from solution A or B
void BipartiteSolution::removeVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	int weight_u = graph->get_weight(u);
	total_weight -= weight_u;

	moveSolutionToFreePartition(u, code);

	vector<int> &neighbors = graph->get_vertex_adjList(u);

    // update tightness
	if(code == 0) {
		for(int neighbor : neighbors) { tightness_A[neighbor]--;}
	}
	else {
        for(int neighbor : neighbors) { tightness_B[neighbor]--;}
	}
}

// add a vertex to solution A or B
void BipartiteSolution::addVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B 
	int weight_u = graph->get_weight(u);
	total_weight += weight_u;

	moveFreeToSolutionPartition(u, code);

	vector<int> &neighbors = graph->get_vertex_adjList(u);

    // update tightness
	if(code == 0) { // placing vertex u in solution A
		for(int neighbor : neighbors) { tightness_A[neighbor]++; }
	}
	else { // placing vertex u in solution B
		for(int neighbor : neighbors) { tightness_B[neighbor]++; }
	}
}

void BipartiteSolution::restartSolution(vector<bool> &vertexInGraph) {
	this->graph = graph;
	int vertex, V = graph->getV();

	solution_size_A = 0;
	solution_size_B = 0;

	free_size_A = V - removed_size_A;
	free_size_B = V - removed_size_B;

	total_weight = 0;	

	// restarting bipartite instance
    for(int idx = 0; idx < V - removed_size_A; idx++) {
        vertex = solution_A[idx];

        tightness_A[vertex] = 0;
        tightness_B[vertex] = 0;		
    }

    for(int idx = 0; idx < V - removed_size_B; idx++) {
        vertex = solution_B[idx];

        tightness_A[vertex] = 0;
        tightness_B[vertex] = 0;		
    }
}

void BipartiteSolution::swapVertices(int vertex_to_add, int vertex_to_remove, int code) { // code == 0 for partition A and code != 0 for partition B	
    // move the vertex_to_remove to the Free Partition
    removeVertex(vertex_to_remove, code);
    // move the vertex_to_add to the solution Partition
    addVertex(vertex_to_add, code);

	
	checkFreePartition();
	checkNonFreePartition();
}

bool BipartiteSolution::swap1_1(int code) { // code == 0 for partition A and code != 0 for partition B
    // improvement is just a variable to get the best neighbor_size (in the free partition) 
    // and the worse neighbor_size (in the solution partition)
	int neighbor_size, improvement = 0;  
    int vertex, vertex_to_remove = -1, best_vertex = -1;

	if(code == 0) {
        for(int iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // initial loop to get the best tightness among the free vertices
			vertex = solution_A[iter];
            neighbor_size = graph->get_vertex_adjList(vertex).size();
			if(neighbor_size > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = neighbor_size;
			}
		}

        for(int iter = 0; iter < solution_size_A; iter++) { // last loop to get the worse tightness of the solution vertex part
            vertex = solution_A[iter];
            neighbor_size = graph->get_vertex_adjList(vertex).size();
            if(neighbor_size < improvement) {
                vertex_to_remove = vertex;
                improvement = neighbor_size;
            }
        }		
	}
	else {
        for(int iter = solution_size_B; iter < solution_size_B + free_size_B; iter++) { // initial loop to get the best tightness among the free vertices
			vertex = solution_B[iter];
            neighbor_size = graph->get_vertex_adjList(vertex).size();
			if(neighbor_size > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = neighbor_size;
			}
		}

        for(int iter = 0; iter < solution_size_B; iter++) { // last loop to get the worse tightness of the solution vertex part
            vertex = solution_B[iter];
            neighbor_size = graph->get_vertex_adjList(vertex).size();
            if(neighbor_size < improvement) {
                vertex_to_remove = vertex;
                improvement = neighbor_size;
            }
        }	
	}

	if(vertex_to_remove > -1) {
        swapVertices(best_vertex, vertex_to_remove, code);
        return true;
    } 

	return false;
}

bool BipartiteSolution::swap1_k(int code) { // code == 0 for partition A and code != 0 for partition B
    int vertex, vertexChosen;
	bool operationSuccessFlag = false;
	vector<int> possibleVerticesToSolution;

	if(code == 0) {
		// initial loop to get the possible vertices to join the solution
        for(unsigned int iter = solution_size_A + free_size_A; iter < solution_A.size() - removed_size_A; iter++) { 
			vertex = solution_A[iter];
			if (tightness_B[vertex] == solution_size_B - 1) {
				possibleVerticesToSolution.push_back(vertex);
			}
		}

		if (possibleVerticesToSolution.size() == 0) { return operationSuccessFlag; }
		
		operationSuccessFlag = true;

		srand(time(0));
		vertexChosen = possibleVerticesToSolution[rand() % possibleVerticesToSolution.size()];

		for(int iter = 0; iter < solution_size_B; iter++) { 
			vertex = solution_B[iter];
			if (!isNeighbor(vertex, vertexChosen)) {
				//cout << "Removed vertex in B: " << vertex << endl; 
				removeVertex(vertex, 1);
			}
		}

		checkFreePartition();
		checkNonFreePartition();

		for(unsigned int iter = 0; iter < possibleVerticesToSolution.size(); iter++) { 
			vertex = possibleVerticesToSolution[iter];
			if (tightness_B[vertex] == solution_size_B && position_A[vertex] >= solution_size_A) {
				//cout << vertex << endl;
				addVertex(vertex, code);
			}
		}
	}
	else {
		// initial loop to get the possible vertices to join the solution
        for(unsigned int iter = solution_size_B + free_size_B; iter < solution_B.size() - removed_size_B; iter++) { 
			vertex = solution_B[iter];
			if (tightness_A[vertex] == solution_size_A - 1) {
				possibleVerticesToSolution.push_back(vertex);
			}
		}

		if (possibleVerticesToSolution.size() == 0) { return operationSuccessFlag; }
		
		operationSuccessFlag = true;

		srand(time(0));
		vertexChosen = possibleVerticesToSolution[rand() % possibleVerticesToSolution.size()];

		for(int iter = 0; iter < solution_size_A; iter++) { 
			vertex = solution_A[iter];
			if (!isNeighbor(vertex, vertexChosen)) {
				//cout << "Removed vertex in A: " << vertex << endl;
				removeVertex(vertex, 0);
			}
		}

		checkFreePartition();
		checkNonFreePartition();

		for(unsigned int iter = 0; iter < possibleVerticesToSolution.size(); iter++) { 
			vertex = possibleVerticesToSolution[iter];
			if (tightness_A[vertex] == solution_size_A && position_B[vertex] >= solution_size_B) {
				addVertex(vertex, code);
			}
		}
	}

	checkFreePartition();
	checkNonFreePartition();

	return operationSuccessFlag;
}

bool BipartiteSolution::swap2_2(int code) { // code == 0 for partition A and code != 0 for partition B
    /*
	int V = graph->getV(), vertex1, vertex2, improvement = 0; // improvement means how much weight the partition will get after the swap
	if(code == 0 && solution_size_A >= 2) {
		for(int iter = solution_size_A + free_size_A; iter < V - removed_size_A; iter++) {
			vertex1 = solution_A[iter];
			if(tightness_A[vertex1] == 1 && tightness_B[vertex1] == solution_size_B) { // possible candidate for swap2_2
				for(int iter2 = iter + 1; iter2 < V - removed_size_A; iter2++) {
					vertex2 = solution_A[iter2];
					if(tightness_A[vertex2] == 1 && tightness_B[vertex2] == solution_size_B) { // another possible candidate for swap2_2
						if((mu_A[vertex1] + mu_A[vertex2]) > 0 && !isNeighbor(vertex1, vertex2) && !sameNeighbor(vertex1, vertex2, code)) { // verify if they are neighbors and if they share the same neighbor in the solution
							// do the swap2_2 for partition A
							oneImprovement(vertex1, code);
							oneImprovement(vertex2, code);
							return true;	
						}
					}	
				}
			}
		}
	}
	else if(code != 0 && solution_size_B >= 2) {
		for(int iter = solution_size_B + free_size_B; iter < V - removed_size_B; iter++) {
			vertex1 = solution_B[iter];
			if(tightness_B[vertex1] == 1 && tightness_A[vertex1] == solution_size_A) { // possible candidate for swap2_2
				for(int iter2 = iter; iter2 < V - removed_size_B; iter2++) {
					vertex2 = solution_B[iter2];
					if(tightness_B[vertex2] == 1 && tightness_A[vertex2] == solution_size_A) { // another possible candidate for swap2_2
						if((mu_B[vertex1] + mu_B[vertex2]) > 0 && !isNeighbor(vertex1, vertex2) && !sameNeighbor(vertex1, vertex2, code)) { // verify if they are neighbors and if they share the same neighbor in the solution
							// do the swap2_2 for partition B
							oneImprovement(vertex1, code);
							oneImprovement(vertex2, code);
							return true;	
						}
					}	
				}
			}
		}
	}
    */
	return false;
}

bool BipartiteSolution::addPairOfVertices() { // verify if can add a pair of vertices to the solution (one to partition A and one to partition B)
	if(free_size_A > 0 && free_size_B > 0) {
		uniform_int_distribution<int> distribution(0, free_size_A - 1);
		int iter = solution_size_A + distribution(BipartiteGenerator);
		int u, vertex1 = -1, vertex2 = -1, loop = free_size_A;
		int vertex_weight, neighbor_weight;
		
		while(loop--) { // check if there is a pair to join the solution
			if(iter < free_size_A) u = solution_A[iter];
			else u = solution_A[(iter % free_size_A) + solution_size_A]; 

			vertex_weight = graph->get_weight(u);

			vector<int> &neighbors = graph->get_vertex_adjList(u);
			for(int neighbor : neighbors) {
				// check if there is an edge between these free vertices
				neighbor_weight = graph->get_weight(neighbor);
				if(position_B[neighbor] >= solution_size_B && position_B[neighbor] < solution_size_B + free_size_B && 0 <= vertex_weight + neighbor_weight) {
					vertex1 = u;
					vertex2 = neighbor;
					break;
				}
			}
			if(vertex1 != -1 && vertex2 != -1) {
				addVertex(vertex1, 0);
				addVertex(vertex2, 1);
				checkFreePartition();
				return true;
			}
			iter++;
		}
	}
	return false;
}

void BipartiteSolution::VND(int K) { // run VND iterations
	srand (time(NULL));
	int k = 1, code;
	while(k <= K) {
		switch(k) {
			case 1:
				if(addPairOfVertices()) k = 1;
				else k++;
				break;
			case 2:
				code = rand() % 2; // try to randomly choose a partition
				if(swap1_1(code)) { k = 1; }
				else if(swap1_1(abs(code - 1))) { k = 1; } // if first partition goes wrong it will try the next partition
				else { k++; }
				break;
			case 3:
				// TODO ANALYZE THE LINE BELOW
				//if (solution_size_A <= 1 || solution_size_B <= 1) { k++; break; }
				code = rand() % 2; // try to randomly choose a partition
				if (swap1_k(code)) { k = 1; }			
				else { k++; }
				//else if (swap1_k(abs(code - 1))) { k = 1; }// if first partition goes wrong it will try the next partition
				
				break;
			default:
				k++;
				break;
		}
	}
}

void BipartiteSolution::rclConstruction(int code, double alpha) { // construct the restricted candidate list for a specific partition and choose a random vertex from it to put in the solution
	int iter, vertex, vertexValue;
	// int vertex_weight;
	int c_min = 100000000, c_max = -1; // variables to get the interval for the quality-based RCL list
	//vector<int> &weight = graph->get_weight_list();
	vector<int> &am = graph->getAmList();
	rclList.clear();
	rclListProbability.clear();

	if(code == 0 && free_size_A != 0) { // partition A
		for(iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // gets the c_min and c_max parameters from the free vertices
			vertex = solution_A[iter];
			//vertex_weight = weight[vertex];
			vertexValue = am[vertex];


			if(vertexValue <= c_min) c_min = vertexValue;
			if(vertexValue >= c_max) c_max = vertexValue;
		}

		for(iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // creates the rcl List based in a quality-based construction
			vertex = solution_A[iter];
			//vertex_weight = weight[vertex];
			vertexValue = am[vertex];

			if((c_min + alpha * (c_max - c_min)) <= ((double) vertexValue) && vertexValue <= c_max) { // condition to get a good quality in the RCL list
				rclList.push_back(vertex);
			}
		}

		// creates the rclListProbability
		createRclProbability();
		
		// generate a random number between [0, rclList.size() - 1] using a linear bias function
		discrete_distribution<int> distribution(rclListProbability.begin(), rclListProbability.end());
		int pos = distribution(BipartiteGenerator);
		vertex = rclList[pos]; 
		addVertex(vertex, 0);
	}
	else if(code != 0 && free_size_B != 0) { // partition B
		for(iter = solution_size_B; iter < solution_size_B + free_size_B; iter++) { // gets the c_min and c_max parameters from the free vertices
			vertex = solution_B[iter];
			//vertex_weight = weight[vertex];
			vertexValue = am[vertex];

			if(vertexValue <= c_min) c_min = vertexValue;
			if(vertexValue >= c_max) c_max = vertexValue;
		}

		for(iter = solution_size_B; iter < solution_size_B + free_size_B; iter++) { // creates the rcl List based in a quality-based construction
			vertex = solution_B[iter];
			//vertex_weight = weight[vertex];
			vertexValue = am[vertex];

			if((c_min + alpha * (c_max - c_min)) <= ((double) vertexValue) && vertexValue <= c_max) { // condition to get a good quality in the RCL list
				rclList.push_back(vertex);
			}
		}
	
		// creates the rclListProbability
		createRclProbability();

		// generate a random number between [0, rclList.size() - 1] using a linear bias function
		discrete_distribution<int> distribution(rclListProbability.begin(), rclListProbability.end());
		int pos = distribution(BipartiteGenerator);
		vertex = rclList[pos]; 
		addVertex(vertex, 1);
	}

	checkFreePartition();
}

void BipartiteSolution::greedyRandomizedConstructive(double p) { 
	while(free_size_A > 0 && free_size_B > 0) { // can generate a solution with an extra vertex in solution A
		rclConstruction(0, p); 
		rclConstruction(1, p); 
		assert(checkIntegrity());
	}
}

void BipartiteSolution::reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight, double timeLimit) {
	int vertex, bicliquePredictedWeight, verticesRemoved = 0;
	double totalTime = 0;
	bool timeFlag = false;
	time_point<system_clock> start, end;
	duration<double> elapsedSeconds;
	
	if(minBicliqueWeight > bestWeight) { return; } // avoid a final loop if there is no possible vertex to be removed

	// start counting time 
	start = system_clock::now();

    for(unsigned int iter = 0; iter < solution_A.size() - removed_size_A; iter++) { // examine all the vertices that are in partition A
        vertex = solution_A[iter];
        bicliquePredictedWeight = predictBicliqueWeight(vertex);

		if(bicliquePredictedWeight < minBicliqueWeight) { minBicliqueWeight = bicliquePredictedWeight; }
		
        if(bicliquePredictedWeight <= bestWeight) { // removes vertex from graph 
            removeVertexFromGraph(vertex);
            moveVertexToRemovedVertices(vertex, 0);
            verticesRemoved++;
        }

		
		end = std::chrono::system_clock::now();
		elapsedSeconds = end - start;
		totalTime = elapsedSeconds.count();

		if(totalTime > timeLimit) {
			timeFlag = true;
			break;
		} 
    }

	// start counting time
	start = system_clock::now();
	totalTime = 0;    

    for(unsigned int iter = 0; iter < solution_B.size() - removed_size_B; iter++) { // examine all the vertices that are in partition B
        vertex = solution_B[iter];
        bicliquePredictedWeight = predictBicliqueWeight(vertex);

		if(bicliquePredictedWeight < minBicliqueWeight) { minBicliqueWeight = bicliquePredictedWeight; }

        if(bicliquePredictedWeight <= bestWeight) { // removes vertex from graph
            removeVertexFromGraph(vertex);
            moveVertexToRemovedVertices(vertex, 1);
            verticesRemoved++;
        }

		
		end = std::chrono::system_clock::now();
		elapsedSeconds = end - start;
		totalTime = elapsedSeconds.count();

		if(totalTime > timeLimit) {
			timeFlag = true;
			break;
		}
    }

	
	if(verticesRemoved > 0 && !timeFlag) reduceGraph(vertexInGraph, bestWeight, minBicliqueWeight, timeLimit);
}

void BipartiteSolution::balanceBiclique() { // remove the worst vertices in the biggest part to balance the Biclique
	int vertex, vertexWeight;
	int diff = abs(solution_size_A - solution_size_B);
	int code = (solution_size_A > solution_size_B) ? 0 : 1;
   	priority_queue<pair<int, int>> pq;

	if (code == 0) { 
		for (int iter = 0; iter < solution_size_A; iter++) {
			vertex = solution_A[iter];
			vertexWeight = graph->get_weight(vertex);
			pq.push(make_pair(vertexWeight, vertex));
		}
	}
	else {
		for (int iter = 0; iter < solution_size_B; iter++) {
			vertex = solution_B[iter];
			vertexWeight = graph->get_weight(vertex);
			pq.push(make_pair(vertexWeight, vertex));
		}
	}

	for (int iter = 0; iter < diff; iter++) {
		vertex = pq.top().second;
		removeVertex(vertex, code); 
		pq.pop();
	}
		
	checkFreePartition();
	checkNonFreePartition();
}