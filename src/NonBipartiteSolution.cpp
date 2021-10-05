#include <iostream>
#include <vector>
#include <algorithm> 
#include <random>
#include <ctime>
#include "NonBipartiteSolution.hpp"
#define NDEBUG
#include <assert.h>

using namespace std;

// changed the variable name to avoid double definition error while compiling
random_device NonBipartiteDevice;
mt19937 NonBipartiteGenerator(NonBipartiteDevice());

NonBipartiteSolution::NonBipartiteSolution(Graph *graph, int partitionA_size, int partitionB_size):Solution() { // initialize all the variables and vectors
	this->graph = graph;
	int idx_weight, V = graph->getV();

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

	mu_A.resize(V);
	mu_B.resize(V);

	total_weight = 0;	
	removed_edges = 0;

    removed_size_A = 0; 
    removed_size_B = 0;

    for(int idx = 0; idx < V; idx++) {
        position_A[idx] = idx;
        position_B[idx] = idx;

        tightness_A[idx] = 0;
        tightness_B[idx] = 0;

        solution_A[idx] = idx;
        solution_B[idx] = idx;

        idx_weight = graph->get_weight(idx); 
        mu_A[idx] = idx_weight;
        mu_B[idx] = idx_weight;
    }
}

// return how many vertices were removed
int NonBipartiteSolution::getRemovedVertices() {
	return removed_size_A;
}

// update the free partiton of solution A and solution B
void NonBipartiteSolution::checkFreePartition() { 
	int vertex, u;

	for(vertex = solution_size_A; vertex < solution_size_A + free_size_A; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] > 0 || tightness_B[u] != solution_size_B) {
			moveFreeToNonFreePartition(u, 0);  
			vertex--;
		}
	}
	for(vertex = solution_size_B; vertex < solution_size_B + free_size_B; vertex++) { // verify the free partition of solution B
		u = solution_B[vertex];
		
		if(tightness_A[u] != solution_size_A || tightness_B[u] > 0) {
			moveFreeToNonFreePartition(u, 1);  
			vertex--;
		}
	}
	
}

// update the non-free partiton of solution A and solution B
void NonBipartiteSolution::checkNonFreePartition() { 
	int vertex, u, V = graph->getV();
	
	for(vertex = solution_size_A + free_size_A; vertex < V - removed_size_A; vertex++) { // verify the non-free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] == 0 && tightness_B[u] == solution_size_B) {
			moveNonFreeToFreePartition(u, 0);
			if(vertex >= solution_size_A + free_size_A) vertex--;  
		}
	}
	for(vertex = solution_size_B + free_size_B; vertex < V - removed_size_B; vertex++) { // verify the non-free partition of solution B
		u = solution_B[vertex];
		if(tightness_B[u] == 0 && tightness_A[u] == solution_size_A) {
			moveNonFreeToFreePartition(u, 1);  
			if(vertex >= solution_size_B + free_size_B) vertex--; 
		}
	}
}

// remove a vertex from solution A or B
void NonBipartiteSolution::removeVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	int weight_u = graph->get_weight(u);
	total_weight -= weight_u;

	moveSolutionToFreePartition(u, code);

	vector<int> &neighbors = graph->get_vertex_adjList(u);

	if(code == 0) {
		for(int neighbor : neighbors) {
			tightness_A[neighbor]--;
			mu_A[neighbor] += weight_u;

			// if the neighbor becomes free
			if(tightness_A[neighbor] == 0 && tightness_B[neighbor] == solution_size_B) {
				moveNonFreeToFreePartition(neighbor, code);
			}
		}
	}
	else {
		for(int neighbor : neighbors) {
			tightness_B[neighbor]--;
			mu_B[neighbor] += weight_u;

			// if the neighbor becomes free
			if(tightness_B[neighbor] == 0 && tightness_A[neighbor] == solution_size_A) {
				moveNonFreeToFreePartition(neighbor, code);
			}
		}
	}
}

// add a vertex to solution A or B
void NonBipartiteSolution::addVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B 
	int weight_u = graph->get_weight(u);
	total_weight += weight_u;

	moveFreeToSolutionPartition(u, code);

	vector<int> &neighbors = graph->get_vertex_adjList(u);

	if(code == 0) { // placing vertex u in solution A
		for(int neighbor : neighbors) {
			// increase the tighness of each neighbor by one
			tightness_A[neighbor]++;

			// update the value of mu from each neighbor in solution A
			mu_A[neighbor] -= weight_u;
		}
	}
	else { // placing vertex u in solution B
		for(int neighbor : neighbors) {
			// increase the tighness of each neighbor by one
			tightness_B[neighbor]++;

			// update the value of mu from each neighbor in solution B
			mu_B[neighbor] -= weight_u;
		}
	}
}

bool NonBipartiteSolution::checkMu() { // check if mu_A and mu_B are correct
	int V = graph->getV(), idx_weight;
	int mu_A_temp[V], mu_B_temp[V]; // initialize two mu_ arrays to recalculate from zero each vertex 

	for(int idx = 0; idx < V; idx++) {
		idx_weight = graph->get_weight(idx);
		mu_A_temp[idx] = idx_weight;
		mu_B_temp[idx] = idx_weight;
	}

	for(int idx = 0; idx < solution_size_A; idx++) { // calculate mu_A_temp
		int u = solution_A[idx];
		int weight_u = graph->get_weight(u);
		vector<int> &neighbors = graph->get_vertex_adjList(u);

		for(int neighbor : neighbors) {
			// update the value of mu_A_temp from each neighbor in solution A
			mu_A_temp[neighbor] -= weight_u;
		}
	}

	for(int idx = 0; idx < solution_size_B; idx++) { // calculate mu_B_temp
		int u = solution_B[idx];
		int weight_u = graph->get_weight(u);
		vector<int> &neighbors = graph->get_vertex_adjList(u);

		for(int neighbor : neighbors) {
			// update the value of mu_B_temp from each neighbor in solution B
			mu_B_temp[neighbor] -= weight_u;
		}
	}

	for(int idx = 0; idx < V; idx++) { // check the integrity of each mu_
		if(mu_A[idx] != mu_A_temp[idx] || mu_B[idx] != mu_B_temp[idx]) {
			cout << "Mu error" << endl;
			return false;
		}
	}

	return true;
}

void NonBipartiteSolution::restartSolution(vector<bool> &vertexInGraph) {
	this->graph = graph;
	int vertex, vertex_weight, V = graph->getV();

	solution_size_A = 0;
	solution_size_B = 0;

	free_size_A = V - removed_size_A;
	free_size_B = V - removed_size_B;

	total_weight = 0;	

	
    for(int idx = 0; idx < V - removed_size_A; idx++) { // removed_size_A is equal to removed_size_B in this type of instance
        vertex = solution_A[idx];
        vertex_weight = graph->get_weight(vertex);
        mu_A[vertex] = vertex_weight;
        mu_B[vertex] = vertex_weight;

        tightness_A[vertex] = 0;
        tightness_B[vertex] = 0;		
    }
}

void NonBipartiteSolution::swapVertices(int vertex, int code) { // code == 0 for partition A and code != 0 for partition B
	vector<int> &neighbors = graph->get_vertex_adjList(vertex);

	if(code == 0) {
		for(int neighbor : neighbors) {
			if(position_A[neighbor] < solution_size_A) {
				// move the neighbor to the nonFreePartition A
				removeVertex(neighbor, code);
				moveFreeToNonFreePartition(neighbor, code);
				// move the vertex to the solution A
				addVertex(vertex, code);
				break;
			}
		}
	}
	else {
		for(int neighbor : neighbors) {
			if(position_B[neighbor] < solution_size_B) {
				// move the neighbor to the nonFreePartition B
				removeVertex(neighbor, code);
				moveFreeToNonFreePartition(neighbor, code);
				// move the vertex to the solution B
				addVertex(vertex, code);
				break;
			}
		}
	}

	checkFreePartition();
	checkNonFreePartition();
}

bool NonBipartiteSolution::swap1_1(int code) { // code == 0 for partition A and code != 0 for partition B
	int V = graph->getV(), vertex, best_vertex, improvement = 0; // improvement means how much weight the partition will get after the swap
	if(code == 0) {
		for(int iter = solution_size_A + free_size_A; iter < V - removed_size_A; iter++) {
			vertex = solution_A[iter];
			if(tightness_A[vertex] == 1 && tightness_B[vertex] == solution_size_B && mu_A[vertex] > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = mu_A[vertex];
			}
		}
	}
	else {
		for(int iter = solution_size_B + free_size_B; iter < V - removed_size_B; iter++) {
			vertex = solution_B[iter];
			if(tightness_B[vertex] == 1 && tightness_A[vertex] == solution_size_A && mu_B[vertex] > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = mu_B[vertex];
			}
		}
	}

	if(improvement > 0) {
		swapVertices(best_vertex, code);
		return true;
	}

	return false;
}

bool NonBipartiteSolution::swap2_2(int code) { // code == 0 for partition A and code != 0 for partition B
	int V = graph->getV(), vertex1, vertex2; // improvement means how much weight the partition will get after the swap
	if(code == 0 && solution_size_A >= 2) {
		for(int iter = solution_size_A + free_size_A; iter < V - removed_size_A; iter++) {
			vertex1 = solution_A[iter];
			if(tightness_A[vertex1] == 1 && tightness_B[vertex1] == solution_size_B) { // possible candidate for swap2_2
				for(int iter2 = iter + 1; iter2 < V - removed_size_A; iter2++) {
					vertex2 = solution_A[iter2];
					if(tightness_A[vertex2] == 1 && tightness_B[vertex2] == solution_size_B) { // another possible candidate for swap2_2
						if((mu_A[vertex1] + mu_A[vertex2]) > 0 && !isNeighbor(vertex1, vertex2) && !sameNeighbor(vertex1, vertex2, code)) { // verify if they are neighbors and if they share the same neighbor in the solution
							// do the swap2_2 for partition A
							swapVertices(vertex1, code);
							swapVertices(vertex2, code);
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
							swapVertices(vertex1, code);
							swapVertices(vertex2, code);
							return true;	
						}
					}	
				}
			}
		}
	}

	return false;
}

// try to add a pair of vertices to the solution (one to partition A and one to partition B)
bool NonBipartiteSolution::addPairOfVertices() { 
	if(free_size_A > 0 && free_size_B > 0) {
		uniform_int_distribution<int> distribution(0, free_size_A - 1);
		int iter = solution_size_A + distribution(NonBipartiteGenerator);
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

void NonBipartiteSolution::VND(int K) { // run VND iterations
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
				if(swap1_1(code)) k = 1;
				else if(swap1_1(abs(code - 1))) k = 1; // if first partition goes wrong it will try the next partition
				else k++;
				break;
			case 3:
				code = rand() % 2; // try to randomly choose a partition
				if(swap2_2(code)) k = 1;			
				else if(swap2_2(abs(code - 1))) k = 1; // if first partition goes wrong it will try the next partition
				else k++;
				break;
			default:
				k++;
				break;
		}
	}
}

// construct the restricted candidate list for a specific partition and choose a random vertex from it to put in the solution
void NonBipartiteSolution::rclConstruction(int code, double alpha) { 
	int iter, vertex, vertex_weight;
	int c_min = 100000000, c_max = -1; // variables to get the interval for the quality-based RCL list
	vector<int> &weight = graph->get_weight_list();
	rclList.clear();
	rclListProbability.clear();

	if(code == 0 && free_size_A != 0) { // partition A
		for(iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // gets the c_min and c_max parameters from the free vertices
			vertex = solution_A[iter];
			vertex_weight = weight[vertex];

			if(vertex_weight <= c_min) c_min = vertex_weight;
			if(vertex_weight >= c_max) c_max = vertex_weight;
		}

		for(iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // creates the rcl List based in a quality-based construction
			vertex = solution_A[iter];
			vertex_weight = weight[vertex];

			if((c_min + alpha * (c_max - c_min)) <= ((double) vertex_weight) && vertex_weight <= c_max) { // condition to get a good quality in the RCL list
				rclList.push_back(vertex);
			}
		}

		// creates the rclListProbability
		createRclProbability();
		
		// generate a random number between [0, rclList.size() - 1] using a linear bias function
		discrete_distribution<int> distribution(rclListProbability.begin(), rclListProbability.end());
		int pos = distribution(NonBipartiteGenerator);
		vertex = rclList[pos]; 
		addVertex(vertex, 0);
	}
	else if(code != 0 && free_size_B != 0) { // partition B
		for(iter = solution_size_B; iter < solution_size_B + free_size_B; iter++) { // gets the c_min and c_max parameters from the free vertices
			vertex = solution_B[iter];
			vertex_weight = weight[vertex];

			if(vertex_weight <= c_min) c_min = vertex_weight;
			if(vertex_weight >= c_max) c_max = vertex_weight;
		}

		for(iter = solution_size_B; iter < solution_size_B + free_size_B; iter++) { // creates the rcl List based in a quality-based construction
			vertex = solution_B[iter];
			vertex_weight = weight[vertex];

			if((c_min + alpha * (c_max - c_min)) <= ((double) vertex_weight) && vertex_weight <= c_max) { // condition to get a good quality in the RCL list
				rclList.push_back(vertex);
			}
		}
	
		// creates the rclListProbability
		createRclProbability();

		// generate a random number between [0, rclList.size() - 1] using a linear bias function
		discrete_distribution<int> distribution(rclListProbability.begin(), rclListProbability.end());
		int pos = distribution(NonBipartiteGenerator);
		vertex = rclList[pos]; 
		addVertex(vertex, 1);
	}

	checkFreePartition();
}

void NonBipartiteSolution::greedyRandomizedConstructive(double p) { 
	while(free_size_A > 0 && free_size_B > 0) { // can generate a solution with an extra vertex in solution A
		rclConstruction(0, p); 
		rclConstruction(1, p); 
		assert(checkIntegrity());
		assert(checkMu());
	}
}

void NonBipartiteSolution::reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight) {
	int vertex, bicliquePredictedWeight, verticesRemoved = 0;

	if(minBicliqueWeight > bestWeight) { return; } // avoid a final loop if there is no possible vertex to be removed

    for(unsigned int iter = 0; iter < solution_A.size() - removed_size_A; iter++) { // examine all the vertices that are in the solution (for non bipartite instances, removed_size_A == removed_size_B)
        vertex = solution_A[iter];
        bicliquePredictedWeight = predictBicliqueWeight(vertex);

		if(bicliquePredictedWeight < minBicliqueWeight) { minBicliqueWeight = bicliquePredictedWeight; }

        if(bicliquePredictedWeight <= bestWeight) { // removes u from graph
            removeVertexFromGraph(vertex);
            vertexInGraph[vertex] = false; 
            moveVertexToRemovedVertices(vertex, 0);
            moveVertexToRemovedVertices(vertex, 1);
            verticesRemoved++;
        }		
    }

	if(verticesRemoved > 0) reduceGraph(vertexInGraph, bestWeight, minBicliqueWeight);
}

void NonBipartiteSolution::balanceBiclique() { // remove the vertex with the worst weight in solution A to balance the Biclique
	int minimum_weight = solution_A[0], vertex_to_remove = solution_A[0], actual_vertex;
	
	for(int iter = 1; iter < solution_size_A; iter++) {
		actual_vertex = solution_A[iter];
		if(mu_A[actual_vertex] < minimum_weight) {
			vertex_to_remove = actual_vertex;
			minimum_weight = mu_A[actual_vertex];
		}
	}

	removeVertex(vertex_to_remove, 0); // remove the worst vertex in solution A	
	checkNonFreePartition();
}