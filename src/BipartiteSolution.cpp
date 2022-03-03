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

	tabuList.resize(V);

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

		/* starts with a big negative number to allow 
		   every vertex to join or be removed 
		   from the solution on first iteration 
		*/

		tabuList[idx] = -1000; 
    }
    
    // partition B
    for(int idx = free_size_A; idx < V; idx++) {
  	    tightness_A[idx] = 0;
        tightness_B[idx] = 0;

        position_A[idx] = idx; // move a vertex from partition B to the removed vertices of solution A
        position_B[idx] = idx - free_size_A;

        solution_A[idx] = idx;
        solution_B[idx - free_size_A] = idx;

		/* starts with a big negative number to allow 
		   every vertex to join or be removed 
		   from the solution on first iteration 
		*/

		tabuList[idx] = -1000;
    }
}

// return how many vertices were removed
int BipartiteSolution::getRemovedVertices() {
	return (removed_size_A - partition_size_B) + (removed_size_B - partition_size_A);
}

// update the solution partiton of solution A and solution B
void BipartiteSolution::checkSolutionPartition() { 
	int vertex, u;

	for(vertex = 0; vertex < solution_size_A; vertex++) { // verify the solution partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] != 0 || tightness_B[u] != solution_size_B) {
			moveSolutionToFreePartition(u, 0);
			moveFreeToNonFreePartition(u, 0);  
			vertex--;
		}
	}
	for(vertex = 0; vertex < solution_size_B; vertex++) { // verify the solution partition of solution B
		u = solution_B[vertex];
		
		if(tightness_B[u] != 0 || tightness_A[u] != solution_size_A) {
			moveSolutionToFreePartition(u, 1);
			moveFreeToNonFreePartition(u, 1);  
			vertex--;
		}
	}
	
}

// update the free partiton of solution A and solution B
void BipartiteSolution::checkFreePartition() { 
	int vertex, u;

	for(vertex = solution_size_A; vertex < solution_size_A + free_size_A; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_B[u] != solution_size_B) {
			moveFreeToNonFreePartition(u, 0);  
			if(vertex >= solution_size_A) vertex--;
		}
	}
	for(vertex = solution_size_B; vertex < solution_size_B + free_size_B; vertex++) { // verify the free partition of solution B
		u = solution_B[vertex];
		
		if(tightness_A[u] != solution_size_A) {
			moveFreeToNonFreePartition(u, 1);  
			if(vertex >= solution_size_B) vertex--;
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

	checkNonFreePartition();
}

// add a vertex to solution A or B
void BipartiteSolution::addVertex(int u, int code, int iteration) { // code == 0 for solution A and code != 0 for solution B 
	int weight_u = graph->get_weight(u);
	total_weight += weight_u;
	tabuList[u] = iteration;

	moveFreeToSolutionPartition(u, code);

	vector<int> &neighbors = graph->get_vertex_adjList(u);

    // update tightness
	if(code == 0) { // placing vertex u in solution A
		for(int neighbor : neighbors) { tightness_A[neighbor]++; }
	}
	else { // placing vertex u in solution B
		for(int neighbor : neighbors) { tightness_B[neighbor]++; }
	}

	checkFreePartition();
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

void BipartiteSolution::swapVertices(int vertex_to_add, int vertex_to_remove, int code, int iteration) { // code == 0 for partition A and code != 0 for partition B	
    // move the vertex_to_remove to the Free Partition
    removeVertex(vertex_to_remove, code);
    // move the vertex_to_add to the solution Partition
    addVertex(vertex_to_add, code, iteration);	
}

bool BipartiteSolution::swap1_1(int code, int iteration) { // code == 0 for partition A and code != 0 for partition B
    // improvement is a variable to get the best possible neighbor_size (in the free partition) 
    // and the worse neighbor_size (in the solution partition)
	int neighbor_size, improvement = 0;  
    int vertex, vertex_to_remove = -1, best_vertex = -1;

	if(code == 0) {
        for(int iter = solution_size_A; iter < solution_size_A + free_size_A; iter++) { // initial loop to get the best tightness among the free vertices
			vertex = solution_A[iter];
            neighbor_size = graph->get_vertex_adjList(vertex).size();
			if(neighbor_size > improvement && validateTabu(vertex, iteration)) { // possible candidate for swap(1,1)
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
			if(neighbor_size > improvement && validateTabu(vertex, iteration)) { // possible candidate for swap(1,1)
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
        swapVertices(best_vertex, vertex_to_remove, code, iteration);
        return true;
    } 

	return false;
}

bool BipartiteSolution::swap1_k(int code, int iteration) { // code == 0 for partition A and code != 0 for partition B
    int vertex, vertexChosen;
	bool operationSuccessFlag = false;
	vector<int> possibleVerticesToSolution;

	if(code == 0 && solution_size_B > 1) {
		// initial loop to get the possible vertices to join the solution
        for(unsigned int iter = solution_size_A + free_size_A; iter < solution_A.size() - removed_size_A; iter++) { 
			vertex = solution_A[iter];
			if (tightness_B[vertex] == solution_size_B - 1 && validateTabu(vertex, iteration)) {
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
				removeVertex(vertex, 1);
			}
		}

		for(unsigned int iter = 0; iter < possibleVerticesToSolution.size(); iter++) { 
			vertex = possibleVerticesToSolution[iter];
			if (tightness_B[vertex] == solution_size_B && position_A[vertex] >= solution_size_A) {
				addVertex(vertex, code, iteration);
			}
		}
	}
	else if(solution_size_A > 1) {
		// initial loop to get the possible vertices to join the solution
        for(unsigned int iter = solution_size_B + free_size_B; iter < solution_B.size() - removed_size_B; iter++) { 
			vertex = solution_B[iter];
			if (tightness_A[vertex] == solution_size_A - 1 && validateTabu(vertex, iteration)) {
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

		for(unsigned int iter = 0; iter < possibleVerticesToSolution.size(); iter++) { 
			vertex = possibleVerticesToSolution[iter];
			if (tightness_A[vertex] == solution_size_A && position_B[vertex] >= solution_size_B) {
				addVertex(vertex, code, iteration);
			}
		}
	}

	possibleVerticesToSolution.clear();

	return operationSuccessFlag;
}

bool BipartiteSolution::addPairOfVertices(int iteration) { // verify if can add a pair of vertices to the solution (one to partition A and one to partition B)
	if(free_size_A > 0 && free_size_B > 0) {
		uniform_int_distribution<int> distribution(0, free_size_A - 1);
		int iter = solution_size_A + distribution(BipartiteGenerator);
		int u, vertex1 = -1, vertex2 = -1, loop = free_size_A;
		int vertex_weight, neighbor_weight;
		
		while(loop--) { // check if there is a pair to join the solution
			if(iter < solution_size_A + free_size_A) u = solution_A[iter];
			else u = solution_A[(iter % (solution_size_A + free_size_A)) + solution_size_A]; 

			vertex_weight = graph->get_weight(u);

			vector<int> &neighbors = graph->get_vertex_adjList(u);
			for(int neighbor : neighbors) {
				// check if there is an edge between these free vertices
				neighbor_weight = graph->get_weight(neighbor);
				if((position_B[neighbor] >= solution_size_B && position_B[neighbor] < solution_size_B + free_size_B) && (0 <= vertex_weight + neighbor_weight)) {
					vertex1 = u;
					vertex2 = neighbor;
					break;
				}
			}
			if(vertex1 != -1 && vertex2 != -1) {
				addVertex(vertex1, 0, iteration);
				addVertex(vertex2, 1, iteration);

				return true;
			}
			iter++;
		}
	}
	return false;
}

/*void BipartiteSolution::tabuIteration(int neighborCode, int iteration) { // run a tabu iteration
	srand (time(NULL));
	int partitionCode;

	switch(neighborCode) {
		case 1:
			addPairOfVertices(iteration);
			break;
		case 2:
			partitionCode = rand() % 2; // try to randomly choose a partition
			swap1_1(partitionCode, iteration);
			break;
		case 3:
			// TODO ANALYZE THE LINE BELOW
			//if (solution_size_A <= 1 || solution_size_B <= 1) { k++; break; }
			partitionCode = rand() % 2; // try to randomly choose a partition
			swap1_k(partitionCode, iteration);			
			break;
		default:
			break;
	}
	
}*/

// create a random solution for tabu-vnd
void BipartiteSolution::randomConstructive() { 
	srand (time(NULL));
	int vertexToAdd, idx, iteration = 0;

	while(free_size_A > 0 && free_size_B > 0) { // can generate a solution with an extra vertex in solution A
		// first randomly choose a free vertex to partition A
		idx = solution_size_A + (rand() % free_size_A);
		//cout << "idx: " << solution_size_A << endl;
		vertexToAdd = solution_A[idx];
		//cout << "Vertex to add: " << vertexToAdd << endl;
		addVertex(vertexToAdd, 0, iteration);
		
		// then select randomly a free vertex to partition B if there is a free vertex
		if(free_size_B > 0) {
			idx = solution_size_B + (rand() % free_size_B);
			//cout << "idx: " << idx << endl;
			vertexToAdd = solution_B[idx];
			//cout << "Vertex to add: " << vertexToAdd << endl;
			addVertex(vertexToAdd, 1, iteration);
		}
	}

	assert(checkIntegrity());
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

			if(iter < solution_size_A) {
				checkSolutionPartition();
			}
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

			if(iter < solution_size_B) {
				checkSolutionPartition();
			}
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

// perturb the solution removing random vertices
void BipartiteSolution::perturb() { 
	srand (time(NULL));
	int verticesToRemove = (rand() % 2) + 1;
	int v, idx;

	while(verticesToRemove && verticesToRemove <= solution_size_A && verticesToRemove <= solution_size_B) {
		// removing random vertex from partition A
		idx = rand() % solution_size_A;
		v = solution_A[idx];
		removeVertex(v, 0);

		// removing random vertex from partition B
		idx = rand() % solution_size_B;
		v = solution_B[idx];
		removeVertex(v, 1);

		verticesToRemove--;
	}
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
}