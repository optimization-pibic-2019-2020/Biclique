#include <iostream>
#include <vector>
#include <algorithm> 
#include <random>
#include <ctime>
#include "Solution.hpp"
//#define NDEBUG
#include <assert.h>
using namespace std;

// check if the solution is maximal
bool Solution::isMaximal(int code) { // code == 0 for solution A and code != 0 for solution B
	if(code == 0) return free_size_A == 0;
	else return free_size_B == 0; 
}

// check if both solutions have the same size
bool Solution::checkBicliqueSize() {
	if(solution_size_A == solution_size_B) return true;

	return false;
}

// return the solution (weight)
int Solution::getTotalWeight() {
	return total_weight;
}

// return how many edges were removed
int Solution::getRemovedEdges() {
	return removed_edges;
}

// check if both vertices are neighbors
bool Solution::isNeighbor(int vertex1, int vertex2) { // checks if two vertices are neighbors
	vector<int> &neighbors = graph->get_vertex_adjList(vertex1);
	for(int neighbor : neighbors) {
		if(neighbor == vertex2) return true;
	}

	return false;
}

// check if both vertices have the same neighbor
bool Solution::sameNeighbor(int vertex1, int vertex2, int code) { // checks if two vertices share the same neighbor in a partition's solution
	int solutionVertex;
	if(code == 0) { // verify partition A
		for(int iter = 0; iter < solution_size_A; iter++) {
			solutionVertex = solution_A[iter];
			if(isNeighbor(vertex1, solutionVertex) != isNeighbor(vertex2, solutionVertex)) return false;
		}
	}
	else { // verify partition B
		for(int iter = 0; iter < solution_size_B; iter++) {
			solutionVertex = solution_B[iter];
			if(isNeighbor(vertex1, solutionVertex) != isNeighbor(vertex2, solutionVertex)) return false;
		}
	}

	return true;
}

// move a vertex from free partition to the solution 
void Solution::moveFreeToSolutionPartition(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	assert(u < graph->getV());

	if(code == 0) {
		// current position of u in the solution_A vector
		int pos_u = position_A[u];

		// new position of u in the solution_A vector
		int new_pos_u = solution_size_A;

		// first vertex of the second partition
		int j = solution_A[solution_size_A];

		// ensures u is in the free partition of the solution vector A
		assert((solution_size_A <= pos_u) && (pos_u < solution_size_A + free_size_A));
		assert(tightness_B[u] == solution_size_B);

		// swap u with the first vertex of the second partition
		swap(solution_A[pos_u], solution_A[new_pos_u]);
		position_A[u] = new_pos_u;
		position_A[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// first partition
		solution_size_A++;
		free_size_A--;
	}
	else {
		// current position of u in the solution_B vector
		int pos_u = position_B[u];

		// new position of u in the solution_B vector
		int new_pos_u = solution_size_B;

		// first vertex of the second partition
		int j = solution_B[solution_size_B];

		// ensures u is in the free partition of the solution vector B
		assert((solution_size_B <= pos_u) && (pos_u < solution_size_B + free_size_B));
		assert(tightness_A[u] == solution_size_A);

		// swap u with the first vertex of the second partition
		swap(solution_B[pos_u], solution_B[new_pos_u]);
		position_B[u] = new_pos_u;
		position_B[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// first partition
		solution_size_B++;
		free_size_B--;
	}
}

// move a vertex from free partition to the non free partition
void Solution::moveFreeToNonFreePartition(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	assert(u < graph->getV());

	if(code == 0) {
		// current position of u in the solution vector A
		int pos_u = position_A[u];

		// new position of u in the solution vector A
		int new_pos_u = solution_size_A + free_size_A - 1;

		// last vertex of the second partition
		int j = solution_A[solution_size_A + free_size_A - 1];

		// ensures u is in the free partition of the solution vector A
		assert((solution_size_A <= pos_u) && (solution_size_A + free_size_A > pos_u));

		// swap u with the last vertex of the second partition
		swap(solution_A[pos_u], solution_A[new_pos_u]);
		position_A[u] = new_pos_u;
		position_A[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// second partition
		free_size_A--;
	}
	else {
		// current position of u in the solution vector B
		int pos_u = position_B[u];

		// new position of u in the solution vector B
		int new_pos_u = solution_size_B + free_size_B - 1;

		// last vertex of the second partition
		int j = solution_B[solution_size_B + free_size_B - 1];

		// ensures u is in the free partition of the solution vector B
		assert((solution_size_B <= pos_u) && (solution_size_B + free_size_B > pos_u));

		// swap u with the last vertex of the second partition
		swap(solution_B[pos_u], solution_B[new_pos_u]);
		position_B[u] = new_pos_u;
		position_B[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// second partition
		free_size_B--;
	}
}

// move a vertex from the solution to the free partition
void Solution::moveSolutionToFreePartition(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	assert(u < graph->getV());

	if(code == 0) {
		// current position of u in the solution vector A
		int pos_u = position_A[u];

		// new position of u in the solution vector A
		int new_pos_u = solution_size_A - 1;

		// last vertex of the first partition
		int j = solution_A[solution_size_A - 1];

		// ensures u is in the solution partition of the solution vector A
		assert(pos_u < solution_size_A);

		// swap u with the last vertex of the second partition
		swap(solution_A[pos_u], solution_A[new_pos_u]);
		position_A[u] = new_pos_u;
		position_A[j] = pos_u;

		// change the boundary between the blocks to make u the first vertex of the
		// second partition
		solution_size_A--;
		free_size_A++;
	}
	else {
		// current position of u in the solution vector B
		int pos_u = position_B[u];

		// new position of u in the solution vector B
		int new_pos_u = solution_size_B - 1;

		// last vertex of the first partition
		int j = solution_B[solution_size_B - 1];

		// ensures u is in the solution partition of the solution vector B
		assert(pos_u < solution_size_B);

		// swap u with the last vertex of the second partition
		swap(solution_B[pos_u], solution_B[new_pos_u]);
		position_B[u] = new_pos_u;
		position_B[j] = pos_u;

		// change the boundary between the blocks to make u the first vertex of the
		// second partition
		solution_size_B--;
		free_size_B++;
	}
} 

// move a vertex from non free partition to the free partition
void Solution::moveNonFreeToFreePartition(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	assert(u < graph->getV());

	if(code == 0) {
		// current position of u in the solution vector A
		int pos_u = position_A[u];

		// new position of u in the solution vector A
		int new_pos_u = solution_size_A + free_size_A;

		// first vertex of the third partition
		int j = solution_A[solution_size_A + free_size_A];

		// ensures u is in the non free partition of the solution vector A
		assert(pos_u >= solution_size_A + free_size_A);

		// swap u with the last vertex of the second partition
		swap(solution_A[pos_u], solution_A[new_pos_u]);
		position_A[u] = new_pos_u;
		position_A[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// second partition
		free_size_A++;
	}
	else {
		// current position of u in the solution vector B
		int pos_u = position_B[u];

		// new position of u in the solution vector B
		int new_pos_u = solution_size_B + free_size_B;

		// first vertex of the third partition
		int j = solution_B[solution_size_B + free_size_B];

		// ensures u is in the non free partition of the solution vector B
		assert(pos_u >= solution_size_B + free_size_B);

		// swap u with the last vertex of the second partition
		swap(solution_B[pos_u], solution_B[new_pos_u]);
		position_B[u] = new_pos_u;
		position_B[j] = pos_u;

		// change the boundary between the blocks to make u the last vertex of the
		// second partition
		free_size_B++;
	}
}

// move a vertex to the removed partition
void Solution::moveVertexToRemovedVertices(int u, int code) { // code == 0 for solution A and code != 0 for solution B
	int V = graph->getV();
	assert(u < V);

	if(code == 0) {
		// current position of u in the solution vector A
		int pos_u = position_A[u];

		// new position of u in the solution vector A
		int new_pos_u = V - removed_size_A - 1;

		// last vertex of the third partition
		int j = solution_A[new_pos_u];

		// ensures u is not a removed vertex of the solution vector A
		assert(pos_u <= new_pos_u);

		// swap u with the last vertex of the third partition
		swap(solution_A[pos_u], solution_A[new_pos_u]);
		position_A[u] = new_pos_u;
		position_A[j] = pos_u;

		// change the boundary between the blocks to make u the first vertex of the forth partition
		removed_size_A++;
	}
	else {
		// current position of u in the solution vector B
		int pos_u = position_B[u];

		// new position of u in the solution vector B
		int new_pos_u = V - removed_size_B - 1;

		// last vertex of the third partition
		int j = solution_B[new_pos_u];

		// ensures u is not a removed vertex of the solution vector B
		assert(pos_u <= new_pos_u);

		// swap u with the last vertex of the third partition
		swap(solution_B[pos_u], solution_B[new_pos_u]);
		position_B[u] = new_pos_u;
		position_B[j] = pos_u;

		// change the boundary between the blocks to make u the first vertex of the forth partition
		removed_size_B++;
	}
}

// checks the integrity of the solution in general
bool Solution::checkIntegrity() {  
	vector<vector<int>> &neighbors = graph->get_adjList();
	int V = graph->getV();

	if(abs(solution_size_A - solution_size_B) > 1) { // both solution sizes cannot be different  
		cout << "abs(solution_size_A - solution_size_B) > 1" << endl;
		return false; 
	}

	for(int idx = 0; idx < solution_size_A; idx++) {
		int vertex_u = solution_A[idx];
		int neighbor_amount_B = 0; // variable to check the amount of neighbors in solution B

		if(tightness_A[vertex_u] > 0 || tightness_B[vertex_u] != solution_size_B) {
			cout << "tA: " << tightness_A[vertex_u] << endl;
			cout << "tb: " << tightness_B[vertex_u] << " != " << solution_size_B << endl;
			cout << "A: " << vertex_u << endl;
			cout << "B: " << endl;
			for(int idx = 0; idx < solution_size_B; idx++) {
				int vertex_b = solution_B[idx];
				cout << vertex_b << endl;
			}
			cout << "Tightness error" << endl;
			return false;
		}


		// verify if every vertex in solution A is disconnected to all the vertices in solution A
		// and if every vertex in solution B is connected to every vertex in solution A 
		for(int neighbor : neighbors[vertex_u]) {
			if(find(solution_A.begin(), solution_A.begin() + solution_size_A, neighbor) != solution_A.begin() + solution_size_A) {
				cout << "There are neighbors in Solution A" << endl;
				return false;
			}
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) != solution_B.begin() + solution_size_B) { 
				neighbor_amount_B++;
			}
		}

		if(neighbor_amount_B != solution_size_B || neighbor_amount_B != tightness_B[vertex_u]) {
			cout << neighbor_amount_B << " != " <<  solution_size_B << endl;
			cout << neighbor_amount_B << " != " <<  tightness_B[vertex_u] << endl;
			cout << "A: " << vertex_u << endl;
			cout << "B: " << endl;
			for(int idx = 0; idx < solution_size_B; idx++) {
				int vertex_b = solution_B[idx];
				cout << vertex_b << endl;
			}
			cout << "The quantity of neighbors in solution B of a vertex in Solution A != solution_size_B\nOR\nThe quantity of neighbors in solution B of a vertex in Solution A != its tightness_B" << endl;
			return false; 
		}
	}

	for(int idx = 0; idx < solution_size_B; idx++) {
		int vertex_u = solution_B[idx];
		int neighbor_amount_A = 0; // variable to check the amount of neighbors in solution A

		if(tightness_B[vertex_u] > 0 || tightness_A[vertex_u] != solution_size_A) {
			cout << "Tightness error" << endl;
			return false;
		}

		// verify if every vertex in solution B is disconnected to all the vertices in solution B
		// and if every vertex in solution A is connected to every vertex in solution B 
		for(int neighbor : neighbors[vertex_u]) {
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) != solution_B.begin() + solution_size_B) {
				cout << "There are neighbors in Solution B" << endl;
				return false;
			}
			if(find(solution_A.begin(), solution_A.begin() + solution_size_A, neighbor) != solution_A.begin() + solution_size_A) { 
				neighbor_amount_A++;
			}
		}

		if(neighbor_amount_A != solution_size_A || neighbor_amount_A != tightness_A[vertex_u]) {
			cout << "The quantity of neighbors in solution A of a vertex in Solution B != solution_size_A\nOR\nThe quantity of neighbors in solution A of a vertex in Solution B != its tightness_A" << endl;
			return false; 
		}	
	}

	for(int idx = solution_size_A; idx < solution_size_A + free_size_A; idx++) { // verify if every free vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] > 0 || tightness_B[vertex] != solution_size_B) {
			cout << "Free Partition A does not meet the requirements" << endl;
			return false;
		}
	}

	for(int idx = solution_size_B; idx < solution_size_B + free_size_B; idx++) { // verify if every free vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] > 0 || tightness_A[vertex] != solution_size_A) {
			cout << "Free Partition B does not meet the requirements" << endl;
			return false;
		}
	}

	for(int idx = solution_size_A + free_size_A; idx < V - removed_size_A; idx++) { // verify if every non solution vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] == 0 && tightness_B[vertex] == solution_size_B) {
			cout << "Non Free Partition A does not meet the requirements" << endl;
			return false;
		}
	}

	for(int idx = solution_size_B + free_size_B; idx < V - removed_size_B; idx++) { // verify if every non solution vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] == 0 && tightness_A[vertex] == solution_size_A) {
			cout << "Non Free Partition B does not meet the requirements" << endl;
			return false;
		}
	}

	return true;
}

// function that creates the probability of each vertex based in linear bias function (1.0 / bias_rank)
void Solution::createRclProbability() { 
	double bias_rank; // variable that represents the rank of each element in rcl
	vector<int> &weight = graph->get_weight_list();
	vector<vector<int>> &neighbors = graph->get_adjList();
	
	for(unsigned int iter = 0; iter < rclList.size(); iter++) {
		bias_rank = 1.0 / ((double) (weight[rclList[iter]] + neighbors[iter].size()));
		rclListProbability.push_back(1.0 / bias_rank);
	}
}

// remove all the edges from the vertex to be removed
void Solution::removeVertexFromGraph(int vertex) { 
	int neighbor;
	vector<int> &vertex_neighbors = graph->get_vertex_adjList(vertex);
	
	removed_edges += vertex_neighbors.size(); // counting the quantity of edges that the removed vertex has

	for(unsigned int i = 0; i < vertex_neighbors.size(); i++) { // first loop to find all the neighbors
		neighbor = vertex_neighbors[i];
		vector<int> &neighbor_adjList = graph->get_vertex_adjList(neighbor);

		for(unsigned int j = 0; j < neighbor_adjList.size(); j++) { // second loop to erase the edge from the neighbor that connects to the removed vertex 
			if(neighbor_adjList[j] == vertex) { 
				graph->removeVertexFromAdjList(neighbor, j); // removes the edge between the removed vertex and its neighbor

				// recalculates neighbor informations due to the vertex remotion
				graph->calculateHIndex(neighbor);  
				graph->calculateAccumulatedSum(neighbor);

				break;
			}
		}
	}

	graph->clearVertexAdjList(vertex); // remove all the edges from the vertex removed
}

// try to predict the best biclique possible with the vertex with h_index and accumulatedSum
int Solution::predictBicliqueWeight(int vertex) { 
	int best_neighbor_weight = 0, neighbor, index, h_index = graph->getVertexHIndex(vertex), vertex_weight = graph->get_weight(vertex);
	vector<int> neighbors = graph->get_vertex_adjList(vertex);
	vector<int> vertex_accumulatedSum = graph->get_vertex_accumulatedSum(vertex);

	if(neighbors.empty()) return 0;
	
	for(unsigned int idx = 0; idx < neighbors.size(); idx++) { // tries to predict the best biclique using h_index and accumulatedSum the neighbors
		neighbor = neighbors[idx];
		vector<int> neighbor_accumulatedSum = graph->get_vertex_accumulatedSum(neighbor);

		if((unsigned int) h_index <= neighbor_accumulatedSum.size()) index = h_index - 1;
		else index = neighbor_accumulatedSum.size() - 1;
		
		if(neighbor_accumulatedSum[index] >= best_neighbor_weight) best_neighbor_weight = neighbor_accumulatedSum[index];
	}
	
	return vertex_accumulatedSum[h_index - 1] + best_neighbor_weight - vertex_weight; 
}

void Solution::restartAm(double beta) {
	int u, t;
	for(int idx = 0; idx < solution_size_A; idx++) { 
		u = solution_A[idx];
		t = solution_B[idx];
		
		graph->setVertexAm(u, beta);
		graph->setVertexAm(t, beta);
	}

	for(long unsigned int idx = solution_size_A; idx < solution_A.size() - removed_size_A; idx++) { 
		u = solution_A[idx];
		graph->setVertexAm(u, 1);
	}

	for(long unsigned int idx = solution_size_B; idx < solution_B.size() - removed_size_B; idx++) { 
		u = solution_B[idx];
		graph->setVertexAm(u, 1);
	}
}

void Solution::updateAm() {
	int u, t, uAmValue, tAmValue;
	for(int idx = 0; idx < solution_size_A; idx++) { 
		u = solution_A[idx];
		uAmValue = graph->getVertexAm(u);

		t = solution_B[idx];
		tAmValue = graph->getVertexAm(t);

		graph->setVertexAm(u, uAmValue + 1);
		graph->setVertexAm(t, tAmValue + 1);
	}
}

// show the current solution
void Solution::printSolution() {
	cout << "{ {";
	for(int idx = 0; idx < solution_size_A; idx++) {
		if(idx == solution_size_A - 1) {
			cout << solution_A[idx] + 1 << "}";
			break;
		}
		cout << solution_A[idx] + 1 << ", ";
	}
	cout << ", {";
	for(int idx = 0; idx < solution_size_B; idx++) {
		if(idx == solution_size_B - 1) {
			cout << solution_B[idx] + 1 << "}";
			break;
		}
		cout << solution_B[idx] + 1 << ", ";
	}
	cout << " }" << endl;
}