#include <bits/stdc++.h>
#include <random>
#include "Solution.hpp"

using namespace std;

Solution::Solution(Graph *graph) { // initialize all the variables and vectors
	int idx_weight;
	this->graph = graph;

	solution_A.resize(graph->getV());
	solution_B.resize(graph->getV());

	solution_size_A = 0;
	solution_size_B = 0;

	free_size_A = graph->getV();
	free_size_B = graph->getV();

	position_A.resize(graph->getV());
	position_B.resize(graph->getV());

	tightness_A.resize(graph->getV());
	tightness_B.resize(graph->getV());

	mu_A.resize(graph->getV());
	mu_B.resize(graph->getV());

	total_weight = 0;	

	for(int idx = 0; idx < this->graph->getV(); idx++) {
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

// check if the solution is maximal
bool Solution::isMaximal(int code) { // code == 0 for solution A and code != 0 for solution B
	if(code == 0) return free_size_A == 0;
	else return free_size_B == 0; 
}

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
		assert((solution_size_A <= pos_u) && (solution_size_A + free_size_A > pos_u));

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
		assert((solution_size_B <= pos_u) && (solution_size_B + free_size_B > pos_u));

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

// update the free partiton of solution A and solution B
void Solution::checkFreePartition() { 
	int vertex, u;

	for(vertex = solution_size_A; vertex < solution_size_A + free_size_A; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] > 0 || tightness_B[u] != solution_size_B) moveFreeToNonFreePartition(u, 0);  
	}
	for(vertex = solution_size_B; vertex < solution_size_B + free_size_B; vertex++) { // verify the free partition of solution B
		u = solution_B[vertex];
		if(tightness_A[u] != solution_size_A || tightness_B[u] > 0) moveFreeToNonFreePartition(u, 1);  
	}
}

// update the non-free partiton of solution A and solution B
void Solution::checkNonFreePartition() { 
	int vertex, u, V = graph->getV();
	
	for(vertex = solution_size_A + free_size_A; vertex < V; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] == 0 && tightness_B[u] == solution_size_B) moveNonFreeToFreePartition(u, 0);  
	}
	for(vertex = solution_size_B + free_size_B; vertex < V; vertex++) { // verify the free partition of solution A
		u = solution_B[vertex];
		if(tightness_B[u] == 0 && tightness_A[u] == solution_size_A) moveNonFreeToFreePartition(u, 1);  
	}
}

// add a vertex to solution A or B
void Solution::addVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B 
	int weight_u = graph->get_weight(u);
	total_weight += weight_u;

	moveFreeToSolutionPartition(u, code);

	vector<int> adjListTemp = graph->get_vertex_adjList(u);

	if(code == 0) { // placing vertex u in solution A
		for(int neighbor : adjListTemp) {
			// increase the tighness of each neighbor by one
			tightness_A[neighbor]++;

			// update the value of mu from each neighbor in solution A
			mu_A[neighbor] -= weight_u;
		}
	}
	else { // placing vertex u in solution B
		for(int neighbor : adjListTemp) {
			// increase the tighness of each neighbor by one
			tightness_B[neighbor]++;

			// update the value of mu from each neighbor in solution B
			mu_B[neighbor] -= weight_u;
		}
	}

	checkFreePartition();
}

// add a random vertex to solution A or solution B
void Solution::addRandomVertex(int code) { // code == 0 for solution A and code != 0 for solution B
	random_device device;
	mt19937 generator(device());

	assert(!isMaximal(code));

	int vertex;

	// generate a random number between [0, free_size_ - 1]
	if(code == 0) {
		uniform_int_distribution<int> distribution(0, free_size_A - 1);
		int free_pos = distribution(generator);
		vertex = solution_A[solution_size_A + free_pos];
	}
	else {
		uniform_int_distribution<int> distribution(0, free_size_B - 1);
		int free_pos = distribution(generator);
		vertex = solution_B[solution_size_B + free_pos];
	}

	addVertex(vertex, code);
	checkFreePartition();
} 

bool Solution::checkIntegrity() { // checks the integrity of the solution 
	vector<vector<int>> adjListTemp = graph->get_adjList();
	int V = graph->getV();

	if(solution_size_A != solution_size_B) return false; // both solution sizes cannot be different  

	for(int idx = 0; idx < solution_size_A; idx++) {
		int vertex_u = solution_A[idx];
		int neighbor_amount_B = 0; // variable to check the amount of neighbors in solution B

		if(tightness_A[vertex_u] > 0 || tightness_B[vertex_u] != solution_size_B) return false;

		// verify if every vertex in solution A is disconnected to all the vertices in solution A
		// and if every vertex in solution B is connected to every vertex in solution A 
		for(int neighbor : adjListTemp[vertex_u]) {
			if(find(solution_A.begin(), solution_A.begin() + solution_size_A, neighbor) != solution_A.begin() + solution_size_A) {
				return false;
			}
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) != solution_B.begin() + solution_size_B) { 
				neighbor_amount_B++;
			}
		}

		if(neighbor_amount_B != solution_size_B || neighbor_amount_B != tightness_B[vertex_u]) return false; 
	}

	for(int idx = 0; idx < solution_size_B; idx++) {
		int vertex_u = solution_B[idx];
		int neighbor_amount_A = 0; // variable to check the amount of neighbors in solution A

		if(tightness_B[vertex_u] > 0 || tightness_A[vertex_u] != solution_size_A) return false;

		// verify if every vertex in solution B is disconnected to all the vertices in solution B
		// and if every vertex in solution A is connected to every vertex in solution B 
		for(int neighbor : adjListTemp[vertex_u]) {
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) != solution_B.begin() + solution_size_B) {
				return false;
			}
			if(find(solution_A.begin(), solution_A.begin() + solution_size_A, neighbor) != solution_A.begin() + solution_size_A) { 
				neighbor_amount_A++;
			}
		}

		if(neighbor_amount_A != solution_size_A || neighbor_amount_A != tightness_A[vertex_u]) return false; 
	}

	for(int idx = solution_size_A; idx < solution_size_A + free_size_A; idx++) { // verify if every free vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] > 0 || tightness_B[vertex] != solution_size_B) {
			return false;
		}
	}

	for(int idx = solution_size_B; idx < solution_size_B + free_size_B; idx++) { // verify if every free vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] > 0 || tightness_A[vertex] != solution_size_A) {
			return false;
		}
	}

	for(int idx = solution_size_A + free_size_A; idx < V; idx++) { // verify if every non solution vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] == 0 && tightness_B[vertex] == solution_size_B) {
			return false;
		}
	}

	for(int idx = solution_size_B + free_size_B; idx < V; idx++) { // verify if every non solution vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] == 0 && tightness_A[vertex] == solution_size_A) {
			return false;
		}
	}

	return true;
}

bool Solution::checkMu() { // check if mu_A and mu_B are correct
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
		vector<int> adjListTemp = graph->get_vertex_adjList(u);

		for(int neighbor : adjListTemp) {
			// update the value of mu_A_temp from each neighbor in solution A
			mu_A_temp[neighbor] -= weight_u;
		}
	}

	for(int idx = 0; idx < solution_size_B; idx++) { // calculate mu_B_temp
		int u = solution_B[idx];
		int weight_u = graph->get_weight(u);
		vector<int> adjListTemp = graph->get_vertex_adjList(u);

		for(int neighbor : adjListTemp) {
			// update the value of mu_B_temp from each neighbor in solution B
			mu_B_temp[neighbor] -= weight_u;
		}
	}

	for(int idx = 0; idx < V; idx++) { // check the integrity of each mu_
		if(mu_A[idx] != mu_A_temp[idx] || mu_B[idx] != mu_B_temp[idx]) return false;
	}

	return true;
}

void Solution::generateRandomSolution()  {
	while(free_size_A != 0 && free_size_B != 0) {
		addRandomVertex(0); // add a random vertex in solution A
		addRandomVertex(1); // add a random vertex in solution B
		assert(checkIntegrity());
		assert(checkMu());
	}
}

void Solution::printSolution()  {
	cout << "Partition A:" << endl;
	for(int idx = 0; idx < solution_size_A; idx++) {
		cout << solution_A[idx] << " ";
	}
	cout << "\nPartition B:" << endl;
	for(int idx = 0; idx < solution_size_A; idx++) {
		cout << solution_B[idx] << " ";
	}
	cout << endl;
}
