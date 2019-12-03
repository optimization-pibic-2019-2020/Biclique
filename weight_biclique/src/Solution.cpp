#include <bits/stdc++.h>
#include <random>
#include "Solution.hpp"

using namespace std;

Solution::Solution(Graph *graph) { // initialize all the variables and vectors
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
		mu_A[idx] = graph->get_weight(idx);
		mu_B[idx] = graph->get_weight(idx);
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

// add a vertex to solution A or B
void Solution::addVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B 
	int weight_u = graph->get_weight(u);
	total_weight += weight_u;

	moveFreeToSolutionPartition(u, code);
	if(code == 0) moveFreeToNonFreePartition(u, 1);
	else moveFreeToNonFreePartition(u, 0);

	vector<int> adjListTemp = graph->get_vertex_adjList(u);

	if(code == 0) { // placing vertex u in solution A
		for(int neighbor : adjListTemp) {
			// increase the tighness of each neighbor by one
			tightness_A[neighbor]++;

			mu_A[neighbor] -= weight_u;

			// if the neighbor is in the free partition, move to non free partition
			int neighbor_pos = position_A[neighbor];
			if ((solution_size_A <= neighbor_pos) && (solution_size_A + free_size_A > neighbor_pos)) {
				moveFreeToNonFreePartition(neighbor, 0);
			}
		}
	}
	else { // placing vertex u in solution B
		for(int neighbor : adjListTemp) {
			// increase the tighness of each neighbor by one
			tightness_B[neighbor]++;

			mu_B[neighbor] -= weight_u;

			// if the neighbor is in the free partition, move to non free partition
			int neighbor_pos = position_B[neighbor];
			if ((solution_size_B <= neighbor_pos) && (solution_size_B + free_size_B > neighbor_pos)) {
				moveFreeToNonFreePartition(neighbor, 1);
			}
		}
	}
}

// add a random vertex to solution A or solution B
void Solution::addRandomVertex(int code) { // code == 0 for solution A and code != 0 for solution B
	default_random_engine generator; 
	uniform_int_distribution<int> distribution;
	int vertex;
	assert(!isMaximal(code));

	// generate a random number between [0, free_size_ - 1]
	if(code == 0) uniform_int_distribution<int> distribution(0, free_size_A - 1);
	else uniform_int_distribution<int> distribution(0, free_size_B - 1);
	int free_pos = distribution(generator);

	if(code == 0) vertex = solution_A[solution_size_A + free_pos];
	else vertex = solution_B[solution_size_B + free_pos];

	addVertex(vertex, code);
} 

bool Solution::checkIntegrity() { // checks the integrity of the solution 
	vector<vector<int>> adjListTemp = graph->get_adjList();

	if(solution_size_A != solution_size_B) return false; // both solution sizes cannot be different  

	for(int idx = 0; idx < solution_size_A; idx++) {
		int vertex_u = solution_A[idx];
		int vertex_t = solution_B[idx];

		if (tightness_A[vertex_u] > 0 || tightness_B[vertex_u] < solution_size_B) return false;
		if (tightness_B[vertex_t] > 0 || tightness_A[vertex_t] < solution_size_A) return false;

		// verify if every vertex in solution A is disconnected to all the vertices in solution A
		// and if every vertex in solution B is connected to every vertex in solution A 
		for(int neighbor : adjListTemp[vertex_u]) {
			if(find(solution_A.begin(), solution_A.begin() + solution_size_A, neighbor) != solution_A.begin() + solution_size_A) {
				return false;
			}
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) == solution_B.begin() + solution_size_B) { 
				return false;
			}
		}

		// verify if every vertex in solution B is disconnected to all the vertices in solution B
		for(int neighbor : adjListTemp[vertex_t]) {
			if(find(solution_B.begin(), solution_B.begin() + solution_size_B, neighbor) != solution_B.begin() + solution_size_B) {
				return false;
			}
		}
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

	for(int idx = solution_size_A + free_size_A; idx < graph->getV(); idx++) { // verify if every non solution vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] == 0 && tightness_B[vertex] == solution_size_B) {
			return false;
		}
	}

	for(int idx = solution_size_B + free_size_B; idx < graph->getV(); idx++) { // verify if every non solution vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] == 0 && tightness_A[vertex] == solution_size_A) {
			return false;
		}
	}

	return true;
}

void Solution::generateRandomSolution()  {

}


/*void checkFreePartition() { 
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

void checkNonFreePartition() { 
	int vertex, u;

	int n = graph->getV();
	for(vertex = solution_size_A + free_size_A; vertex < n; vertex++) { // verify the free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] == 0 && tightness_B[u] == solution_size_B) moveNonFreeToFreePartition(u, 0);  
	}
	for(vertex = solution_size_B + free_size_B; vertex < n; vertex++) { // verify the free partition of solution A
		u = solution_B[vertex];
		if(tightness_B[u] == 0 && tightness_A[u] == solution_size_A) moveNonFreeToFreePartition(u, 1);  
	}
}*/
