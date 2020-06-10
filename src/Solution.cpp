#include <bits/stdc++.h>
#include <random>
#include "Solution.hpp"
#define NDEBUG
#include <assert.h>
using namespace std;

random_device device;
mt19937 generator(device());

Solution::Solution(Graph *graph) { // initialize all the variables and vectors
	this->graph = graph;
	int idx_weight, V = graph->getV();

	solution_A.resize(V);
	solution_B.resize(V);

	solution_size_A = 0;
	solution_size_B = 0;

	free_size_A = V;
	free_size_B = V;

	position_A.resize(V);
	position_B.resize(V);

	tightness_A.resize(V);
	tightness_B.resize(V);

	mu_A.resize(V);
	mu_B.resize(V);

	total_weight = 0;	

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

// check if the solution is maximal
bool Solution::isMaximal(int code) { // code == 0 for solution A and code != 0 for solution B
	if(code == 0) return free_size_A == 0;
	else return free_size_B == 0; 
}

bool Solution::checkBicliqueSize() {
	if(solution_size_A == solution_size_B) return true;

	return false;
}

int Solution::getTotalWeight() {
	return total_weight;
}

bool Solution::isNeighbor(int vertex1, int vertex2) { // checks if two vertices are neighbors
	vector<int> &neighbors = graph->get_vertex_adjList(vertex1);
	for(int neighbor : neighbors) {
		if(neighbor == vertex2) return true;
	}

	return false;
}

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
void Solution::checkNonFreePartition() { 
	int vertex, u, V = graph->getV();
	
	for(vertex = solution_size_A + free_size_A; vertex < V; vertex++) { // verify the non-free partition of solution A
		u = solution_A[vertex];
		if(tightness_A[u] == 0 && tightness_B[u] == solution_size_B) {
			moveNonFreeToFreePartition(u, 0);
			if(vertex >= solution_size_A + free_size_A) vertex--;  
		}
	}
	for(vertex = solution_size_B + free_size_B; vertex < V; vertex++) { // verify the non-free partition of solution B
		u = solution_B[vertex];
		if(tightness_B[u] == 0 && tightness_A[u] == solution_size_A) {
			moveNonFreeToFreePartition(u, 1);  
			if(vertex >= solution_size_B + free_size_B) vertex--; 
		}
	}
}

// remove a vertex from solution A or B
void Solution::removeVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B
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
void Solution::addVertex(int u, int code) { // code == 0 for solution A and code != 0 for solution B 
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

// add a random vertex to solution A or solution B
void Solution::addRandomVertex(int code) { // code == 0 for solution A and code != 0 for solution B
	if(isMaximal(code)) return;
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

	for(int idx = solution_size_A + free_size_A; idx < V; idx++) { // verify if every non solution vertex in solution_A is correct
		int vertex = solution_A[idx];
		if(tightness_A[vertex] == 0 && tightness_B[vertex] == solution_size_B) {
			cout << "Non Free Partition A does not meet the requirements" << endl;
			return false;
		}
	}

	for(int idx = solution_size_B + free_size_B; idx < V; idx++) { // verify if every non solution vertex in solution_B is correct
		int vertex = solution_B[idx];
		if(tightness_B[vertex] == 0 && tightness_A[vertex] == solution_size_A) {
			cout << "Non Free Partition B does not meet the requirements" << endl;
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

void Solution::restartSolution() {
	this->graph = graph;
	int idx_weight, V = graph->getV();

	total_weight = 0;	

	solution_size_A = 0;
	solution_size_B = 0;

	free_size_A = V;
	free_size_B = V;


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

void Solution::oneImprovement(int vertex, int code) { // code == 0 for partition A and code != 0 for partition B
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

bool Solution::swap1_1(int code) { // code == 0 for partition A and code != 0 for partition B
	int V = graph->getV(), vertex, best_vertex, improvement = 0; // improvement means how much weight the partition will get after the swap
	if(code == 0) {
		for(int iter = solution_size_A + free_size_A; iter < V; iter++) {
			vertex = solution_A[iter];
			if(tightness_A[vertex] == 1 && tightness_B[vertex] == solution_size_B && mu_A[vertex] > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = mu_A[vertex];
			}
		}
	}
	else {
		for(int iter = solution_size_B + free_size_B; iter < V; iter++) {
			vertex = solution_B[iter];
			if(tightness_B[vertex] == 1 && tightness_A[vertex] == solution_size_A && mu_B[vertex] > improvement) { // possible candidate for swap(1,1)
				best_vertex = vertex;
				improvement = mu_B[vertex];
			}
		}
	}

	if(improvement > 0) {
		oneImprovement(best_vertex, code);
		return true;
	}

	return false;
}

bool Solution::swap2_2(int code) { // code == 0 for partition A and code != 0 for partition B
	int V = graph->getV(), vertex1, vertex2, improvement = 0; // improvement means how much weight the partition will get after the swap
	if(code == 0 && solution_size_A >= 2) {
		for(int iter = solution_size_A + free_size_A; iter < V; iter++) {
			vertex1 = solution_A[iter];
			if(tightness_A[vertex1] == 1 && tightness_B[vertex1] == solution_size_B) { // possible candidate for swap2_2
				for(int iter2 = iter; iter2 < V; iter2++) {
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
		for(int iter = solution_size_B + free_size_B; iter < V; iter++) {
			vertex1 = solution_B[iter];
			if(tightness_B[vertex1] == 1 && tightness_A[vertex1] == solution_size_A) { // possible candidate for swap2_2
				for(int iter2 = iter; iter2 < V; iter2++) {
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

	return false;
}

bool Solution::addFirstVertex() { // verify if can add a pair of vertices to the solution (one to partition A and one to partition B)
	if(free_size_A > 0 && free_size_B > 0) {
		uniform_int_distribution<int> distribution(0, free_size_A - 1);
		int iter = solution_size_A + distribution(generator);
		int u, vertex1 = -1, vertex2 = -1, loop = free_size_A;
		
		while(loop--) { // check if there is a pair to join the solution
			if(iter < free_size_A) u = solution_A[iter];
			else u = solution_A[(iter % free_size_A) + solution_size_A]; 

			vector<int> &neighbors = graph->get_vertex_adjList(u);
			for(int neighbor : neighbors) {
				// check if there is an edge between these free vertices
				if(position_B[neighbor] >= solution_size_B && position_B[neighbor] < solution_size_B + free_size_B && 0 <= mu_A[u] + mu_B[neighbor]) {
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

void Solution::VND(int K) { // run VND iterations
	srand (time(NULL));
	int k = 1, code;
	while(k <= K) {
		switch(k) {
			case 1:
				if(addFirstVertex()) k = 1;
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

void Solution::createRclProbability() { // function that creates the probability of each vertex based in linear bias function (1.0 / bias_rank)
	double bias_rank;
	vector<int> &weight = graph->get_weight_list();
	for(int iter = 0; iter < rclList.size(); iter++) {
		bias_rank = 1.0 / ((double) weight[rclList[iter]]);
		rclListProbability.push_back(1.0 / bias_rank);
	}
}

void Solution::rclConstruction(int code, double alpha) { // construct the restricted candidate list for a specific partition and choose a random vertex from it to put in the solution
	if(code == 0 && free_size_A != 0) { // partition A
		int iter, vertex, vertex_weight;
		int c_min = 100000000, c_max = -1; // variables to get the interval for the quality-based RCL list
		vector<int> &weight = graph->get_weight_list();
		rclList.clear();
		rclListProbability.clear();

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
		int pos = distribution(generator);
		vertex = rclList[pos]; 
		addVertex(vertex, 0);
		checkFreePartition();

		/*cout << "\nRCL probabilities (Bias Function)" << endl;
		for(int i = 0; i < rclListProbability.size(); i++) {
			cout << "Vertex: " << rclList[i] << " Weight: " << (rclList[i] % 200) + 1 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
		}*/
	}
	else if(code != 0 && free_size_B != 0) { // partition B
		int iter, vertex, vertex_weight;
		int c_min = 100000000, c_max = -1; // variables to get the interval for the quality-based RCL list
		vector<int> &weight = graph->get_weight_list();
		rclList.clear();
		rclListProbability.clear();

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
		int pos = distribution(generator);
		vertex = rclList[pos]; 
		addVertex(vertex, 1);
		checkFreePartition();

		/*cout << "\nRCL probabilities (Bias Function)" << endl;
		for(int i = 0; i < rclListProbability.size(); i++) {
			cout << "Vertex: " << rclList[i] << " Weight: " << (rclList[i] % 200) + 1 << " Probability: " << distribution.probabilities()[i] * 100 << "%" << endl;
		}*/
	}
}

void Solution::greedyRandomizedConstructive(double p) { 
	while(free_size_A != 0 && free_size_B != 0) { // can generate a solution with an extra vertex in solution A
		rclConstruction(0, p); 
		rclConstruction(1, p); 
		assert(checkIntegrity());
		assert(checkMu());
	}
}

void Solution::balanceBiclique() { // remove the vertex with the worst weight in solution A to balance the Biclique
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