#include <bits/stdc++.h>
#include <random>
#include "Graph.hpp"

using namespace std;

class Solution {
private:
	Graph *graph; // graph instance

	/* 
	Each solution_ vector is partitioned into three blocks: first vertices in the solution, then 
	the free vertices (vertices that are not adjacent to any vertex in the actual solution_ vector) and 
	the non-solution vertices that are not free 
	*/

	// There will be two solution_ vector

	vector<int> solution_A;
	vector<int> solution_B;

	// size of each solution verticies partition

	int solution_size_A;
	int solution_size_B;

	// size of the free vertices partition

	int free_size_A;
	int free_size_B;

	// for each vertex, the number of adjacent vertices that are in each solution

	vector<int> tightness_A;
	vector<int> tightness_B;

	// position of each vertex in each solution_ vector

	vector<int> position_A;
	vector<int> position_B;

	// weight of each vertex i minus the sum of the weights of its neighbors that are in each solution

	vector<int> mu_A;
	vector<int> mu_B;
	
	// current biclique weight
	
	int total_weight;

	// vector of the restricted candidate list for the greedy randomized construction
	
	vector<int> rclList;

public:
	Solution(Graph *graph);
	int getTotalWeight();
	bool isNeighbor(int vertex1, int vertex2);
	bool sameNeighbor(int vertex1, int vertex2, int code);
	bool checkBicliqueSize();
	void checkFreePartition();
	void checkNonFreePartition();
	void moveFreeToSolutionPartition(int u, int code);
	void moveFreeToNonFreePartition(int u, int code);
 	void moveSolutionToFreePartition(int u, int code);
	void moveNonFreeToFreePartition(int u, int code); 
	void addVertex(int u, int code);
	void removeVertex(int u, int code);
	bool isMaximal(int code);
	void addRandomVertex(int code);
	bool checkIntegrity();
	bool checkMu();
	void generateRandomSolution();
	void restartSolution();
	void oneImprovement(int vertex, int code);
	bool swap1_1(int code);
	bool swap2_2(int code);
	bool addFirstVertex();
	bool addBestVertex();
	void VND(int K);
	void shake(double z);
	void greedyRandomizedConstructive(double p);
	void rclConstruction(int code, double p);
	void balanceBiclique();
	void printSolution();
};