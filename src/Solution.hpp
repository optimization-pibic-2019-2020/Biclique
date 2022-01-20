#include <vector>

#ifndef GRAPH_HPP
#define GRAPH_HPP
    
#include "Graph.hpp"

#endif

using namespace std;

class Solution {
protected:
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

	// size of the removed vertices partition 
	int removed_size_A;
	int removed_size_B;

	// total of removed edges from the graph
	int removed_edges;

	// for each vertex, the number of adjacent vertices that are in each solution

	vector<int> tightness_A;
	vector<int> tightness_B;

	// position of each vertex in each solution_ vector

	vector<int> position_A;
	vector<int> position_B;
	
	// current biclique weight
	
	int total_weight;

	// vector of the restricted candidate list for the greedy randomized construction
	
	vector<int> rclList;

	// vector that represents the probability of each element in the rcl

	vector<double> rclListProbability;

public:
	int getTotalWeight();
	int getRemovedEdges();
	bool isNeighbor(int vertex1, int vertex2);
	bool sameNeighbor(int vertex1, int vertex2, int code);
	bool checkBicliqueSize();
	void moveFreeToSolutionPartition(int u, int code);
	void moveFreeToNonFreePartition(int u, int code);
 	void moveSolutionToFreePartition(int u, int code);
	void moveNonFreeToFreePartition(int u, int code); 
	void moveVertexToRemovedVertices(int u, int code);
	bool isMaximal(int code);
	void addRandomVertex(int code);
	bool checkIntegrity();
	bool addBestVertex();
	void createRclProbability();
	void removeVertexFromGraph(int vertex);
	int predictBicliqueWeight(int vertex);
	void printSolution();
	void checkFreePartition();
	void checkNonFreePartition();
	void addVertex(int u, int code);
	void removeVertex(int u, int code);
	void restartSolution(vector<bool> &vertexInGraph);
	void swapVertices(int vertex, int code);
	bool swap1_1(int code);
	bool swap2_2(int code);
	void VND(int K);
	void greedyRandomizedConstructive(double p);
	void reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight);
	void restartAm(double beta);
	void updateAm();
	void balanceBiclique();
};