#ifndef SOLUTION_HPP
#define SOLUTION_HPP
    
#include "Solution.hpp"

#endif

class BipartiteSolution: public Solution {
private:
	// size of each partition 
	int partition_size_A;
	int partition_size_B; 

public:
 	BipartiteSolution(Graph *graph, int partitionA_size, int partitionB_size);
	int getRemovedVertices();
    void checkFreePartition();
    void checkNonFreePartition();
    void removeVertex(int u, int code);
	void addVertex(int u, int code);
	bool checkMu();
    void restartSolution(vector<bool> &vertexInGraph);
	void swapVertices(int vertex_to_add, int vertex_to_remove, int code);
	bool swap1_1(int code);
	bool swap2_2(int code);
	bool swap1_k(int code);
	bool addPairOfVertices();
    void VND(int K);
	void rclConstruction(int code, double p);
    void greedyRandomizedConstructive(double p);
	void reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight, double timeLimit);
	void balanceBiclique();
};