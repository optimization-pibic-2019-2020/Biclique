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
	void checkSolutionPartition();
    void checkFreePartition();
    void checkNonFreePartition();
    void removeVertex(int u, int code);
	void addVertex(int u, int code, int iteration);
	bool checkMu();
    void restartSolution(vector<bool> &vertexInGraph);
	void swapVertices(int vertex_to_add, int vertex_to_remove, int code, int iteration);
	bool swap1_1(int code, int iteration);
	bool swap2_2(int code);
	bool swap1_k(int code, int iteration);
	bool addPairOfVertices(int iteration);
	void randomConstructive();
	void rclConstruction(int code, double alpha, int iteration);
	void greedyRandomizedConstructive(double p, int iteration);
	void moveVertexToNonFreePartition(int u, int code);
	void reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight, double timeLimit);
	void perturb();
	void balanceBiclique();
};