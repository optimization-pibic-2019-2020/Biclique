#ifndef SOLUTION_HPP
#define SOLUTION_HPP
    
#include "Solution.hpp"

#endif



class NonBipartiteSolution: public Solution {
private:
	// weight of each vertex i minus the sum of the weights of its neighbors that are in each solution

	vector<int> mu_A;
	vector<int> mu_B;

public:
    NonBipartiteSolution(Graph *graph, int partitionA_size, int partitionB_size);
	int getRemovedVertices();
    void checkFreePartition();
    void checkNonFreePartition();
    void removeVertex(int u, int code);
	void addVertex(int u, int code);
	bool checkMu();
    void restartSolution(vector<bool> &vertexInGraph);
	void swapVertices(int vertex, int code);
	bool swap1_1(int code);
	bool swap2_2(int code);
	bool addPairOfVertices();
    void VND(int K);
	void rclConstruction(int code, double p);
    void greedyRandomizedConstructive(double p);
	void reduceGraph(vector<bool> &vertexInGraph, int bestWeight, int minBicliqueWeight);
	void balanceBiclique();
};

