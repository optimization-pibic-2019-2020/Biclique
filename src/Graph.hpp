#include <bits/stdc++.h>

using namespace std;

class Graph {
private:	
	int v; // number of vertices
	int e; // number of edges
	vector<int> weight; // vector that contains the vertice´s weight
	vector<vector<int>> adjList; // contains the all the edges

public:
	Graph(int v, int e);
	int getV();
	int getE();
	int get_weight(int u);
	vector<int> &get_weight_list();
	vector<int> &get_vertex_adjList(int u);
	vector<vector<int>> &get_adjList();
	void showGraphInformations();
	void showAdjList();
	void showWeight();
	void readEdges();
	void readWeight();
	void addEdges(int u, int t);
	void removeEdges(int u, int t);
	void sort();

};