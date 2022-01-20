#include <vector>

using namespace std;

class Graph {
private:	
	int v; // number of vertices
	int e; // number of edges
	vector<int> weight; // vector that contains the vertex weight
	vector<int> am; // vector that contains the value of adaptative memory for each vertex
	vector<int> h_index; // vector that contains the h-index for each vertex
	vector<vector<int>> accumulatedSum; // vector that contais the accumulated sum for each vertex based on the largest neighbors weights
	vector<vector<int>> adjList; // contains the all the edges

public:
	Graph(int v, int e);
	Graph(Graph &graph);
	int getV();
	int getE();
	int repeatedEdge(int u, int t);
	int getVertexAm(int u);
	void setVertexAm(int u, double newValue);
	void initializeAm();
	int get_weight(int u);
	int getVertexHIndex(int u);
	vector<int> &get_vertex_accumulatedSum(int u);
	vector<int> &getAmList();
	vector<int> &get_weight_list();
	vector<int> &get_vertex_adjList(int u);
	vector<vector<int>> &get_adjList();
	vector<int> &getHIndex();
	vector<vector<int>> &getAccumulatedSum();
	void calculateAccumulatedSum(int vertex);
	void initializeAccumulatedSum();
	void calculateHIndex(int vertex);
	void initializeHIndex();
	void removeVertexFromAdjList(int vertex, int position);
	void clearVertexAdjList(int vertex);
	void showGraphInformations();
	void showAdjList();
	void showWeight();
	void readEdges();
	void readWeight();
	void addEdges(int u, int t);
	void removeEdges(int u, int t);
	void sort();

};