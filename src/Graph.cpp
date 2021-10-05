#include <iostream>
#include <vector>
#include <algorithm> 
#include "Graph.hpp"

using namespace std;

Graph::Graph(int v, int e) {
	this->v = v;
	this->e = e;
	adjList.resize(v);
	weight.resize(v);
	h_index.resize(v);
	accumulatedSum.resize(v);
}

Graph::Graph(Graph &graph) { // copy constructor
	v = graph.getV();
	e = graph.getE();

	adjList.resize(v);
	adjList = graph.get_adjList();

	weight.resize(v);
	weight = graph.get_weight_list();

	h_index.resize(v);
	h_index = graph.getHIndex();

	accumulatedSum.resize(v);
	accumulatedSum = graph.getAccumulatedSum();
}

int Graph::getV() {
	return v;
}

int Graph::getE() {
	return e;
}

int Graph::repeatedEdge(int u, int t) {
	for(unsigned int i = 0; i < adjList[u].size(); i++) {
		if(adjList[u][i] == t) { return 1; }
	}

	return 0;
}

int Graph::get_weight(int u) {
	return weight[u];
}

int Graph::getVertexHIndex(int u) {
	return h_index[u];
}

vector<int> &Graph::get_vertex_accumulatedSum(int u) {
	return accumulatedSum[u];
}

vector<int> &Graph::get_weight_list() {
	return weight;
}

vector<int> &Graph::get_vertex_adjList(int u) {
	return adjList[u];
}

vector<vector<int>> &Graph::get_adjList() {
	return adjList;
}

vector<int> &Graph::getHIndex() {
	return h_index;
}

vector<vector<int>> &Graph::getAccumulatedSum() {
	return accumulatedSum;
}

void Graph::calculateAccumulatedSum(int vertex) { // calculate the accumulated sum for each vertex based on the largest neighbors weights
	vector<int> neighbors = get_vertex_adjList(vertex);
	vector<int> neighbors_weight;
	neighbors_weight.resize(neighbors.size());


	if(!neighbors.empty()) { // vertex has neighbors
		int neighbor_weight, accumulatedVertexSum = 0;
		accumulatedSum[vertex].resize(neighbors.size());

		for(unsigned int idx = 0; idx < neighbors_weight.size(); idx++) {
			neighbor_weight = get_weight(neighbors[idx]);
			neighbors_weight[idx] = neighbor_weight;
		}

		std::sort(neighbors_weight.begin(), neighbors_weight.end(), greater<int>()); // sort neighbors weight in descending order

		for(unsigned int idx = 0; idx < neighbors_weight.size(); idx++) { // start to fill the accumulated sum vector for the designated vertex
			accumulatedVertexSum += neighbors_weight[idx];
			accumulatedSum[vertex][idx] = accumulatedVertexSum; 
		}
	}
	else { // vertex does not have any neighbor
		accumulatedSum[vertex].resize(1);
		accumulatedSum[vertex][0] = 0;
	}
}

void Graph::initializeAccumulatedSum() {
	int V = getV(); 

	for(int idx = 0; idx < V; idx++) { calculateAccumulatedSum(idx); }
}

void Graph::calculateHIndex(int vertex) { // calculate h-index for each vertex
	vector<int> neighbors = get_vertex_adjList(vertex);
	vector<int> hList; // vector used to find the h-index

	if(!neighbors.empty()) {
		unsigned int neighbor_vertex, neighbor_degree, count = 0; 
		hList.resize(neighbors.size() + 1);
		fill(hList.begin(), hList.end(), 0);

		for(unsigned int idx = 0; idx < neighbors.size(); idx++) {
			neighbor_vertex = neighbors[idx];
			neighbor_degree = get_vertex_adjList(neighbor_vertex).size();
			
			if(neighbor_degree >= neighbors.size()) hList[neighbors.size()]++;
			else hList[neighbor_degree]++;
		}

		for(unsigned int idx = neighbors.size(); idx > 0; idx--) {
			count += hList[idx];

			if(count >= idx) { // h-index was found
				h_index[vertex] = idx;
				hList.clear();		
				return;
			}
		}
	}
	else {
		h_index[vertex] = 0;
		hList.clear();
	}
}

void Graph::initializeHIndex() {
	int V = getV(); 

	for(int idx = 0; idx < V; idx++) { calculateHIndex(idx); }
}

void Graph::removeVertexFromAdjList(int vertex, int position) {
	adjList[vertex].erase(adjList[vertex].begin() + position); 
}

void Graph::clearVertexAdjList(int vertex) {
	adjList[vertex].clear();
}

void Graph::addEdges(int u, int t) {
	adjList[u].push_back(t);
	adjList[t].push_back(u);
	return;
}

void Graph::removeEdges(int u, int t) {
	return;
}

void Graph::sort() {
	for(int idx = 0; idx < v; idx++) {
		std::sort(adjList[idx].begin(), adjList[idx].end());
	}
}

void Graph::readEdges() {
	int u, t;
	for(int idx = 0; idx < e; idx++) {
		cin >> u >> t;

		if(u != t && !repeatedEdge(u - 1, t - 1)) { addEdges(u - 1, t - 1); }
	}
}

void Graph::readWeight() {
	for(int idx = 0; idx < v; idx++) {
		cin >> weight[idx];
	}
}

void Graph::showGraphInformations() {
	cout << "Total Vertices: " << getV() << " " << "Total Edges: " << getE() << endl;
	showWeight();
	showAdjList();
}

void Graph::showAdjList() {
	for(int idx = 0; idx < v; idx++) {
		cout << "Vertex: " << idx << endl;
		for(unsigned int idy = 0; idy < adjList[idx].size(); idy++) {
			cout << "Edge: " << idx << " -> " << adjList[idx][idy] << endl; 
		}
		cout << endl;
	}
}

void Graph::showWeight() {
	for(int idx = 0; idx < v; idx++) {
		cout << "Vertex: " << idx << " Weight: " << weight[idx] << endl;
	}
	cout << endl;
}