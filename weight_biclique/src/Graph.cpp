#include <bits/stdc++.h>
#include "Graph.hpp"


using namespace std;

Graph::Graph(int v, int e) {
	this->v = v;
	this->e = e;
	adjList.resize(v);
	weight.resize(v);
}

int Graph::getV() {
	return v;
}

int Graph::getE() {
	return e;
}

int Graph::get_weight(int u) {
	return weight[u];
}

vector<int> Graph::get_vertex_adjList(int u) {
	return adjList[u];
}

vector<vector<int>> Graph::get_adjList() {
	return adjList;
}

void Graph::readEdges() {
	int u, t;
	for(int idx = 0; idx < e; idx++) {
		cin >> u >> t;
		addEdges(u, t);
	}
}

void Graph::readWeight() {
	int u_weight;
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
		for(int idy = 0; idy < adjList[idx].size(); idy++) {
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