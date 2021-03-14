#include <bits/stdc++.h>

using namespace std;


int main(int argc, char* argv[]) {
    int V = 1000, E = 0, chosenNumber, edgePercentage = 1;
    edgePercentage = atoi(argv[1]);
    vector<pair<int, int>> edges;
    srand((unsigned) time(0));

    for(int i = 0; i < V; i++) {
        for(int j = 0; j < V; j++) {
            chosenNumber = rand() % 10;
            if(i != j && chosenNumber < edgePercentage) {
                edges.push_back(make_pair(i + 1, j + 1));
                E++;
            }
        }
    }

    cout << V << " " << E << endl;

    for(int i = 0; i < V; i++) {
        cout << (i % 200) + 1 << endl;
    }

    for(int i = 0; i < edges.size(); i++) {
        cout << edges[i].first << " " << edges[i].second << endl;
    }

    return 0;
}