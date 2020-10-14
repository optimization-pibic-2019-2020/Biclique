#include <bits/stdc++.h>

using namespace std;

int main() {
    long long int edges = 0, partition1_size = -1, partition2_size = -1, v1, v2;
    char virgula;

    while(cin >> v1 >> virgula >> v2) {
        if(v1 > partition1_size) partition1_size = v1;
        if(v2 > partition2_size) partition2_size = v2;
        edges++;
    }

    cout << "Partition 1: " << partition1_size << endl;
    cout << "Partition 2: " << partition2_size << endl;
    cout << "Edges: " << edges << endl;
    return 0;
}