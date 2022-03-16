
#ifndef DIGRAPH_H
#define DIGRAPH_H

#include "unordered_map"
#include "unordered_set"
#include <utility>     
#include <vector>
using namespace std;

class DiGraph {
public:
    DiGraph();
    DiGraph(const DiGraph& orig);
    virtual ~DiGraph();

    void addEdge(int u, int v, vector<double> w);

    bool existEdge(int u, int v);

    vector<double> returnWeight(int u, int v);

    void print();

    unordered_map<int, vector<double>>& operator[](int u) {
        return adj[u];
    }

    int nvtx, nedge;

    unordered_set<int> V;
    unordered_map<int, unordered_map<int, vector<double>> > adj;
    
};

#endif /* DIGRAPH_H */

