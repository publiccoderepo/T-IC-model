#ifndef DTARRS_H
#define DTARRS_H

#include <unordered_set>
#include <chrono>
#include <random>
#include <queue>
#include <unordered_map>
#include <iostream>
#include "DiGraph.h"
#include "Hypergraph.h"
#include <list>
using namespace std;

class DTARRS {
public:
    DTARRS(DiGraph & GT);
    virtual ~DTARRS();

    uniform_real_distribution<double> * real_rand;
    uniform_int_distribution<int> * dice_rand;
    mt19937 * mt_rand;

    vector<int> vtxarr;

    unordered_set<int> visited;
    queue<int> Q;

    DiGraph & GT;
    Hypergraph H;

    unordered_set<int> * generateRRS(int t1, int t2);
    void generateHypergraph(int t1, int t2, int N);
//    map<int, unordered_set<int>>  *computeMaxSpread(int timwindowIncre, int largestTimewindow);
    list<unordered_set<int>>  computeMaxSpread(int timwindowIncre, int largestTimewindow);
//    map<int, unordered_set<int>*> *computeTopK(int timwindowIncre, int largestTimewindow);
    list<unordered_set<int>>  computeTopK(int timwindowIncre, int largestTimewindow);
    list<unordered_set<int>>  computeTopK(int timwindowIncre, int largestTimewindow, string summaryResult);
};

#endif /* DTARRS_H */

