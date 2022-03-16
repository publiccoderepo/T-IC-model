

#ifndef DTIC_H
#define DTIC_H

#include <pthread.h>
#include "DiGraph.h"
#include "Hypergraph.h"
#include "DTARRS.h"

#include <chrono>
#include <random>
#include <cstdlib>
#include <vector>
#include <queue>
#include <iostream>                  
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;

double DTIC(int t1, int t2, unordered_set<int> & A, DiGraph & G, int iterNum);
int DTIC2(int t1, int t2, unordered_set<int> & S, unordered_set<int> & T, DiGraph & G);
void DTIC_P(int t1, int t2, unordered_set<int> * A, DiGraph * G, int nthreads, int nIters);
vector<int> DTIC3(int t1, int t2, unordered_set<int> & SI, unordered_set<int> & TI1, unordered_set<int> & TI2, DiGraph & G);
vector<int> new_DTIC3(int t1, int t2, unordered_set<int> & SI, unordered_set<int> & TI1, unordered_set<int> & TI2, unordered_set<int> & TI3, map<int, unordered_set<int>> &IterrandomMethod,DiGraph & G);
//vector<int> new_DTIC3(int t1, int t2, unordered_set<int> & SI, unordered_set<int> & TI1, unordered_set<int> & TI2, unordered_set<int> & TI3, unordered_set<int> & TI4, DiGraph & G);
vector<int> new_DTIC4(int t1, int t2, unordered_set<int> & SI, unordered_set<int> & TI1, DiGraph & G);
map<int,unordered_set<int>> greedy_solution_set(int nodeNum, int selectK, int batchNum, int sampleNum, int t1, int t2, DiGraph & G);
double greedIC(unordered_set<int> selectSeed, int candidateNode, int iterNum, int t1, int t2, DiGraph & G);
vector<int> allbaseline_DTIC3(int t1, int t2,unordered_set<int> & seed_S, unordered_set<int> & rsm, unordered_set<int> & esm, unordered_set<int> & maxd, unordered_set<int> & greedyone, map<int, unordered_set<int>> &IterrandomMethod,DiGraph & G);
vector<int> dynamic_TIMMU_DTIC(int t1, int t2,unordered_set<int> & seed_S, unordered_set<int> & TIMMU_set, map<int, unordered_set<int>>& DYNAMIC_set,DiGraph & G);
map<int, map<int,unordered_set<int>>> longgreedy_solution_set(int nodeNum, int selectK, int batchNum, int sampleNum, int t1, int t2, DiGraph & G);
vector<int> greedy_DTIC(int t1, int t2,unordered_set<int> & seed_S, unordered_set<int> & rsm, unordered_set<int> & esm, unordered_set<int> & maxd, map<int, unordered_set<int>> & greedyone, map<int, unordered_set<int>> &IterrandomMethod,DiGraph & G);



vector<int> alldynamic_TIMMU_DTIC(int t1, int t2,unordered_set<int> & seed_S, unordered_set<int> & TI1, unordered_set<int> & TI2, unordered_set<int> & TI3, map<int, unordered_set<int>> &IterrandomMethod,unordered_set<int> & TIMMU_set, map<int, unordered_set<int>>& DYNAMIC_set,DiGraph & G);

vector<int> allonlydynamic_DTIC(int t1, int t2,unordered_set<int> & seed_S, unordered_set<int> & TI1, unordered_set<int> & TI2, unordered_set<int> & TI3, map<int, unordered_set<int>> &IterrandomMethod, map<int, unordered_set<int>>& DYNAMIC_set,DiGraph & G);


//vector<int> onlydynamic_DTIC(int t1, int t2,unordered_set<int> & seed_S, map<int, unordered_set<int>>& DYNAMIC_set,DiGraph & G);


#endif /* DTIC_H */

