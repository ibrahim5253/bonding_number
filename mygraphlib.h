#ifndef MYGRAPHLIB_H
#define MYGRAPHLIB_H

#include<bits/stdc++.h>
#include<igraph.h>

using std::vector;
using std::set;
using std::string;

bool is_bonded(const vector<vector<int>>& G, const set<int>& S, int d);

int max_bonding(const vector<vector<int>>& G, const set<int>& S);

void read_graph(igraph_t*, char*, string, int);

void do_clustering(igraph_t*, igraph_vector_t*, const string&);

bool is_upper_bound(const vector<vector<int>>&, int);

int get_upper_bound(const vector<vector<int>>&, int, int);

int upper_bound_via_shortest_paths(const vector<vector<int>>&);

void generate_barabasi_albert(vector<vector<int>>&, int, int, int);

#include "mygraphlib.cpp"
#endif
