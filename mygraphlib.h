#ifndef MYGRAPHLIB_H
#define MYGRAPHLIB_H

#include<bits/stdc++.h>
#include<igraph.h>

using std::vector;
using std::set;
using std::string;
using std::pair;
using vi=vector<int>;
using vvi=vector<vi>;

bool is_bonded(const vector<vector<int>>& G, const set<int>& S, int d);

int max_bonding(const vector<vector<int>>& G, const set<int>& S);

void read_graph(igraph_t*, char*, string, int);

void do_clustering(igraph_t*, igraph_vector_t*, const string&);

pair<bool, int> is_upper_bound(const vector<vector<int>>&, int);

int get_upper_bound(const vector<vector<int>>&, int, int);

int upper_bound_via_shortest_paths(const vector<vector<int>>&);

void generate_barabasi_albert(vector<vector<int>>&, int, int, int);

class dsu {
    private:
        vi par, rank, csize;
        int nc=0;
    public:
        dsu(int n) : par(n), rank(n, 0), csize(n, 1), nc(n) {
            iota(par.begin(), par.end(), 0);
        }  
        int find_set(int x)
        {
            if (x != par[x])
                par[x] = find_set(par[x]);
            return par[x];
        }
        void union_(int x, int y) {link(find_set(x), find_set(y));}
        void link(int x, int y) {
            if (x==y) return;
            if (rank[x] > rank[y])
                par[y] = x, csize[x] += csize[y];
            else {
                par[x] = y, csize[y] += csize[x];
                if (rank[x] == rank[y])
                    ++rank[y];
            }
            --nc;
        }
        int comp_size(int x) {
            if (x != par[x])
                par[x] = find_set(par[x]);
            return csize[par[x]];
        }
        int num_comp() {return nc;}
};

int max_bonding_dsu(const vvi& G, const set<int>& s, dsu& D);
    

#include "mygraphlib.cpp"
#endif
