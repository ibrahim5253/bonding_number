#include <bits/stdc++.h>
#include <igraph.h>
#include "mygraphlib.h"

using namespace std;

using vvi=vector<vector<int>>;

vvi g;

int n, m;
int num_steps = 100, num_edges = 5;

int main(int argc, char** argv)
{
    num_steps = atoi(argv[1]), num_edges = atoi(argv[2]);
    cin>>n>>m;
    g.resize(n);

    for (int i=0; i<m; ++i) {
        int u, v;cin>>u>>v;
        g[u].pb(v), g[v].pb(u);
    }
    auto seed = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now().time_since_epoch()).count();
    srand(seed);

    generate_barabasi_albert(g, num_steps, num_edges, m);

    long long max_deg=0;

    for (int i=0; i<g.size(); ++i)
        max_deg = max(max_deg, (long long)g[i].size());

    cout<<g.size()<<" "<<get_upper_bound(g, 1, g.size())<<" "<<max_deg<<" "
        <<pow(max_deg+m,3)/(2.0*g.size()*m*(m+1))<<"\n";



    //cout<<"Graph Size: "<<g.size()<<endl;
    /*for (int i=0; i<g.size(); ++i) {
        for (auto& x : g[i])
            if (x>i) cout<<i<<" "<<x<<"\n";
    }*/

    /*
    int min_degree = n+101;

    for (auto& v: g)
        min_degree = min(min_degree, (int)v.size());*/

    //cout << "Lower bound (based on min degree): " << (min_degree-1)/2 <<endl;
    //cout << "Upper bound: "<< get_upper_bound(g, 0, n+101) << endl;
    return 0;
}
