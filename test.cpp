#include "mygraphlib.h"
#include <bits/stdc++.h>
#include <igraph.h>

using namespace std;

void makeset(const vector<set<int>>& v, set<int>& S)
{
    for (auto& x: v)
        S.insert(x.begin(), x.end());
}

int main(int argc, char** argv)
{
	ios::sync_with_stdio(false);
    igraph_t g;

    read_graph(&g, argv[1], argv[2], atoi(argv[3]));

    igraph_bool_t res;
    igraph_is_connected(&g, &res, IGRAPH_STRONG);
    cout<<"Is connected: "<<res<<endl;
    
    int n = igraph_vcount(&g);
    
    
    igraph_vector_t membership;
    igraph_vector_init(&membership, n);
    
    do_clustering(&g, &membership, (argc>4?argv[4]:"greedy"));

    cout<<"Clustering done."<<endl;

    int n_clusters = 0;
    vector<set<int>> cluster(igraph_vcount(&g));
    for (int i=0; i<igraph_vector_size(&membership); ++i) {
        cluster[VECTOR(membership)[i]].insert(i);
        if (VECTOR(membership)[i] >= n_clusters) ++n_clusters;
    }

    cout<<"Number of clusters: "<<n_clusters<<endl;
    int nv=0;
    for (int i=0; i<n_clusters; ++i)
        cout<<"Cluster "<<i<<": "<<cluster[i].size()<<endl, nv+=cluster[i].size();
    cout<<"N="<<nv<<endl;
    
    int max_deg=0, min_deg=n+1;

    vector<vector<int>> G(n);
    for (int i=0; i<n; ++i) {
        igraph_vector_t neigh;
        igraph_vector_init(&neigh, 0);
        igraph_neighbors(&g, &neigh, i, IGRAPH_ALL);
        int num_neighbors= igraph_vector_size(&neigh);
        G[i].reserve(num_neighbors);
        max_deg = max(max_deg, num_neighbors);
        min_deg = min(min_deg, num_neighbors);
        for (int j=0; j<igraph_vector_size(&neigh); ++j) {
            G[i].push_back(VECTOR(neigh)[j]);
        }
        igraph_vector_destroy(&neigh);
    }

   cout<<"Graph converted to adjacency list format."<<endl;
   cout<<"Max degree: "<<max_deg<<endl;

    igraph_destroy(&g);
    igraph_vector_destroy(&membership);

    int d1=0;
    for (int i=0; i<n_clusters; ++i)
        d1 = max(d1, max_bonding(G, cluster[i]));
    cout<<"Lower bound(based on clustering, original): "<<d1+1<<endl;
    
    /*
    int d=0;
    srand(time(NULL));
    for (int iter=1; iter < (1<<n_clusters)-1; ++iter) {
        set<int> S;
        for (int i=0; i<n_clusters; ++i)
            if ( (iter>>i)&1 ) S.insert(cluster[i].begin(), cluster[i].end());
        if (S.empty()) continue;
        d = max(d, max_bonding(G, S));
        cout<<iter<<" "<<d<<endl;
    }
        

    cout<<"Lower bound(based on clustering): "<<d+1<<endl;
    */

/*    int d2=0;
    for (int i=0; i<1000; ++i) 
        d2 = max(d2, rand_partition_lb(G));
   
    cout<<"Lower bound(based on rand partition): "<<d2<<endl;*/
    cout<<"Upper bound: "<<is_upper_bound(G,3)<<endl;//get_upper_bound(G, 1, max_deg)<<endl;

    set<int> v;
    set<pii> ed;
    for (int i=0; i<G.size(); ++i) {
        bool f=false;
        for (auto &x : G[i]) if (i<x and (G[i].size() <= 3 or G[x].size() <= 3)){
            v.insert(i+1), v.insert(x+1), ed.insert({i+1, x+1});
        }
    }
    for (auto& i: v) cout<<i<<"\n";
    for (auto& p: ed) cout<<p.ff<<" "<<p.ss<<"\n";

//    cout<<"Upper bound via shortest paths: "<<upper_bound_via_shortest_paths(G)<<endl;

	return 0;
}
