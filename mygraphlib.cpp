#include<bits/stdc++.h>
#include"mygraphlib.h"
#include<igraph.h>

#define pb push_back
#define mp make_pair
#define ff first
#define ss second
using namespace std;
using ll=long long;
using pii=pair<int, int>;
using pll=pair<ll,ll>;
using vi=vector<int>;
using vll=vector<ll>;
using vvi=vector<vi>;

bool is_bonded(const vector<vector<int>>& G, const set<int>& S, int d)
{
    int maxD = max_bonding(G, S);
    return d<=maxD;
}

int rand_partition_lb(const vvi& G)
{
    set<int> S;
    while(1) {
        for (int i=0; i<G.size(); ++i)
            if (rand()%2 == 0) S.insert(i);
        if(not S.empty()) break;
    }
    return max_bonding(G,S)+1;
}

int max_bonding(const vector<vector<int>>& G, const set<int>& S)
{
    int ans=INT_MAX;
    for (int i=0, n=G.size(); i<n; ++i) {
        int deg=0;
        if (S.count(i)) {
            for (auto& j: G[i])
                if (S.count(j)) ++deg;
        }
        else {
            for (auto& j: G[i])
                if (!S.count(j)) ++deg;
        }
        if (deg < G[i].size()) ans = min(ans, deg);
    }
    return ans;
}

void read_graph(igraph_t* g, char* fname, string format, int directed)
{
    FILE* fp = fopen(fname, "r");
    if (format=="gml")
        igraph_read_graph_gml(g, fp);
    else if (format=="edgelist")
        igraph_read_graph_edgelist(g, fp, 0, directed);
    fclose(fp);
}

void do_clustering(igraph_t* g, igraph_vector_t* membership, const string& method)
{
    igraph_vector_t modularity;
    igraph_matrix_t merges;

    igraph_vector_init(&modularity,0);
    igraph_matrix_init(&merges,0,0);

    int n = igraph_vcount(g);

    if (method=="optimal") 
        igraph_community_optimal_modularity(g, 0, membership, 0);

    else if (method=="greedy") {
        igraph_community_fastgreedy(g, NULL, &merges, &modularity, membership);
        long i=igraph_vector_which_max(&modularity);
        igraph_community_to_membership(&merges, n, i, membership, 0);
    }

    else if (method=="multilevel") 
        igraph_community_multilevel(g, NULL, membership, NULL, NULL);

    igraph_vector_destroy(&modularity);
    igraph_matrix_destroy(&merges);
}

void dfs(const vvi& G, int x, vector<int>& vis, vi& label, vi& cc, int k, vvi& reverse_cc, int dp=1)
{
   //cout<<"In dfs. x: "<<x<<" d: "<<dp<<endl; 
   vis[x]=true;
   cc[x]=k;
   reverse_cc[k].pb(x);
   //cout<<"Assigned."<<endl;
   
   for (auto& i: G[x]) {
       if (!vis[i] and label[i] and label[i]+label[x]<4)
           /*cout<<i<<endl,*/ dfs(G, i, vis, label, cc, k, reverse_cc, dp+1);
   }
}

pair<bool, int> is_upper_bound(const vvi& G, int d)
{
    int n = G.size();
    dsu D(n);

    for (int i=0; i<n; ++i) {
        if (G[i].size() <= d)
            for (auto& j: G[i])
                D.union_(i, j);
    }

    if (D.num_comp() == 0) return {false, -1};

    list<int> T;
    for (int i=0; i<n; ++i)
        if (D.comp_size(i) == 1) T.pb(i);

    while (1) {
        bool f = false;
        for (auto it = T.begin(); it != T.end(); ){
            auto i = *it;
            map<int, int> cnt;
            for (auto& j: G[i])
                ++cnt[D.find_set(j)];
            int k = G[i].size()-d+1;
            for (auto& p: cnt)
                if (p.ss >= k) D.union_(i, p.ff);
            if (D.comp_size(i) > 1) it = T.erase(it), f = true;
            else ++it;
        }
        if (not f) break;
    }

    if (D.num_comp() == 1) return {true, 0};

    set<int> conn_comp;
    for (int i=0; i<n; ++i) conn_comp.insert(D.find_set(i));

    int k = conn_comp.size();
    cout<<"Number of components: "<<k<<endl;
    if (k > 60) return {false, 0};

    vi v(conn_comp.begin(), conn_comp.end());

    int M = 0;
    for (int i=1; i<(1ll<<k)-1; ++i) {
        set<int> s;
        for (int j=0; j<k; ++j)
            if ( (i>>j)&1 ) s.insert(v[j]); 
        M = max(M, max_bonding_dsu(G, s, D));
    }

    if (M < d) return {true, M};
    else return {false, M};
}

int get_upper_bound(const vvi& G, int l, int r, int& lb)
{
    while (l<r) {
        int m = l+(r-l)/2;
        auto p = is_upper_bound(G, m);
        lb = max(lb, p.ss);
        if (p.ff)
            r = m;
        else
            l = m+1;
    }
    return l;
}

void compute_weights(const vvi& G, vvi& D)
{
    for (int i=0, j=G.size(); i<j; ++i)
        for (auto& x: G[i])
            D[i][x] = min(G[i].size(), G[x].size());
}

int upper_bound_via_shortest_paths(const vvi& G)
{
    int n=G.size();
    vvi D(n, vi(n, 100));
    for (int i=0; i<n; ++i)
        D[i][i]=0;

    compute_weights(G, D);

    for (int k=0; k<n; ++k)
        for (int i=0; i<n; ++i)
            for (int j=0; j<n; ++j)
                D[i][j] = min(D[i][j], max(D[i][k], D[k][j]));

    /*
    for (int i=0; i<n; ++i) {
        cout<<i<<": ";
        for (int j=0; j<n; ++j)
            cout<<D[i][j]<<" ";
        cout<<endl;
    }
    */

    int ans=0;
    for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
            ans = max(ans, D[i][j]);

    return ans;
}

void generate_barabasi_albert(vvi &G, int n, int m, int ne)
{
    if (n==0) return;
    vi belongs(2*ne, -1);

    int j=0;
    for (int i=0; i<G.size(); ++i) {
        int d = G[i].size();
        while (d--) belongs[j++] = i;
    }

    int idx = G.size();
    
    //vector<bool> added(G.size(), false);    
    vector<int> added(G.size(), 0);

    vi edges;

    for (int i=0; i<m; ++i) {
        int x = rand() % (2*ne);
        while (added[belongs[x]]) x = rand() % (2*ne);
        edges.pb(belongs[x]), added[belongs[x]] = 1, G[belongs[x]].pb(idx);
    }

    G.pb(edges);

    //cout<<G.size()<<" "<<get_upper_bound(G, 1, G.size())/pow(G.size(),0.5)<<"\n";

    generate_barabasi_albert(G, n-1, m, ne+m);
}

int max_bonding_dsu(const vvi& G, const set<int>& S, dsu& D)
{
    int ans=INT_MAX;
    for (int i=0, n=G.size(); i<n; ++i) {
        int deg=0;
        if (S.count(D.find_set(i))) {
            for (auto& j: G[i])
                if (S.count(D.find_set(j))) ++deg;
        }
        else {
            for (auto& j: G[i])
                if (not S.count(D.find_set(j))) ++deg;
        }
        if (deg < G[i].size()) ans = min(ans, deg);
    }
    return ans;
}


