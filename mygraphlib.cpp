#include<bits/stdc++.h>
#include"mygraphlib.h"
#include<igraph.h>

#define MOD 1000000007
#define pb push_back
#define mp make_pair
#define umap unordered_map
#define ff first
#define ss second
using namespace std;
using ll=long long;
using ull=unsigned long long;
using pii=pair<int, int>;
using pll=pair<ll,ll>;
using vi=vector<int>;
using vll=vector<ll>;
using pill=pair<int,ll>;
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

bool is_upper_bound(const vvi& G, int d)
{
    int n = G.size();
    vi label(n, 0);

    int r_size=0, max_deg=0;
    for (int i=0; i<n; ++i) {
        if (G[i].size() <= d)
            label[i]=1, ++r_size;
        max_deg = max(max_deg, (int)G[i].size());
    }

    if (r_size == 0) return 0;

    cout<<"R: ";
    for (int i=0; i<n; ++i) 
        if (label[i]==1) cout<<i+1<<" ";
    cout<<endl;

    cout<<"S: ";
    for (int i=0; i<n; ++i) {
        if (label[i]==1) continue;
        for (auto&x : G[i])
            if (label[x]==1) {
                label[i]=2;
                cout<<i+1<<" ";
                break;
            }
    }
    cout<<endl;

//    //cout<<"reached here. R_size: "<<r_size<<endl;

    vi cc(n,0);
    vector<int> vis(n, false);
    vvi v(max_deg+1, vi());
    vvi reverse_cc(1, vi());
    
    int k=0;
    for (int i=0; i<n; ++i) {
        //cout<<"i: "<<i<<endl;
        if (label[i] and !vis[i]) reverse_cc.pb(vi()), dfs(G, i, vis, label, cc, ++k, reverse_cc);
        v[G[i].size()].pb(i);
    }

//    //cout<<"cc complete. #CC: "<<k<<endl;

    while (1) {
        bool f = false;
        for (int i=1; i<=max_deg-d; ++i) 
            for (auto& x: v[i+d]) {
                if (label[x]) continue;
                vi neighbors(k+1,0);
                for (auto& y: G[x])
                    ++neighbors[cc[y]];
                vi tbc;
                for (int j=1; j<=k; ++j)
                    if (neighbors[j]>i) tbc.pb(j);
                if (!tbc.empty()) {
                    cout<<"Adding "<<x<<"..."<<endl;
                    label[x]=2, cc[x]=tbc.front(), f=true, reverse_cc[tbc.front()].pb(x);
                    for (int j=1; j<tbc.size(); ++j) {
                        int comp = tbc[j];
                        for (auto& y: reverse_cc[comp])
                            cc[y]=tbc.front();
                        reverse_cc[tbc.front()].insert(reverse_cc[tbc.front()].end(), 
                                reverse_cc[comp].begin(), reverse_cc[comp].end());
                        reverse_cc[comp].clear();
                    }
                }
            }
        if (!f) break;
    }

    return reverse_cc[1].size() ==  n;
}

int get_upper_bound(const vvi& G, int l, int r)
{
    while (l<r) {
        int m = l+(r-l)/2;
        if (is_upper_bound(G, m))
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

