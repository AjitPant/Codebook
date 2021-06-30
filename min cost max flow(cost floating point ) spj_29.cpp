const double INF = 1e9;
struct MinimumCostMaximumFlow 
{
    typedef int Flow;
    typedef double Cost;
    static constexpr Cost infiniteDistance=1e90;
    static constexpr Cost EPS=1e-7;
    static const Flow infiniteFlow=1e9;
    struct Edge
    {
        int u,v;
        Flow f,c;
        Cost w;
        Edge(int u,int v,Flow f,Flow c,Cost w):u(u),v(v),f(f),c(c),w(w){}
    };  
    vector<Edge> e;
    vector<vector<int> > g;  
    int n,source,sink,*prev;
    Cost *dist;
    MinimumCostMaximumFlow(int n):n(n)
    {
        dist=(Cost*)malloc(sizeof(Cost)*n);
        prev=(int*)malloc(sizeof(int)*n);
        g.resize(n);
    }
    ~MinimumCostMaximumFlow()
    {
        free(dist);
        free(prev);
        g.clear();
    }
    void add(int u,int v,Flow c,Cost w)
    {
        g[u].push_back(e.size());
        e.push_back(Edge(u,v,0,c,w));
        g[v].push_back(e.size());
        e.push_back(Edge(v,u,0,0,-w));
    }   
    pair<Cost,Flow> minimumCostMaximumFlow(int source,int sink)
    {
        this->source=source;
        this->sink=sink;
        Flow flow=0;
        Cost cost=0;
        while(bellmanFord())
        {
            int u=sink;
            Flow pushed=infiniteFlow;
            Cost pushCost=0;
            while(u!=source)
            {   
                int id=prev[u];
                pushed=min(pushed,e[id].c-e[id].f);
                pushCost+=e[id].w;
                u=e[id].u;
            }
            u=sink;
            while(u!=source)
            {
                int id=prev[u];
                e[id].f+=pushed;
                e[id^1].f-=pushed;   
                u=e[id].u;
            }
            flow+=pushed;
            cost+=pushCost*pushed;
        }
 
        return make_pair(cost,flow);
    }
    bool bellmanFord()
    {
        for(int i=0;i<n;++i)
            dist[i]=infiniteDistance;
        dist[source]=0;
 
        for(int k=0;k<n;++k)
        {
            bool update=false;
            for(int id=0;id<(int)e.size();++id)
            {
                int u=e[id].u;
                int v=e[id].v;
                if(dist[u]+EPS>=infiniteDistance)
                    continue;
                Cost w=e[id].w;
                if(e[id].f<e[id].c && dist[v]>dist[u]+w+EPS)
                {
                    dist[v]=dist[u]+w;
                    prev[v]=id;
                    update=true;
                }   
            }
            if(!update)
                break;
        }
        return dist[sink]+EPS<infiniteDistance;
    }
    void displayEdges()
    {
        cout<<"******"<<endl;
        for(int i=0;i<(int)e.size();i+=2)
            cout<<e[i].u<<" "<<e[i].v<<" "<<e[i].f<<" "<<e[i].c<<" "<<e[i].w<<"\n";
        cout<<"******"<<endl;
    }
};
