#include <vector>
#include <unordered_map>
#include "PCenter.h"

using namespace std;

class BaseCircle {

    PCenter *thePc;

    struct Node {
        int id;
        vector<int> from;
        vector<int> to;

        explicit Node(int nodeId) : id(nodeId) {}
    };

    struct Edge {
        int edgeId;
        int from;
        int to;

        Edge(int edgeId, int from, int to) : edgeId(edgeId), from(from), to(to) {}

        Edge(int from, int to) : edgeId(-100), from(from), to(to) {}

        bool operator==(const Edge &e) const {
            return from == e.from && to == e.to;
        }
    };

    struct hash_name {
        size_t operator()(const Edge *e) const {
            return hash<int>()(e->from) * 10000 + hash<int>()(e->to);
        }
    };

    struct edge_equal {
        bool operator()(const Edge *e1, const Edge *e2) const {
            return (e1->from == e2->from && e1->to == e2->to);
        }
    };

public:
    void init(int N, vector<vector<int>> adj);

    vector<int> shouldAdd;
    vector<int> shouldErase;
    vector<int> detectingNodes;
    vector<int> remainNodes;

    vector<int> inDegree;
    vector<int> visit;

    string inputFile;

    int N, M, Type;

    vector<unordered_map<int, int> > connect;
    vector<Node *> nodeVec;
    vector<Node *> subNodeVec;

    vector<Edge *> edges;
    unordered_map<Edge *, int, hash_name, edge_equal> getEdgeIdByNodes;
    unordered_map<int, Edge *> getTwoNodesByEdgeId;

    vector<vector<int>> cycleVector;

    vector<vector<int>> BasicCycleRes;
    int addedCycleCount = 0;

    vector<bool> lastSol;
    vector<int> lastCenters;

    vector<int> father;

    void dfsVisit(int nodeId);

    vector<vector<int>> tmpCycleVec;

    unordered_map<int, vector<int>> addMap;

public:
    void read(string inputFile);

    void dfs4Cycle(int level);

    void buildSubGraph(vector<bool> &firstSolution);

    void findCycleByDFS();

    void addedDeleteIncludeCycle(vector<vector<int>> cycles, vector<vector<int>> &originalCycles, bool isBasic = false);

    void updateSubGraph(vector<bool> firstSolution, vector<int> centers);

    // 随机数 上限时间 初始反馈集
    vector<int> runFVSP(int seed, int limitTime, vector<int> upperBoundNodes);

    void updateShouldAdd(vector<bool> &fvsSolution, vector<int> &nowCenters);

    int dfsDetect(int nodeId);
};