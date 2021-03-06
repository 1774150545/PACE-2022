#ifndef TOPO_TOPO_H
#define TOPO_TOPO_H

#include "myList.h"

#include <csignal>
#include <random>
#include <unordered_map>
#include "cmath"
#include "time.h"
#include <queue>
#include <stack>
#include <set>
#include <chrono>

typedef pair<int, int> Edge;

using namespace std;

class Config {
private:
    Config() {
        flag = false;
    };
    static Config *cfg;
public:
    bool flag;

    static Config *getCfg() {
        return cfg;
    }

    void setFlag() {
        flag = true;
    }

    bool TLE() {
        return flag;
    }
};

struct Node {
    int nodeId;
    std::vector<int> from; // pred
    std::vector<int> to; // succ


    int inFrom(int theNodeId) {
        for (int i = 0; i < from.size(); i++) {
            if (from[i] == theNodeId) {
                return i;
            }
        }
        return -1;
    }

    int inTo(int theNodeId) {
        for (int i = 0; i < to.size(); i++) {
            if (to[i] == theNodeId) {
                return i;
            }
        }
        return -1;
    }

    // ----------- cache -----------
    int legalNext;
    int legalPre;
    int insertNextNode;
    int insertPreNode;
    int rightLimitNode;
    int leftLimitNode;

    int nextConflictCount;
    int preConflictCount;
    int *nextConflictNodes;
    int *preConflictNodes;

    Node(int nodeId) {
        this->nodeId = nodeId;
        nextConflictNodes = nullptr;
        preConflictNodes = nullptr;
    }

    void init(int fromSz, int toSz) {
        legalNext = false;
        legalPre = false;
        insertNextNode = -1;
        rightLimitNode = -1;
        insertPreNode = -1;
        leftLimitNode = -1;
        nextConflictCount = 0;
        preConflictCount = 0;
        nextConflictNodes = new int[toSz];
        preConflictNodes = new int[fromSz];
    }

    Node() {}

    ~Node() {
        delete[] nextConflictNodes;
        delete[] preConflictNodes;
    }
};

struct oneMoveRes {
    int moveType;
    int conflictCount;
    int baseNode;

    oneMoveRes(int type, int count, int node) {
        moveType = type;
        conflictCount = count;
        baseNode = node;
    }
};

class Topo {
public:
    std::mt19937 pseudoRandNumGen;

    int rand(int lb, int ub) { return std::uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }

    int rand(int ub) { return std::uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

    double rand01Double() { return std::uniform_real_distribution<double>(0, 1)(pseudoRandNumGen); }

    void initRand(int seed) { pseudoRandNumGen = std::mt19937(seed); }

    int originN, N, M, Type;
    int edgeCount;
    double GlobalStartTime;

    // ??????????????????
    std::vector <std::unordered_map<int, int>> connect;

    vector<Node *> nodeVec;

    std::vector<Node *> nodeVecTemp;

    // ------------------------------------------ ??????????????????????????? -------------------------------------------
    int CheckInterval = 100000;
    double CostTimeCheckInterVal = 15;
    int MinBestConfCountChange = 2;
    bool checkBetterPerIter = false; // ?????????????????????????????????
    int EnableCheckTime = 500;
    int addTTimes = 0;

    // ----------------------------------------Jump out of local optimum---------------------------------------
    double delPercent = 0.1;
    double increaseT = 0.1;
    bool enableJumpOut = false;

    myList cList;

    int bestConfCount;
    vector<int> bestUnNumbered; // best FVS

    int *unNumbered;
    int *unNumberedPos;
    int unNumberedCount;

    vector<int> fvSolution;

    void getFVSolution();

    int *conflictNodes;
    int conflictCount;

    pair <lType, lType> getLabelRange(int nodeId);
    // --------------------------------------------------------------------------------------------------------------

    oneMoveRes findOneMove(int chosen, int moveType = -1);

    void doMove(int chosen, int moveType, int insertBaseNode);

    int
    FVSLocalSearch(double T0 = 0.6, double alpha = 0.99, int maxMvt = 10000, int maxFail = 50, bool runOnce = false);

    // ????????????label????????? -1???????????????
    int getNextPosNode(int nodeId);

    // ?????????????????????lable????????? -1???????????????
    int getPrePosNode(int nodeId);

    // ???vector????????????????????????
    void eraseUnNumbered(int node);

    // ???vector????????????????????????
    void pushUnNumbered(int node);

    void init();

    // ---------------------------------------reduction-----------------------------------------------
    bool selfLoop();

    bool clearSelfLoop(); // ??????selfLoop
    int mustInFvsCount;
    vector<int> mustInFVS; // ???FVS?????????
    vector<int> isReduced;

    void reduceVertex(int nodeId, bool selfLoop = false);

    bool doGraphReduction2();

    void reduction();

    // 6. pie
    bool pie();

    vector <Edge> acyclicEdges();

    vector <Edge> piEdges();

    vector<int> vertexToStronglyConnectedComponentNumberNR();

    void dfsNR(int nodeId);

    int tarjanCurIndex;
    vector<int> tarjanIndex;
    vector<int> tarjanAncestor;
    int currentCpmponent;
    vector<int> sccsByNum;

    stack<int> st; // ???
    stack<int> trace; // dsf??????????????????

    // 7.core
    vector<int> core();

    bool isClique(vector<int> &nodes);

    // 8.dome
    bool dome();

    bool isDominate(Edge e);

    void reduceEdge(Edge e);

    //9.???????????????
    template<typename Key, typename Value> using HashMap = std::unordered_map<Key, Value>;

    template<typename T = int>
    T Min(const T &value1, const T &value2) { return std::min(value1, value2); }

    template<typename T = int>
    T Max(const T &value1, const T &value2) { return std::max(value1, value2); }

    template<typename Container>
    void Sort(const Container &value1, const Container &value2) { return sort(value1, value2); }

    static constexpr double process_max_time = 50;
    static constexpr double viaCliqueMaxNodeNum = 550;
    template<typename T> using Vec = std::vector<T>;
    using ID = int;

    void maxCliqueReduce();

    void viaClique(Vec<ID> &ComponentSub);

    bool allIIEdgeJudge(Vec<ID> &ComponentSub);

    // ???????????????nodes
    unordered_map<int, int> originToReducedMap;
    unordered_map<int, int> reducedToOriginMap;


    // ------------------------------------- ????????? --------------------------------------
    int initSolutionSelect = 40000; // ??????????????????????????????
    //double singleRatio = 0.85;
    double singleRatio = 0.6;

    template<typename Container, typename T>
    typename Container::iterator Find(Container &container, const T &value) {
        return std::find(container.begin(), container.end(), value);
    }

    void InitSolution();

    bool allIIEdgeComponent(Vec<ID> &subComponet, double &sRatio);

    bool MVCInitSolution(Vec<ID> &subComponet, Vec<ID> &MVCSolution);

    bool InitTopoOrder(Vec<int> &MVCSolutionINT, Vec<int> &newMVCSolutionINT, Vec<bool> &isFVS);

    //-------------------------------- online ------------------
    void loadInput(istream &is);

    void loadInputLocal(string &filename);

    // ---------------------------------------------- cache ---------------------------------------------------
    bool getFromCache;
    bool startCache; // ??????cache??????
    int CacheCheckTime = 300; // ??????????????????cache?????????
    int AddTTimesLimit = 70; // ?????????????????????????????????????????????????????????cache
    bool forbidCache; // ????????????????????????
    void runFVSP(int randSeed);

    ~Topo() {
        delete[] unNumbered;
        delete[] unNumberedPos;
        delete[] conflictNodes;
    }
};


#endif //TOPO_TOPO_H
