#include <random>
#include <chrono>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <queue>

using NodeId = int;
using EdgeId = NodeId;
using Nodes = std::vector<NodeId>;

using Centers = Nodes;

class PCenter {
public:
    // 初始时间
    double startTime = clock();
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

    double endTime;
    double breakTime;
    std::vector<std::vector<int>> BasicCycles;

    void initTimer(long long secTimeout) { breakTime = startTime + secTimeout * CLOCKS_PER_SEC; }

    bool isTimeout() {
        double nowTime = clock();
        return breakTime < (nowTime - 0.2 * CLOCKS_PER_SEC); // 提前0.2s
    }

    int nodeNum;
    int cycleNum;
    int addedCycleCount;
    int centerNum;

    std::vector<Nodes> coverages;
    std::vector<Nodes> canBeCoverd;

    std::vector<int> weights;
    std::vector<int> deltas;
    Centers nowCenters;
    Centers bestCenters;
    std::vector<int> tabuList;
    std::vector<int> servedCentersCount;
    std::vector<NodeId> servedBy;
    std::vector<bool> solution;

    int currentIter;
    int centerIndex;

    std::vector<int> unservedNodes;

    int lastUnservedCounts;
    int minUnservedCounts;

    std::vector<std::vector<int>> candidates;
    int candidateCount;
    std::vector<int> openCandidates;
    std::vector<int> closeCandidates;
    int closeCandidateCount = -1;
    int openCandidateCount = -1;

    std::vector<bool> alreadyBeServed;

    std::mt19937 pseudoRandNumGen;

    int rand(int lb, int ub) { return std::uniform_int_distribution<int>(lb, ub - 1)(pseudoRandNumGen); }

    int rand(int ub) { return std::uniform_int_distribution<int>(0, ub - 1)(pseudoRandNumGen); }

    const int P1 = 5; // 0~10
    const int P2 = 10; // 1~30
    const int P3 = 4;

    std::string fileName;
    std::string outCsvPath;
    std::queue<std::string> resStringQueue;

public:
    int sumIter = 0;

    void upperBoundGreedy(std::vector<int> &upperBoundNOdes);

    void initDelta();

    void openNode(NodeId nodei);

    void closeNode(NodeId nodej);

    std::pair<int, int> findMove();

    std::pair<int, int> findGreedyMove();

    void tryToOpen(NodeId nodei);

    void restoreDeltas(int nodei);

    bool listInclude(Nodes nodes1, Nodes nodes2);

    void openNodeAndDecreaseCount(NodeId openNode);

    void move(int nodei, int nodej);

    bool checkAndUpdate();

    void initRand(int seed) { pseudoRandNumGen = std::mt19937(seed); }

    bool check() {
        std::vector<bool> temp;
        temp.resize(cycleNum);
        for (int i = 0; i < cycleNum; ++i) {
            temp[i] = false;
        }
        for (int i = 0; i < centerNum; ++i) {
            //int c = bestCenters[i];
            int c = nowCenters[i];
            int sc = coverages[c].size();
            for (int j = 0; j < sc; ++j) {
                if (servedCentersCount[coverages[c][j]] < 1) {
                    std::cerr << " err" << servedCentersCount[coverages[c][j]] << std::endl;
                }
            }

            for (int j = 0; j < sc; ++j) {
                temp[coverages[c][j]] = true;
            }
        }

        for (int i = 0; i < cycleNum; ++i) {
            if (!temp[i])return false;
        }
        return true;
    }

    // 初始化
    void initFVSP(int typeASize, std::vector<std::vector<int>> &basicCycle, std::vector<int> upperBoundNodes, int seed);

    // 开始迭代
    void runIteration();

    void downSizing();

    // 发现还有循环，更新cover和canBe
    void updateCover(std::vector<std::vector<int>> cycles);
};
