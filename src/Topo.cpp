#include <fstream>
#include <algorithm>
#include <iostream>
#include "Topo.h"
#include "FastVC.h"
#include "NuMVC.h"
#include "BaseCircle.h"

using namespace std;

void Topo::loadInput(istream &is) {
    int type;
    string s;
    getline(is, s);
    while (s.rfind('%', 0) == 0) {
        getline(is, s);
    }
    istringstream iss(s);
    iss >> N >> M >> type;
    connect.resize(N);
    nodeVec.resize(N);

    for (int i = 0; i < N; i++) {
        nodeVec[i] = new Node(i);
    }

    for (int i = 0; i < N; i++) {
        getline(is, s);
        while (s.rfind('%', 0) == 0) {
            getline(is, s);
        }

        istringstream iss(s);
        int n1 = i;

        int n2;
        while (iss >> n2) {
            n2 = n2 - 1;

            Node *n1Node = nodeVec[n1];
            Node *n2Node = nodeVec[n2];

            n1Node->to.push_back(n2);
            n2Node->from.push_back(n1);
            connect[n1][n2] = 1;
        }
    }
};

void Topo::eraseUnNumbered(int node) {
    int oldPos = unNumberedPos[node];
    unNumberedCount--;
    unNumbered[oldPos] = unNumbered[unNumberedCount];
    unNumberedPos[unNumbered[unNumberedCount]] = oldPos;
    unNumberedPos[node] = -1;
}

void Topo::pushUnNumbered(int node) {
    unNumbered[unNumberedCount] = node;
    unNumberedPos[node] = unNumberedCount;
    unNumberedCount++;
}

int Topo::FVSLocalSearch(double T0, double alpha, int maxMvt, int maxFail, bool runOnce) {
    int totalRunCount = 0;
    int lastTotalRunCount = 0;
    int lastBestConfCount = 0;
    double lastCostTime = 0;

    bool fail = true;
    int moveCount = 0;
    int chosen;
    int chosenIndex;
    int moveType = -1;
    int moveDelta;

    double T = T0;
    int temperatureIncreaseCount = 0;

    bool firstRun = true;

    while (true) {
        int failCount = 0;
        int lastBest = -1;

        if (!firstRun) {
            if (runOnce) {
                return N - bestConfCount;
            }
            if (!enableJumpOut) {
                //std::cerr << "enable jumpOut" << std::endl;
                enableJumpOut = true;
            }

            int delNum = (int) (cList.len * delPercent);
            double addT = increaseT;

            for (int i = 0; i < delNum; i++) {
                int delId = rand(N);
                while (cList.linkNodes[delId]->label == -1) {
                    delId = rand(N);
                }

                cList.removeNode(cList.linkNodes[delId]);
                pushUnNumbered(delId);
            }

            temperatureIncreaseCount++;
            lastBest = cList.len;
            T += addT;
            addTTimes++;

            if (startCache)
                for (int i = 0; i < N; i++) {
                    nodeVec[i]->legalPre = nodeVec[i]->legalNext = false;
                }
        }

        while (failCount < maxFail) {
            moveCount = 0;
            fail = true;
            while (moveCount < maxMvt) {
                chosenIndex = rand(unNumberedCount);
                chosen = unNumbered[chosenIndex];

                moveType = rand(2);

                bool noFrom = false;
                bool noTo = false;
                bool noRoot = false;
                int insertNextNode;
                int insertPreNode;

                int limitRightNode; // insertNextNode
                int limitLeftNode; // insertPreNode

                bool getFromCache = false;

                Node *node = nodeVec[chosen];

                if (startCache) {
                    if (moveType == 1) {
                        if (node->legalNext) {
                            getFromCache = true;
                            insertNextNode = node->insertNextNode;
                            conflictCount = node->nextConflictCount;
                            for (int i = 0; i < conflictCount; i++) {
                                conflictNodes[i] = node->nextConflictNodes[i];
                            }
                        }
                    } else if (moveType == 0) {
                        if (node->legalPre) {
                            getFromCache = true;
                            insertPreNode = node->insertPreNode;
                            conflictCount = node->preConflictCount;
                            for (int i = 0; i < conflictCount; i++) {
                                conflictNodes[i] = node->preConflictNodes[i];
                            }
                        }
                    }
                }

                if (!getFromCache) {
                    if (moveType == 1) {
                        insertNextNode = getNextPosNode(chosen);
                        if (insertNextNode == -1) {
                            noFrom = true;
                            moveType = 0;
                        }
                    }
                    if (moveType == 0) {
                        insertPreNode = getPrePosNode(chosen);
                        if (insertPreNode == -1) {
                            noTo = true;
                            if (!noFrom) {
                                insertNextNode = getNextPosNode(chosen);
                                if (insertNextNode == -1) {
                                    noFrom = true;
                                } else
                                    moveType = 1;
                            }
                            if (noFrom && noTo) {
                                moveType = 3;
                                if (cList.rootNode == nullptr) {
                                    noRoot = true;
                                    moveType = 4;
                                }
                            }
                        }
                    }

                    Node *newNode = nodeVec[chosen];
                    conflictCount = 0;
                    if (moveType == 1) {
                        lType tempL = LMAX;
                        linkNode *lkInsertNextNode = cList.linkNodes[insertNextNode];
                        for (int toNodeId: newNode->to) {
                            linkNode *lkToNode = cList.linkNodes[toNodeId];
                            if (lkToNode->label != -1 && lkToNode->label <= lkInsertNextNode->label) {
                                conflictNodes[conflictCount] = toNodeId;
                                conflictCount++;

                                if (lkToNode->label < tempL) {
                                    tempL = lkToNode->label;
                                    limitRightNode = toNodeId;
                                }
                            }
                        }
                    } else if (moveType == 0) {
                        lType tempL = -1;
                        linkNode *lkInsertPreNode = cList.linkNodes[insertPreNode];
                        for (int fromNodeId: newNode->from) {
                            linkNode *lkFromNode = cList.linkNodes[fromNodeId];
                            if (lkFromNode->label != -1 && lkFromNode->label >= lkInsertPreNode->label) {
                                conflictNodes[conflictCount] = fromNodeId;
                                conflictCount++;

                                if (lkFromNode->label > tempL) {
                                    tempL = lkFromNode->label;
                                    limitLeftNode = fromNodeId;
                                }
                            }
                        }
                    } else if (moveType == 3) {
                        conflictCount = 0;
                    } else if (moveType == 4) {
                        conflictCount = 0;
                    }

                    if (startCache) {
                        if (conflictCount >= 2) {
                            if (moveType == 1) {
                                node->legalNext = true;
                                node->insertNextNode = insertNextNode;
                                node->nextConflictCount = conflictCount;
                                for (int i = 0; i < conflictCount; i++) {
                                    node->nextConflictNodes[i] = conflictNodes[i];
                                }
                                node->rightLimitNode = limitRightNode;

                            } else if (moveType == 0) {
                                node->legalPre = true;
                                node->insertPreNode = insertPreNode;
                                node->preConflictCount = conflictCount;
                                for (int i = 0; i < conflictCount; i++) {
                                    node->preConflictNodes[i] = conflictNodes[i];
                                }
                                node->leftLimitNode = limitLeftNode;
                            }
                        }
                    }
                }

                moveDelta = conflictCount - 1;

                bool accept = moveDelta <= 0 || exp(-moveDelta / double(T)) > rand01Double();

                if (accept) {
                    if (moveType == 0) {
                        doMove(chosen, moveType, insertPreNode);
                    } else {
                        doMove(chosen, moveType, insertNextNode);
                    }

                    moveCount++;

                    if (cList.len > lastBest) {
                        lastBest = cList.len;
                        fail = false;
                    }

                    if (!checkBetterPerIter) {
                        if (totalRunCount - lastTotalRunCount > CheckInterval) {
                            if (bestConfCount < cList.len) {
                                bestConfCount = cList.len;
                                std::cerr << "updateBest " << N - cList.len + mustInFvsCount << std::endl;
                                bestUnNumbered.clear();
                                for (int i = 0; i < unNumberedCount; i++) {
                                    bestUnNumbered.push_back(unNumbered[i]);
                                }
                            }
                        }
                    } else {
                        if (bestConfCount < cList.len) {
                            bestConfCount = cList.len;
                            bestUnNumbered.clear();
                            for (int i = 0; i < unNumberedCount; i++) {
                                bestUnNumbered.push_back(unNumbered[i]);
                            }
                            std::cerr << "updateBest " << N - cList.len + mustInFvsCount << std::endl;
                        }

                    }

                    if (totalRunCount - lastTotalRunCount > CheckInterval) {
                        double endTime = clock();
                        double costTime = (double(endTime) - GlobalStartTime) / CLOCKS_PER_SEC;

                        if (!forbidCache && !startCache && costTime >= CacheCheckTime && addTTimes < AddTTimesLimit) {
                            if (costTime >= CacheCheckTime && addTTimes < AddTTimesLimit) {
                                std::cerr << "start cache!" << std::endl;
                                startCache = true;
                                forbidCache = true;
                            } else {
                                startCache = false;
                                forbidCache = true;
                            }
                        }

                        if (!checkBetterPerIter && costTime - lastCostTime >= CostTimeCheckInterVal) {
                            //std::cerr << " check " << bestConfCount <<" "<< lastBestConfCount << std::endl;
                            if (bestConfCount - lastBestConfCount <= MinBestConfCountChange) {
                                //std::cerr << "start:  check per iter!" << std::endl;
                                checkBetterPerIter = true;
                            }
                            lastCostTime = costTime;
                            lastBestConfCount = bestConfCount;
                        }

                        if (!checkBetterPerIter && costTime > EnableCheckTime) {
                            //std::cerr << "enable check per iter!" << std::endl;
                            checkBetterPerIter = true;
                        }
                        lastTotalRunCount = totalRunCount;
                    }

                    if (unNumberedCount == 0)
                        return 0;
                }

                totalRunCount++;

                if (Config::getCfg()->TLE()) {
                    break;
                }
            }

            if (fail) {
                failCount++;
                if (enableJumpOut) {
                    T = T * alpha;
                }
            } else {
                failCount = 0;
            }

            if (!enableJumpOut) {
                T = T * alpha;
            }

            if (Config::getCfg()->TLE())
                break;
        }

        firstRun = false;

        if (Config::getCfg()->TLE()) {
            break;
        }
    }

    return N - bestConfCount;
}

oneMoveRes Topo::findOneMove(int chosen, int moveType) {
    if (moveType == -1) {
        moveType = rand(2);
    }

    bool noFrom = false;
    bool noTo = false;
    bool noRoot = false;
    int insertNextNode;
    int insertPreNode;
    int baseNode = -1;

    if (moveType == 1) {
        insertNextNode = getNextPosNode(chosen);
        baseNode = insertNextNode;
        if (insertNextNode == -1) {
            noFrom = true;
            moveType = 0;
        }
    }
    if (moveType == 0) {
        insertPreNode = getPrePosNode(chosen);
        baseNode = insertPreNode;
        if (insertPreNode == -1) {
            noTo = true;
            if (!noFrom) {
                insertNextNode = getNextPosNode(chosen);
                baseNode = insertNextNode;
                if (insertNextNode == -1) {
                    noFrom = true;
                } else
                    moveType = 1;
            }
            if (noFrom && noTo) {
                moveType = 3;
                if (cList.rootNode == nullptr) {
                    noRoot = true;
                    moveType = 4;
                }
            }
        }
    }

    Node *newNode = nodeVec[chosen];
    conflictCount = 0;

    if (moveType == 1) {
        linkNode *lkInsertNextNode = cList.linkNodes[insertNextNode];
        for (int toNodeId: newNode->to) {
            linkNode *lkToNode = cList.linkNodes[toNodeId];
            if (lkToNode->label != -1 && lkToNode->label <= lkInsertNextNode->label) {
                conflictNodes[conflictCount] = toNodeId;
                conflictCount++;
            }
        }
    } else if (moveType == 0) {
        linkNode *lkInsertPreNode = cList.linkNodes[insertPreNode];
        for (int fromNodeId: newNode->from) {
            linkNode *lkFromNode = cList.linkNodes[fromNodeId];
            if (lkFromNode->label != -1 && lkFromNode->label >= lkInsertPreNode->label) {
                conflictNodes[conflictCount] = fromNodeId;
                conflictCount++;
            }
        }
    } else if (moveType == 3) {
        conflictCount = 0;
    } else if (moveType == 4) {
        conflictCount = 0;
    }

    return {moveType, conflictCount, baseNode};
}

void Topo::doMove(int chosen, int moveType, int insertBaseNode) {
    Node *node = nodeVec[chosen];

    auto &linkNodes = cList.linkNodes;
    if (moveType == 1) {
        cList.insertAfter(linkNodes[insertBaseNode], linkNodes[chosen]);
    } else if (moveType == 0) {
        cList.insertBefore(linkNodes[insertBaseNode], linkNodes[chosen]);
    } else if (moveType == 3) {
        if (rand(2) == 0) {
            cList.insertAfter(cList.rootNode, linkNodes[chosen]);
        } else {
            cList.insertBefore(cList.rootNode, linkNodes[chosen]);
        }
    } else if (moveType == 4) {
        cList.insertFirst(linkNodes[chosen]);
    }

    if (startCache) {
        for (int toId: node->to) {
            if (unNumberedPos[toId] == -1) continue;
            Node *toNode = nodeVec[toId];

            if (toNode->legalNext || toNode->legalPre) {
                auto lRange = getLabelRange(toId);
                if (lRange.first >= linkNodes[chosen]->label || lRange.second <= linkNodes[chosen]->label) {
                    toNode->legalPre = false;
                    toNode->legalNext = false;
                }
            }
        }
        for (int fromId: node->from) {
            if (unNumberedPos[fromId] == -1) continue;
            Node *fromNode = nodeVec[fromId];

            if (fromNode->legalNext || fromNode->legalPre) {
                auto lRange = getLabelRange(fromId);
                if (lRange.first >= linkNodes[chosen]->label || lRange.second <= linkNodes[chosen]->label) {
                    fromNode->legalPre = false;
                    fromNode->legalNext = false;
                }
            }
        }

        node->legalPre = false;
        node->legalNext = false;

        for (int i = 0; i < conflictCount; i++) {
            int delId = conflictNodes[i];
            Node *delNode = nodeVec[delId];
            for (int toId: delNode->to) {
                if (unNumberedPos[toId] == -1) continue;
                Node *toNode = nodeVec[toId];
                if (toNode->legalNext || toNode->legalPre) {
                    auto lRange = getLabelRange(toId);
                    if (lRange.first >= linkNodes[chosen]->label || lRange.second <= linkNodes[chosen]->label) {
                        toNode->legalPre = false;
                        toNode->legalNext = false;
                    }
                }
            }
            for (int fromId: delNode->from) {
                if (unNumberedPos[fromId] == -1) continue;
                Node *fromNode = nodeVec[fromId];
                if (fromNode->legalNext || fromNode->legalPre) {
                    auto lRange = getLabelRange(fromId);
                    if (lRange.first >= linkNodes[chosen]->label || lRange.second <= linkNodes[chosen]->label) {
                        fromNode->legalPre = false;
                        fromNode->legalNext = false;
                    }
                }
            }

            delNode->legalPre = delNode->legalNext = false;
        }
    }

    for (int i = 0; i < conflictCount; i++) {
        int delNode = conflictNodes[i];
        cList.removeNode(linkNodes[delNode]);
    }

    eraseUnNumbered(chosen);

    for (int i = 0; i < conflictCount; i++) {
        pushUnNumbered(conflictNodes[i]);
    }
}

// 缓存相关
pair<lType, lType> Topo::getLabelRange(int nodeId) {
    Node *node = nodeVec[nodeId];
    if (node->legalPre || node->legalNext) {
        auto &linkNodes = cList.linkNodes;
        if (node->legalPre) {
            lType leftL = linkNodes[node->leftLimitNode]->label;
            lType rightL = linkNodes[node->insertPreNode]->label;
            return {leftL, rightL};
        } else if (node->legalNext) {
            lType leftL = linkNodes[node->insertNextNode]->label;
            lType rightL = linkNodes[node->rightLimitNode]->label;
            return {leftL, rightL};
        }
    } else {
        return {-1, LMAX};
    }
    return {-1, LMAX};
}

int Topo::getNextPosNode(int nodeId) {
    int insertAfterNodeId = -1;
    lType tempLabel = -1;
    Node *newNode = nodeVec[nodeId];
    for (int fromNodeId: newNode->from) {
        linkNode *fromLKNode = cList.linkNodes[fromNodeId];
        if (fromLKNode->label != -1 && fromLKNode->label > tempLabel) {
            tempLabel = fromLKNode->label;
            insertAfterNodeId = fromNodeId;
        }
    }

    return insertAfterNodeId;
}

int Topo::getPrePosNode(int nodeId) {
    int insertPreNode = -1;
    lType tempLabel = LMAX;
    Node *newNode = nodeVec[nodeId];
    for (int toNodeId: newNode->to) {

        linkNode *toLKNode = cList.linkNodes[toNodeId];
        if (toLKNode->label != -1 && toLKNode->label < tempLabel) {
            tempLabel = toLKNode->label;
            insertPreNode = toNodeId;
        }
    }

    return insertPreNode;
}


void Topo::init() {
    for (int i = 0; i < N; i++) {
        Node *node = nodeVec[i];
        node->init(node->from.size(), node->to.size());
    }

    cList.init(N);

    unNumbered = new int[N];
    unNumberedPos = new int[N];
    unNumberedCount = N;

    conflictNodes = new int[N];
    bestConfCount = 0;

    for (int i = 0; i < N; i++) {
        unNumbered[i] = i;
        unNumberedPos[i] = i;
    }

    enableJumpOut = false;
    checkBetterPerIter = false;

    // 缓存
    //forbidCache = true;
    forbidCache = false;
    startCache = false;

    bestUnNumbered.clear();
    for (int i = 0; i < unNumberedCount; i++) {
        bestUnNumbered.push_back(unNumbered[i]);
    }
}


bool Topo::selfLoop() {
    bool hasSelfLoop = false;
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;
        if (connect[i][i]) {
            hasSelfLoop = true;
            reduceVertex(i, true);
            mustInFVS[i] = true;
        }
    }
    return hasSelfLoop;
}


void Topo::getFVSolution() {
    fvSolution.clear();
    for (int i = 0; i < originN; i++) {
        if (mustInFVS[i]) {
            fvSolution.push_back(i);
        }
    }

    //cerr<<" bestUnnumeredSize:" << bestUnNumbered.size() << endl;
    for (int nodeId: bestUnNumbered) {
        int originId = reducedToOriginMap[nodeId];
        fvSolution.push_back(originId);
    }
}


void handle(int sigNum) {
    Config::getCfg()->setFlag();
}


/********************************minimum vertex cover ********************************************/
void Topo::InitSolution() {
    Vec<ID> MVCSolution;
    Vec<ID> subComponet;
    for (int i = 0; i != N; ++i) {
        subComponet.push_back(i);
    }
    if (!MVCInitSolution(subComponet, MVCSolution))
        return;

    ofstream result;

    if (MVCSolution.empty()) {
        return;
    } else {
        Vec<bool> isFVS(N, 0);
        for (int nodeId: MVCSolution) {
            isFVS[nodeId] = 1;
        }
        vector<int> TopoNodes;
        for (int i = 0; i < N; i++) {
            if (isFVS[i] != 1) {
                TopoNodes.push_back(i);
            }
        }
        //初始拓扑序列
        Vec<int> newTopoNodes;
        if (!InitTopoOrder(TopoNodes, newTopoNodes, isFVS)) {
            return;
        }

        unNumberedCount = 0;
        for (int i = 0; i < N; i++) {
            unNumbered[i] = -1;
            unNumberedPos[i] = -1;
        }
        for (int nodeId: MVCSolution) {
            pushUnNumbered(nodeId);
        }
        //拓扑序列更新      
        cList.setInitSol(newTopoNodes);

        // 更新best
        bestConfCount = cList.len;
        bestUnNumbered.clear();
        for (int i = 0; i < unNumberedCount; i++) {
            bestUnNumbered.push_back(unNumbered[i]);
        }
    }
}

bool Topo::allIIEdgeComponent(Vec<ID> &subComponet, double &sRatio) {
    ID singleEdgeNum = 0;
    ID totalEdgeNum = 0;

    for (auto src = 0; src != N; ++src) {
        Node *Neighbour = nodeVec[src];
        totalEdgeNum += Neighbour->to.size();
        for (auto dst = Neighbour->to.begin(); dst != Neighbour->to.end(); ++dst) {
            if (Find(nodeVec[*dst]->to, src) == nodeVec[*dst]->to.end())singleEdgeNum++;
        }
    }
    sRatio = singleEdgeNum * 1.0 / totalEdgeNum;
    if (sRatio > singleRatio)return false;
    return true;
}

bool Topo::InitTopoOrder(Vec<int> &MVCSolutionINT, Vec<int> &newMVCSolutionINT, Vec<bool> &isFVS) {
    Vec<int> in_degree(N, 0);
    queue<int> p;
    int ConNum = 0;
    for (auto src = MVCSolutionINT.begin(); src != MVCSolutionINT.end(); ++src) {
        Node *newNode = nodeVec[*src];
        for (auto dst = newNode->to.begin(); dst != newNode->to.end(); ++dst) {
            if (isFVS[*dst] == 1)continue;
            in_degree[*dst]++;
        }
    }
    for (auto src = MVCSolutionINT.begin(); src != MVCSolutionINT.end(); ++src) {
        if (isFVS[*src] != 1 && in_degree[*src] == 0)
            p.push(*src);
    }
    while (!p.empty()) {
        ID u = p.front();
        p.pop();
        newMVCSolutionINT.push_back(u);
        ConNum++;
        Node *newNode = nodeVec[u];
        for (const int dst: newNode->to) {
            if (isFVS[dst] == 1)continue;
            --in_degree[dst];
            if (in_degree[dst] == 0) {
                p.push(dst);
            }
        }
    }
    if (ConNum != MVCSolutionINT.size())return false;
    return true;
}

bool Topo::MVCInitSolution(Vec<ID> &subComponet, Vec<ID> &MVCSolution) {
    Vec<ID> subGraph;
    HashMap<ID, ID> subGraph_map;
    Vec<int> mvc_solution;
    for (auto src = subComponet.begin(); src != subComponet.end(); ++src) {
        subGraph.push_back(*src);
        subGraph_map[*src] = subGraph.size() - 1;
    }
    double sRatio = 0.0;
    bool allTTEdgeComponent = allIIEdgeComponent(subComponet, sRatio);
    if (!allTTEdgeComponent) {
        for (auto src = subComponet.begin(); src != subComponet.end(); ++src) {
            MVCSolution.push_back(*src);
        }
        return false;
    }

    int InitRunTime = 0;
    if (sRatio == 0.0) {
        double currTime = (double(clock()) - GlobalStartTime) / CLOCKS_PER_SEC;
        InitRunTime = (EnableCheckTime - currTime) * 0.3 * 1000;//使用1/3的时间产生初始解ms*1000
        InitRunTime = min(InitRunTime, 200000);//min(,200s)
    } else if (sRatio < 0.1) {
        double currTime = (double(clock()) - GlobalStartTime) / CLOCKS_PER_SEC;
        InitRunTime = (EnableCheckTime - currTime) * 0.1 * 1000;//使用1/10的时间产生初始解
        InitRunTime = min(InitRunTime, 50000);//min(,50s)
    } else if (sRatio < 0.35) {
        if (subComponet.size() < initSolutionSelect) {
            InitRunTime = 20000;//20s
        } else {
            InitRunTime = 30000;
        }
    } else {
        if (subComponet.size() < initSolutionSelect) {
            InitRunTime = 15000;//10s
        } else {
            InitRunTime = 20000;
        }
    }
    if (InitRunTime < 0) {
        InitRunTime = 10000;
    }

    bool flag1, flag2;
    if (subComponet.size() < initSolutionSelect) {

        Vec<NuMVC::Edge> MVCEdge;
        for (auto src = 0; src != subComponet.size(); ++src) {
            Node *Neighbour = nodeVec[subGraph[src]];
            for (auto dst = Neighbour->to.begin(); dst != Neighbour->to.end(); ++dst) {
                if (src >= subGraph_map[*dst])continue;
                MVCEdge.push_back({src, subGraph_map[*dst]});
            }
        }
        NuMVC m(0, std::chrono::milliseconds(InitRunTime));
        flag2 = m.NuMVC_minVertexCovSolver(MVCEdge, subComponet.size(), mvc_solution);
    } else {
        Vec<FastVC::Edge> MVCEdge;
        for (auto src = 0; src != subComponet.size(); ++src) {
            Node *Neighbour = nodeVec[subGraph[src]];
            for (auto dst = Neighbour->to.begin(); dst != Neighbour->to.end(); ++dst) {
                if (src >= subGraph_map[*dst])continue;
                MVCEdge.push_back({src, subGraph_map[*dst]});
            }
        }
        FastVC f(0, std::chrono::milliseconds(InitRunTime));//30s
        flag2 = f.FastVC_minVertexCovSolver(MVCEdge, subComponet.size(), mvc_solution);
    }
    if (!flag2)return false;
    for (auto i = mvc_solution.begin(); i != mvc_solution.end(); ++i) {
        MVCSolution.push_back(subGraph[*i]);
    }
    return true;
}

void Topo::runFVSP(int randSeed) {
    GlobalStartTime = clock();

    EnableCheckTime = 400;

    initRand(randSeed);

    loadInput(cin);

    reduction();

    init();

    if (N != 0) {
        int mvCount = 10000;
        double avgEdge = double(edgeCount) / N;

        if (avgEdge >= 50) {
            mvCount = min(8000, N / 5);
        }

        if (N <= 3000 && edgeCount <= 20000) {
            FVSLocalSearch(1.5, 0.99, mvCount, 50, true); // 只运行一次
            BaseCircle bc;
            vector<vector<int>> adj(N); // 图的邻接关系
            for (int i = 0; i < N; i++) {
                adj[i] = nodeVec[i]->to;
            }
            bc.init(N, adj);
            // 随机数、上限时间、初始的反馈集
            bestUnNumbered = bc.runFVSP(rand(100000), 55, bestUnNumbered);

            bestConfCount = N - bestUnNumbered.size();

            double endTime = clock();
            double costTime = (double(endTime - GlobalStartTime)) / CLOCKS_PER_SEC;
        } else {
            InitSolution();
        }

        FVSLocalSearch(1.5, 0.995, mvCount, 50);
    }

    // 获得解
    getFVSolution();
}


Config *Config::cfg = new Config();

int main(int argc, char *argv[]) {
    // 控制信号
    signal(SIGTERM, handle);

    int randSeed = 7654321;

    Topo topo;
    topo.runFVSP(randSeed);

    double endTime = clock();
    double costTime = (double(endTime - topo.GlobalStartTime)) / CLOCKS_PER_SEC;

    int fvsSize = topo.fvSolution.size();
    //cerr << "FVS size:" << fvsSize << endl;

    for (int nodeId: topo.fvSolution) {
        std::cout << nodeId + 1 << std::endl;
    }

}