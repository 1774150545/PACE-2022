//
#include "Topo.h"
#include "LMCMaxClique.h"

bool Topo::doGraphReduction2() {
    bool TryAllRelabel = true;

    if (N < 100000 && N != 0 && int(M / N) > 90) {
        TryAllRelabel = false;
    }

    int sumCount = 0;
    int tempCount = 0;

    while (true) {
        int reduceCountPre5 = 0;
        bool hasSelfLoop = false;

        do {
            hasSelfLoop = false;
            reduceCountPre5 = 0;
            // 1.IN0
            for (int i = 0; i < N; i++) {
                if (isReduced[i]) continue;
                if (nodeVec[i]->from.empty()) {
                    reduceVertex(i);
                    tempCount++;
                    reduceCountPre5++;
                }
            }

            // 2.Out0
            tempCount = 0;
            for (int i = 0; i < N; i++) {
                if (isReduced[i]) continue;
                if (nodeVec[i]->to.empty()) {
                    tempCount++;
                    reduceCountPre5++;
                    sumCount++;
                    reduceVertex(i);
                }
            }

            // 3.selfLoop
            tempCount = 0;
            hasSelfLoop = clearSelfLoop();


            tempCount = 0;
            // 4.In1
            for (int i = 0; i < N; i++) {
                int delId = i;
                if (isReduced[i]) continue;
                Node *delNode = nodeVec[delId];
                if (delNode->from.size() == 1) {
                    if (connect[delId][delId]) {
                        continue;
                    }

                    tempCount++;
                    sumCount++;
                    reduceCountPre5++;

                    int shouldFromId = delNode->from[0];
                    Node *shouldFromNode = nodeVec[shouldFromId];

                    delNode->from.erase(delNode->from.begin());
                    int tempPos = shouldFromNode->inTo(delId);

                    shouldFromNode->to.erase(shouldFromNode->to.begin() + tempPos);
                    connect[shouldFromId][delId] = 0;

                    for (int toId: delNode->to) {
                        Node *toNode = nodeVec[toId];

                        int theDelPos = toNode->inFrom(delId);

                        int isConnected = connect[shouldFromId][toId];

                        if (isConnected) {
                            int theShouldFromPos = toNode->inFrom(shouldFromId);

                            toNode->from.erase(toNode->from.begin() + theDelPos);
                            connect[delId][toId] = 0;
                        } else {
                            toNode->from[theDelPos] = shouldFromId;

                            shouldFromNode->to.push_back(toId);

                            connect[shouldFromId][toId] = 1;
                            connect[delId][toId] = 0; // new
                        }
                    }

                    isReduced[delId] = true;

                    delNode->to.clear();
                }
            }

            // 5.out1
            bool temp = clearSelfLoop();
            hasSelfLoop = temp || hasSelfLoop;
            tempCount = 0;
            for (int i = 0; i < N; i++) {
                int delId = i;
                if (isReduced[i]) continue;
                Node *delNode = nodeVec[delId];
                if (delNode->to.size() == 1) {
                    if (connect[delId][delId]) {
                        continue;
                    }
                    tempCount++;
                    sumCount++;
                    reduceCountPre5++;

                    int shouldToId = delNode->to[0];
                    Node *shouldToNode = nodeVec[shouldToId];

                    delNode->to.erase(delNode->to.begin());
                    int tempPos = shouldToNode->inFrom(delId);
                    shouldToNode->from.erase(shouldToNode->from.begin() + tempPos);
                    connect[delId][shouldToId] = 0;

                    for (int fromId: delNode->from) {
                        Node *fromNode = nodeVec[fromId];
                        int theDelPos = fromNode->inTo(delId);

                        int isConnected = connect[fromId][shouldToId];

                        if (isConnected) {
                            int theShouldToPos = fromNode->inTo(shouldToId);
                            fromNode->to.erase(fromNode->to.begin() + theDelPos);
                            connect[fromId][delId] = 0;
                        } else {
                            fromNode->to[theDelPos] = shouldToId;
                            shouldToNode->from.push_back(fromId);
                            connect[fromId][shouldToId] = 1;
                        }
                    }

                    isReduced[delId] = true;
                    delNode->from.clear();
                    delNode->to.clear();
                }
            }
        } while (reduceCountPre5 != 0 || hasSelfLoop);

        bool haveReduced = false;
        if (TryAllRelabel) {
            hasSelfLoop = false;

            // 进行 6、7、8
            hasSelfLoop = clearSelfLoop();

            // 6.PIE
            haveReduced = pie();

            bool temp = clearSelfLoop();
            hasSelfLoop = temp || hasSelfLoop;

            // 7.CORE
            int reduceCountCore = 0;
            vector<int> mustReduce = core();
            for (int nodeId: mustReduce) {
                mustInFVS[nodeId] = true;
            }
            haveReduced = haveReduced || (!mustReduce.empty());
            temp = clearSelfLoop();
            hasSelfLoop = temp || hasSelfLoop;

            // 8.dome
            haveReduced = haveReduced || dome();
        }


        if (haveReduced || hasSelfLoop) continue;
        else {
            //std::cerr << "run maxClique!" << std::endl;
            //maxCliqueReduce();
            //std::cerr <<" run maxClique finish!" <<endl;
            break;
        }

    }

    return false;
}

/*
* Strongly connected operate
*/
void Topo::maxCliqueReduce() {
    ID ComponentNumber = 0;
    Vec<ID> ComponentSubID(N, 0);
    Vec<Vec<ID>> Component;
    vertexToStronglyConnectedComponentNumberNR();
    Component.resize(currentCpmponent);
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;
        int componentId = sccsByNum[i];
        Component[componentId].push_back(i);
    }
    int deleteNum = 0;
    int num = 0;
    for (auto i: Component) {
        viaClique(i);
        ++num;
    }
}


bool Topo::allIIEdgeJudge(Vec<ID> &ComponentSub) {
    for (auto src = ComponentSub.begin(); src != ComponentSub.end(); ++src) {
        if (isReduced[*src]) continue;
        Node *theNode = nodeVec[*src];
        for (auto dst = theNode->to.begin(); dst != theNode->to.end(); ++dst) {
            if (isReduced[*dst])continue;
            if (!connect[*dst][*src]) {
                return false;
            }
        }

    }
    return true;
}

void Topo::viaClique(Vec<ID> &ComponentSub) {
    Vec<ID> subGraph;
    HashMap<ID, ID> subGraph_map;
    ID newNodeNum = ComponentSub.size();
    ID newEdgeNum = 0;
    if (!allIIEdgeJudge(ComponentSub))return;
    Vec<LMCMaxClique::NewEdge> newEdges;
    Vec<bool> Solver(ComponentSub.size() + 1, 0);
    if (newNodeNum > viaCliqueMaxNodeNum)return;
    for (auto j = ComponentSub.begin(); j != ComponentSub.end(); ++j) {
        subGraph.push_back(*j);
        subGraph_map[*j] = subGraph.size() - 1;
    }
    Vec<Vec<bool>> sign(ComponentSub.size(), Vec<bool>(ComponentSub.size(), 0));
    for (auto j = 0; j != ComponentSub.size(); ++j) {
        if (isReduced[subGraph[j]]) { continue; }
        Vec<ID> visited(ComponentSub.size(), 0);
        Node *theNode = nodeVec[subGraph[j]];
        for (auto dst = theNode->to.begin(); dst != theNode->to.end(); ++dst) {
            if (isReduced[*dst])continue;
            visited[subGraph_map[*dst]] = 1;
        }
        for (auto dst = 0; dst != ComponentSub.size(); ++dst) {
            if (isReduced[subGraph[dst]])continue;
            if (j >= dst)continue;
            if (visited[dst] == 1)continue;
            newEdges.push_back({j + 1, dst + 1});
            newEdgeNum++;
        }
    }
    LMCMaxClique m;
    Vec<int> mvc_solution;
    bool flag = false;
    int tt = 0;
    flag = m.maxCliqueSolver(newEdges, newNodeNum, Solver, process_max_time);
    if (flag == 1) {
        for (auto j = 1; j != ComponentSub.size() + 1; ++j) {
            if ((Solver[j]) == 1) {
                ID t = j - 1;
                isReduced[subGraph[t]] = true;

            } else {
                ID t = j - 1;
                isReduced[subGraph[t]] = true;
                mustInFVS[subGraph[t]] = true;
                ++tt;
            }
        }
    }
}

void Topo::reduction() {
    int oldN = N;

    isReduced.resize(N, false);
    mustInFVS.resize(N, false);

    while (true) {
        bool hasReduced = doGraphReduction2();
        if (!hasReduced) {
            int afterRedN = 0;
            for (int i = 0; i < isReduced.size(); i++) {
                if (!isReduced[i]) afterRedN++;
            }
            edgeCount = 0;
            for (int i = 0; i < N; i++) {
                if (!isReduced[i]) {
                    edgeCount += nodeVec[i]->from.size();
                    edgeCount += nodeVec[i]->to.size();
                }
            }
            edgeCount /= 2;
            break;
        }
    }

    int afterN = 0;
    for (int i = 0; i < N; i++) {
        if (!isReduced[i]) {
            originToReducedMap[i] = afterN;
            reducedToOriginMap[afterN] = i;
            afterN++;
        }
    }

    nodeVecTemp.resize(afterN);
    for (int i = 0; i < afterN; i++) {
        nodeVecTemp[i] = new Node(i);
    }

    for (int i = 0; i < N; i++) {
        if (!isReduced[i]) {
            // 旧节点
            Node *theNode = nodeVec[i];
            int newId = originToReducedMap[i];
            // 新节点
            Node *newNode = nodeVecTemp[newId];

            for (int fromId: theNode->from) {
                int newFrom = originToReducedMap[fromId];
                newNode->from.push_back(newFrom);
            }
            for (int toId: theNode->to) {
                int newTo = originToReducedMap[toId];
                newNode->to.push_back(newTo);
            }
        }
    }

    // 清空旧的
    for (int i = 0; i < N; i++) {
        nodeVec[i]->from.clear();
        nodeVec[i]->to.clear();
        delete nodeVec[i];
    }

    for (int i = 0; i < afterN; i++) {
        nodeVec[i] = nodeVecTemp[i]; // 指针
    }

    mustInFvsCount = 0;
    for (int i = 0; i < N; i++) {
        if (mustInFVS[i])
            mustInFvsCount++;
    }

    // 更新N
    originN = N;
    N = afterN;
    return;
}

void Topo::reduceVertex(int delId, bool selfLoop) {
    isReduced[delId] = true;
    Node *delNode = nodeVec[delId];
    // 额外处理一下selfLoop
    if (selfLoop) {
        reduceEdge({delId, delId});
    }

    for (int fromId: delNode->from) {
        Node *fromNode = nodeVec[fromId];
        int tempPos = fromNode->inTo(delId);
        fromNode->to.erase(fromNode->to.begin() + tempPos);
        connect[fromId][delId] = 0;
    }

    for (int toId: delNode->to) {
        Node *toNode = nodeVec[toId];
        int tempPos = toNode->inFrom(delId);
        toNode->from.erase(toNode->from.begin() + tempPos);
        connect[delId][toId] = 0;
    }

    // 也可以不清空
    delNode->from.clear();
    delNode->to.clear();
}

bool Topo::pie() {
    vector<Edge> pies = piEdges();

    for (auto edge: pies) {
        int nodeId1 = edge.first;
        int nodeId2 = edge.second;
        Node *node1 = nodeVec[nodeId1];
        Node *node2 = nodeVec[nodeId2];

        // 删除to
        int tempPos1 = node1->inTo(nodeId2);
        node1->to.erase(node1->to.begin() + tempPos1);
        int tempPos2 = node2->inTo(nodeId1);
        node2->to.erase(node2->to.begin() + tempPos2);
    }

    // 子图构造完毕，寻找强连通分量，并寻找不在循环中的边
    auto delEdges = acyclicEdges();

    // 真实删除这些边
    for (auto edge: delEdges) {
        int nodeId1 = edge.first;
        int nodeId2 = edge.second;
        // 简洁化
        reduceEdge({nodeId1, nodeId2});
    }

    // 恢复双向边
    for (auto edge: pies) {
        int nodeId1 = edge.first;
        int nodeId2 = edge.second;
        Node *node1 = nodeVec[nodeId1];
        Node *node2 = nodeVec[nodeId2];

        // 恢复to
        node1->to.push_back(nodeId2);
        node2->to.push_back(nodeId1);
    }

    if (delEdges.size() != 0)
        return true;

    return false;
}

// 不参与循环的边
vector<Edge> Topo::acyclicEdges() {
    vector<Edge> acyclicEdges;

    auto nodeToSCC = vertexToStronglyConnectedComponentNumberNR();

    // 遍历所有的边
    for (int nodeId = 0; nodeId < N; nodeId++) {
        if (isReduced[nodeId]) continue;
        Node *theNode = nodeVec[nodeId];
        for (int toId: theNode->to) {
            if (nodeToSCC[nodeId] != nodeToSCC[toId]) {
                acyclicEdges.emplace_back(nodeId, toId); // 添加进去这条边
            }
        }
    }

    return acyclicEdges;
}

vector<Edge> Topo::piEdges() {
    // pi边
    vector<Edge> piEdges;

    for (int nodeId = 0; nodeId < N; nodeId++) {
        if (isReduced[nodeId]) continue;
        Node *theNode = nodeVec[nodeId];
        for (int toId: theNode->to) {
            if (toId <= nodeId) continue;

            if (connect[toId][nodeId]) {
                piEdges.emplace_back(toId, nodeId);
            }
        }
    }

    return piEdges;
}


vector<int> Topo::vertexToStronglyConnectedComponentNumberNR() {
    tarjanIndex.clear();
    tarjanAncestor.clear();
    tarjanAncestor.clear();
    sccsByNum.clear();

    tarjanCurIndex = 0;
    tarjanIndex = vector<int>(N);
    tarjanAncestor = vector<int>(N);
    currentCpmponent = 0;

    sccsByNum = vector<int>(N);
    for (int i = 0; i < N; i++) {
        tarjanIndex[i] = -1;
    }

    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;
        if (tarjanIndex[i] == -1) {
            dfsNR(i);
        }
    }

    return sccsByNum;
}

void Topo::dfsNR(int nodeId) {
    st.push(nodeId);
    unordered_map<int, int> tempParent;
    tempParent[nodeId] = nodeId;
    while (!st.empty()) {
        int v = st.top(); // peek
        if (tarjanIndex[v] == -1) { // 新的节点
            trace.push(v);
            tarjanIndex[v] = tarjanAncestor[v] = ++tarjanCurIndex;
        }
        bool flag = false;
        for (int toId: nodeVec[v]->to) {
            if (isReduced[toId]) { continue; }
            tempParent[toId] = v;
            if (tarjanIndex[toId] == -1) {
                flag = true;
                st.push(toId);
            } else {  // 已经被访问
                tarjanAncestor[v] = min(tarjanAncestor[v], tarjanAncestor[toId]);
            }
        }

        if (!flag) {
            v = st.top();
            st.pop();
            tarjanAncestor[tempParent[v]] = min(tarjanAncestor[tempParent[v]], tarjanAncestor[v]);
            if (tarjanAncestor[v] == tarjanIndex[v]) {
                // 是一个连通分量
                int cur = -1;
                do {
                    cur = trace.top();
                    trace.pop();
                    tarjanAncestor[cur] = INT32_MAX;
                    sccsByNum[cur] = currentCpmponent;
                } while (cur != v);
                currentCpmponent++;
            }
        }

    }
}


vector<int> Topo::core() {
    vector<int> piVertexes; // 有pi边的点的集合
    vector<int> mustReduce;

    // 获取 piVertex
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;

        bool isPiVertex = true;
        for (int toId: nodeVec[i]->to) {
            if (!connect[toId][i]) { // || (toId < i)  指定顺序避免重复
                isPiVertex = false;
                break;
            }
        }

        if (isPiVertex) {
            piVertexes.push_back(i);
        }
    }

    // 最小的团是一个孤立点，也会把它删除

    for (int piV: piVertexes) {
        vector<int> potentialClique = nodeVec[piV]->to;
        potentialClique.push_back(piV);

        if (isClique(potentialClique)) {
            for (int nodeId: potentialClique) {
                reduceVertex(nodeId);
            }
            potentialClique.pop_back();
            mustReduce.insert(mustReduce.end(), potentialClique.begin(), potentialClique.end());
        }
    }

    return mustReduce;
}

bool Topo::isClique(vector<int> &nodes) {
    bool isClique = true;

    for (int node1: nodes) {
        for (int node2: nodes) {
            if (node1 == node2) continue;
            if (!connect[node1][node2] || !connect[node1][node2]) {
                isClique = false;
                break;
            }
        }
    }

    return isClique;
}

bool Topo::dome() {
    bool hasReduce = false;
    int reducedEdgeCount = 0;

    vector<Edge> toReduce;

    // 遍历所有的边
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;
        for (int toId: nodeVec[i]->to) {

            Edge edge = {i, toId};
            if (isDominate(edge)) { // 是支配边，则删除
                toReduce.push_back(edge);
                hasReduce = true;
                reducedEdgeCount++;
            }
        }
    }

    for (auto edge: toReduce) {
        reduceEdge(edge);
    }
    return hasReduce;
}


bool Topo::isDominate(Edge e) {
    int source = e.first;
    int target = e.second;

    Node *sourceNode = nodeVec[source];
    Node *targetNode = nodeVec[target];

    if (connect[target][source])
        return false;

    bool case2 = true;
    for (int nonPiTo: nodeVec[target]->to) {
        if (connect[nonPiTo][target]) continue;
        if (!connect[source][nonPiTo]) {
            case2 = false;
            break;
        }
    }
    if (case2) return true;

    bool case1 = true;
    for (int nonPiFrom: nodeVec[source]->from) {
        if (connect[source][nonPiFrom]) continue; // 是pi边
        if (!connect[nonPiFrom][target]) {
            case1 = false;
            break;
        }
    }
    if (case1) return true;

    return false;
}

void Topo::reduceEdge(Edge e) {
    int source = e.first;
    int target = e.second;

    Node *sourceNode = nodeVec[source];
    Node *targetNode = nodeVec[target];

    // 删除这条边
    int tempPos1 = sourceNode->inTo(target);
    sourceNode->to.erase(sourceNode->to.begin() + tempPos1);
    int tempPos2 = targetNode->inFrom(source);
    targetNode->from.erase(targetNode->from.begin() + tempPos2);

    connect[source][target] = 0;
}


bool Topo::clearSelfLoop() {
    bool hasSelfLoop = selfLoop();
    while (selfLoop());

    return hasSelfLoop;
}