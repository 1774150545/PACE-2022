#include "Topo.h"

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

            // ?????? 6???7???8
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
            // ?????????
            Node *theNode = nodeVec[i];
            int newId = originToReducedMap[i];
            // ?????????
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

    // ????????????
    for (int i = 0; i < N; i++) {
        nodeVec[i]->from.clear();
        nodeVec[i]->to.clear();
        delete nodeVec[i];
    }

    for (int i = 0; i < afterN; i++) {
        nodeVec[i] = nodeVecTemp[i]; // ??????
    }

    mustInFvsCount = 0;
    for (int i = 0; i < N; i++) {
        if (mustInFVS[i])
            mustInFvsCount++;
    }

    // ??????N
    originN = N;
    N = afterN;
    return;
}

void Topo::reduceVertex(int delId, bool selfLoop) {
    isReduced[delId] = true;
    Node *delNode = nodeVec[delId];
    // ??????????????????selfLoop
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

    // ??????????????????
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

        // ??????to
        int tempPos1 = node1->inTo(nodeId2);
        node1->to.erase(node1->to.begin() + tempPos1);
        int tempPos2 = node2->inTo(nodeId1);
        node2->to.erase(node2->to.begin() + tempPos2);
    }

    // ???????????????????????????????????????????????????????????????????????????
    auto delEdges = acyclicEdges();

    // ?????????????????????
    for (auto edge: delEdges) {
        int nodeId1 = edge.first;
        int nodeId2 = edge.second;
        // ?????????
        reduceEdge({nodeId1, nodeId2});
    }

    // ???????????????
    for (auto edge: pies) {
        int nodeId1 = edge.first;
        int nodeId2 = edge.second;
        Node *node1 = nodeVec[nodeId1];
        Node *node2 = nodeVec[nodeId2];

        // ??????to
        node1->to.push_back(nodeId2);
        node2->to.push_back(nodeId1);
    }

    if (delEdges.size() != 0)
        return true;

    return false;
}

// ?????????????????????
vector<Edge> Topo::acyclicEdges() {
    vector<Edge> acyclicEdges;

    auto nodeToSCC = vertexToStronglyConnectedComponentNumberNR();

    // ??????????????????
    for (int nodeId = 0; nodeId < N; nodeId++) {
        if (isReduced[nodeId]) continue;
        Node *theNode = nodeVec[nodeId];
        for (int toId: theNode->to) {
            if (nodeToSCC[nodeId] != nodeToSCC[toId]) {
                acyclicEdges.emplace_back(nodeId, toId); // ?????????????????????
            }
        }
    }

    return acyclicEdges;
}

vector<Edge> Topo::piEdges() {
    // pi???
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
        if (tarjanIndex[v] == -1) { // ????????????
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
            } else {  // ???????????????
                tarjanAncestor[v] = min(tarjanAncestor[v], tarjanAncestor[toId]);
            }
        }

        if (!flag) {
            v = st.top();
            st.pop();
            tarjanAncestor[tempParent[v]] = min(tarjanAncestor[tempParent[v]], tarjanAncestor[v]);
            if (tarjanAncestor[v] == tarjanIndex[v]) {
                // ?????????????????????
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
    vector<int> piVertexes; // ???pi??????????????????
    vector<int> mustReduce;

    // ?????? piVertex
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;

        bool isPiVertex = true;
        for (int toId: nodeVec[i]->to) {
            if (!connect[toId][i]) { // || (toId < i)  ????????????????????????
                isPiVertex = false;
                break;
            }
        }

        if (isPiVertex) {
            piVertexes.push_back(i);
        }
    }

    // ???????????????????????????????????????????????????

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

    // ??????????????????
    for (int i = 0; i < N; i++) {
        if (isReduced[i]) continue;
        for (int toId: nodeVec[i]->to) {

            Edge edge = {i, toId};
            if (isDominate(edge)) { // ???????????????????????????
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
        if (connect[source][nonPiFrom]) continue; // ???pi???
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

    // ???????????????
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