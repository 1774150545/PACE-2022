#ifndef TOPO_MYLIST_H

#include <vector>
#include <iostream>
#include <unordered_map>

typedef long long lType;
#define LMAX INT64_MAX

using namespace std;

struct linkNode {
    int nodeId;
    lType label;

    linkNode *preNode;
    linkNode *nextNode;

    linkNode() {
        this->label = -1;
        preNode = nullptr;
        nextNode = nullptr;
    };

    linkNode(int nodeId) {
        this->nodeId = nodeId;
        this->label = -1;
        preNode = nullptr;
        nextNode = nullptr;
    }

    ~linkNode() {
        preNode = nullptr;
        nextNode = nullptr;
    }
};


class myList {
public:
    int len;

    vector<linkNode *> linkNodes;

    myList() = default;

    void init(int N) {
        for (int i = 0; i < N; i++) {
            linkNodes.push_back(new linkNode(i));
        }
        len = 0;
    }

    myList(int N) {
        init(N);
    }

    void setInitSol(vector<int> &sol);

    void insertFirst(linkNode *baseNode);

    void insertBase(linkNode *baseNode);


    bool insertAfter(linkNode *baseNode, linkNode *newNode);

    bool insertBefore(linkNode *baseNode, linkNode *newNode);

    void removeNode(linkNode *node);

    vector<int> getInOrder();

    linkNode *rootNode{};

    void adjustLabels();

    void rebuildLeftRightLinks(linkNode *node);

    lType getLabelAfter(linkNode *baseNode);

    lType getLabelBefore(linkNode *baseNode);

    ~myList() {
        for (int i = 0; i < linkNodes.size(); ++i) {
            delete linkNodes[i];
        }
    }
};


#define TOPO_MYLIST_H
#endif //TOPO_MYLIST_H
