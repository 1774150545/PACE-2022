#include "myList.h"

using namespace std;

void myList::setInitSol(vector<int> &sol) {
    if (sol.empty()) return;

    int firstId = sol[0];
    insertFirst(linkNodes[firstId]);

    for (int i = 1; i < sol.size(); i++) {
        int nodeId = sol[i];
        int lastId = sol[i - 1];
        insertAfter(linkNodes[lastId], linkNodes[nodeId]);
    }

    int mid = len / 2;
    int temp = 0;
    while (temp < mid) {
        rootNode = rootNode->nextNode;
        temp++;
    }

    adjustLabels();
}

void myList::insertFirst(linkNode *baseNode) {
    baseNode->label = LMAX / 2;
    rootNode = baseNode;
    len = 0;
    len++;
}

bool myList::insertAfter(linkNode *baseNode, linkNode *newNode) {
    lType newLabel = getLabelAfter(baseNode);

    if (newLabel == -1) {
        adjustLabels();
        newLabel = getLabelAfter(baseNode);
    }

    if (baseNode->nextNode == nullptr) {
        baseNode->nextNode = newNode;
        newNode->preNode = baseNode;
    } else {
        newNode->nextNode = baseNode->nextNode;
        baseNode->nextNode->preNode = newNode;

        baseNode->nextNode = newNode;
        newNode->preNode = baseNode;
    }
    len++;
    newNode->label = newLabel;

    return true;
}

bool myList::insertBefore(linkNode *baseNode, linkNode *newNode) {
    lType newLabel = getLabelBefore(baseNode);

    if (newLabel == -1) {
        adjustLabels();
        newLabel = getLabelBefore(baseNode);
    }

    if (baseNode->preNode == nullptr) {
        baseNode->preNode = newNode;
        newNode->nextNode = baseNode;
    } else {
        baseNode->preNode->nextNode = newNode;
        newNode->preNode = baseNode->preNode;
        baseNode->preNode = newNode;
        newNode->nextNode = baseNode;
    }

    len++;
    newNode->label = newLabel;
    return true;
}

void myList::adjustLabels() {
    lType avg = LMAX / (len + 1);

    linkNode *leftMost = rootNode;
    while (leftMost->preNode != nullptr) {
        leftMost = leftMost->preNode;
    }

    linkNode *tempNode = leftMost;
    int count = 1;
    while (tempNode != nullptr) {
        tempNode->label = avg * count;
        tempNode = tempNode->nextNode;
        count++;
    }
}

void myList::rebuildLeftRightLinks(linkNode *node) {
    bool hasLeft = node->preNode != nullptr;
    bool hasRight = node->nextNode != nullptr;

    if (hasLeft) {
        node->preNode->nextNode = node->nextNode;
    }
    if (hasRight) {
        node->nextNode->preNode = node->preNode;
    }
}

void myList::removeNode(linkNode *node) {
    len--;

    if (node == rootNode) {
        if (node->preNode != nullptr)
            rootNode = node->preNode;
        else if (node->nextNode != nullptr)
            rootNode = node->nextNode;
        else
            rootNode = nullptr;
    }

    rebuildLeftRightLinks(node);

    node->preNode = nullptr;
    node->nextNode = nullptr;

    node->label = -1;
}

vector<int> myList::getInOrder() {
    vector<int> ordered;

    linkNode *leftMost = rootNode;
    while (leftMost->preNode != nullptr) {
        leftMost = leftMost->preNode;
    }

    linkNode *tempNode = leftMost;
    while (tempNode != nullptr) {
        ordered.push_back(tempNode->nodeId);
        tempNode = tempNode->nextNode;
    }

    return ordered;
}

lType myList::getLabelAfter(linkNode *baseNode) {
    if (baseNode->nextNode == nullptr) {
        if (LMAX - baseNode->label < 4) {
            return -1;
        }
        return baseNode->label / 2 + LMAX / 2;
    } else {
        if (baseNode->nextNode->label - baseNode->label < 4) {
            return -1;
        }
        return baseNode->nextNode->label / 2 + baseNode->label / 2;
    }
}

lType myList::getLabelBefore(linkNode *baseNode) {
    if (baseNode->preNode == nullptr) {
        if (baseNode->label - 0 < 4) {
            return -1;
        }

        return baseNode->label / 2;
    } else {
        if (baseNode->label - baseNode->preNode->label < 4) {
            return -1;
        }
        return baseNode->label / 2 + baseNode->preNode->label / 2;
    }
}
