#include "TreeIterator.h"

// imports from bpp-core
using namespace bpp;
using namespace std;


/******************************************************************************/

void TreeIterator::init()
{
    auto nodes= tree_.getAllNodes();
    for (const auto& node:nodes) 
    {
        nodeToVisited_[tree_.getNodeIndex(node)] = false;
        nodeToSonVisited_[tree_.getNodeIndex(node)] = false;
    }
}

/******************************************************************************/

shared_ptr<PhyloNode> TreeIterator::begin()
{
    nodeToVisited_[tree_.getNodeIndex(currNode_)] = true;
    if (tree_.hasFather(currNode_))
    {
        unsigned int fatherNodeIndex = tree_.getNodeIndex(tree_.getFatherOfNode(currNode_));
        nodeToSonVisited_[fatherNodeIndex] = true;
        std::vector<unsigned int> sonsIndices = tree_.getSons(fatherNodeIndex);
        for (size_t i=0; i<sonsIndices.size(); ++i) 
        {
            if (tree_.getNodeIndex(currNode_) == sonsIndices[i])
                nodeToLastVisitedSonIndex_[fatherNodeIndex] = i;
        }
    }
    return currNode_;
}

/******************************************************************************/

TreeIterator& TreeIterator::operator++()    // prefix ++
{
    this->next();
    return *this;
}

/******************************************************************************/

shared_ptr<PhyloNode> PostOrderTreeIterator::getLeftMostPredecessor(shared_ptr<PhyloNode> startNode)
{
    size_t nextPos = 0;
    if (nodeToSonVisited_[tree_.getNodeIndex(startNode)])
        nextPos = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(startNode)] + 1;
    // if startNode has no predecessors ->  move to him
    if (nextPos >= tree_.getNumberOfSons(startNode))
    {
        return startNode;
    }

    shared_ptr<PhyloNode> node = tree_.getSons(startNode)[nextPos];
    while (tree_.getNumberOfSons(node) > 0)
    {
        node = tree_.getSons(node)[0];
    }
    return node;
}

/******************************************************************************/

shared_ptr<PhyloNode> PostOrderTreeIterator::next()
{   // order: (left, right, parent)

    // stop condition: currNode_ is root
    if (!tree_.hasFather(currNode_))
    {
        return NULL;
    }

    // by the time you visit currNode_, all the nodes in its subtree were already visited
    shared_ptr<PhyloNode> father = tree_.getFatherOfNode(currNode_);
    size_t numOfBrothers =  tree_.getNumberOfSons(father); // get the number of brothers of currNode_
    size_t lastVisitedBrotherPos = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)]; // since currNode is already visited, its father must have at least one visited son (which is currNode_)
    // if all the brothers were already visited -> now visit the father
    if (lastVisitedBrotherPos == numOfBrothers-1)
    {
        currNode_ = father;
    }
        // else -> visit the leftmost predecessor next brother on the list
    else
    {
        size_t nextSiblingPos = lastVisitedBrotherPos + 1;
        currNode_ = tree_.getSons(father)[nextSiblingPos]; // the next brother is not visited yet if its has a non-visited leftmost predecessor
        currNode_ = getLeftMostPredecessor(currNode_);
    }

    // update the returned node to be visited and its father's last visited son, if needed
    nodeToVisited_[tree_.getNodeIndex(currNode_)] = true;
    if (tree_.hasFather(currNode_))
    {
        father = tree_.getFatherOfNode(currNode_);
        if (!nodeToSonVisited_[tree_.getNodeIndex(father)])
        {
            nodeToSonVisited_[tree_.getNodeIndex(father)] = true;
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = 0;
        }
        else
        {
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] + 1;
        }
    }
    return currNode_;
}

/******************************************************************************/

shared_ptr<PhyloNode> PreOrderTreeIterator::next()
{   // order: (parent, left, right)

    vector<unsigned int> leafIds = tree_.getAllLeavesIndexes();
    bool hasVisitedSons = nodeToSonVisited_[tree_.getNodeIndex(currNode_)];
    size_t numOfSons = tree_.getNumberOfSons(currNode_);

    // stop condition: the node is the rightmost leaf of the tree
    if (tree_.getNodeIndex(currNode_) == leafIds[leafIds.size()-1])
    {
        return NULL;
    }

    // by the time you visit currNode_, all its ancestors and left brothers (if exist) were already visited
    if (!hasVisitedSons && numOfSons > 0)
    {
        currNode_ = tree_.getSons(currNode_)[0]; // this somehow modifies the value of nodeToSonVisited_[tree_.getNodeIndex(father)] to true... how?
    }
        // as long as there are still sons to visit -> visit the leftmost unvisited son
    else if (hasVisitedSons && nodeToLastVisitedSonIndex_[tree_.getNodeIndex(currNode_)] < numOfSons-1)
    {
        // if this point is reached, the current node already has visited sons, or it has no sons at all and we are done
        currNode_ = tree_.getSons(currNode_)[nodeToLastVisitedSonIndex_[tree_.getNodeIndex(currNode_)]+1];
        // else -> traverse to the leftmost brother which is right to currNode_ (also occurs when currNode has no sons
    }
    else
    {
        currNode_ = tree_.getFatherOfNode(currNode_);
        size_t lastVisitedSonPos = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(currNode_)]; // the father of the original currNode_ must have at least one visited child which is the original currNode_
        numOfSons = tree_.getNumberOfSons(currNode_);
        while (lastVisitedSonPos == numOfSons-1)
        {
            currNode_ = tree_.getFatherOfNode(currNode_); // the father node must have visited sons as the currNode is a visited son of his
            lastVisitedSonPos = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(currNode_)];
            numOfSons = tree_.getNumberOfSons(currNode_);
        }
        currNode_ = tree_.getSons(currNode_)[lastVisitedSonPos+1];    // need to remove +1??
    }
    // update the returned node to be visited and its father's last visited son, if needed
    if (tree_.hasFather(currNode_))
    {
        shared_ptr<PhyloNode> father = tree_.getFatherOfNode(currNode_);
        if (nodeToSonVisited_[tree_.getNodeIndex(father)])
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] + 1;
        else
        {
            nodeToSonVisited_[tree_.getNodeIndex(father)] = true;
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = 0;
        }
    }
    nodeToVisited_[tree_.getNodeIndex(currNode_)] = true;
    return currNode_;
}

/******************************************************************************/

shared_ptr<PhyloNode> InOrderTreeIterator::doStep(shared_ptr<PhyloNode> node)
{
    // if the node has unvisited left sons -> visit the leftmost unvisited son
    size_t lastVisitedSon;
    if (nodeToSonVisited_[tree_.getNodeIndex(node)])
    {
        lastVisitedSon = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(node)];
    }
    size_t numOfSons = tree_.getNumberOfSons(node);
    bool is_visited =  nodeToVisited_[tree_.getNodeIndex(node)];

    // if the node has unvisited left sons - go to the leftmost unvisited son
    if (!nodeToSonVisited_[tree_.getNodeIndex(node)] && numOfSons > 0)
    {
        return tree_.getSons(node)[0];
    }
    else if (nodeToSonVisited_[tree_.getNodeIndex(node)] && lastVisitedSon < static_cast<size_t>(floor(numOfSons/2)-1))
    {
        return tree_.getSons(node)[lastVisitedSon+1];
    }

    // if the node last visited its last left son or is a not yet visited leaf / parent of a single node -> go to it
    else if (numOfSons > 1 && lastVisitedSon == static_cast<size_t>(floor(numOfSons/2)-1) && !is_visited)
    {
        return node;
    }

    // if the node has unvisited right sons -> go to the leftmost unvisited right son
    else if (nodeToSonVisited_[tree_.getNodeIndex(node)] && lastVisitedSon < numOfSons-1)
    {
        return tree_.getSons(node)[lastVisitedSon+1];
    }

    // else - the entire subtree of the node has been scanned - move to its father
    else
    {
        return tree_.getFatherOfNode(node);
    }
}


/******************************************************************************/

shared_ptr<PhyloNode> InOrderTreeIterator::next()
{   // order: (left (0,..,n/2-1), parent, right (n/2,....,n))

    vector<unsigned int> leafIds = tree_.getAllLeavesIndexes();

    // stop condition: the node is the rightmost leaf of the tree
    if (tree_.getNodeIndex(currNode_) == leafIds[leafIds.size()-1])
    {
        return NULL;
    }

    // while curNode still has unvisited left sons -> do another step
    while (nodeToVisited_[tree_.getNodeIndex(currNode_)]
           || (!nodeToSonVisited_[tree_.getNodeIndex(currNode_)] && !tree_.isLeaf(currNode_)) || (nodeToSonVisited_[tree_.getNodeIndex(currNode_)] && nodeToLastVisitedSonIndex_[tree_.getNodeIndex(currNode_)] < static_cast<size_t>(floor(tree_.getNumberOfSons(currNode_)/2)-1) && !tree_.isLeaf(currNode_)))
    {
        currNode_ = doStep(currNode_);
    }

    // update the returned node to be visited and its father's last visited son, if needed
    if (tree_.hasFather(currNode_))
    {
        shared_ptr<PhyloNode> father = tree_.getFatherOfNode(currNode_);
        if (nodeToSonVisited_[tree_.getNodeIndex(father)])
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] + 1;
        else
        {
            nodeToSonVisited_[tree_.getNodeIndex(father)] = true;
            nodeToLastVisitedSonIndex_[tree_.getNodeIndex(father)] = 0;
        }
    }
    nodeToVisited_[tree_.getNodeIndex(currNode_)] = true;
    return currNode_;
}
