#include "StochasticMapping.h"
#include "../Tree/PhyloBranch.h"
#include "../Simulation/MutationProcess.h"
#include "RewardMappingTools.h"
#include "Reward.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "../TreeIterator.h"
#include "../Io/Nhx.h"
#include "../Io/Newick.h"
#include "../Model/RateDistribution/ConstantRateDistribution.h"
#include "../Model/MarkovModulatedSubstitutionModel.h"
#include "../Likelihood/DRHomogeneousTreeLikelihood.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric> // to sum over items in a vector

using namespace bpp;
using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, unsigned int numOfMappings) :
  likelihood_(drl),
  tree_(make_shared<const PhyloTree>(likelihood_.get()->getSubstitutionProcess().getParametrizablePhyloTree())),
  fractionalProbabilities_(),
  ConditionalProbabilities_(),
  numOfMappings_(numOfMappings)
{
  ComputeConditionals();
}

/******************************************************************************/

vector<TreeMapping> StochasticMapping::generateStochasticMappings()
{
  vector<TreeMapping> mappings;
  mappings.clear();
  mappings.resize(numOfMappings_);
  for (unsigned int i = 0; i < numOfMappings_; ++i)
  {
    TreeMapping mapping;

    // simulate a set of ancestral states, based on the fractional likelihoods from step 1
    sampleAncestrals(mapping);

    // simulate mutational history of each lineage of the phylogeny, conditional on the ancestral states
    sampleMutationsGivenAncestrals(mapping);

    // add the mapping to the vector of mapping
    mappings[i] = mapping;
  }
  return mappings;
}

/******************************************************************************/

double StochasticMapping::getDistance(const TreeMapping& mapping1, const TreeMapping& mapping2)
{
    double distance = 0;
    auto nodes = tree_.get()->getAllNodes();
    for (auto node: nodes)
    {
        if (tree_.get()->hasFather(node))
        {
            unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(node));
            distance += getPathsDistance(mapping1.at(branchId), mapping2.at(branchId));
        }
    }
    distance /= static_cast<double>(nodes.size());
    return distance;
}

/******************************************************************************/

double StochasticMapping::getPathsDistance(const MutationPath& path1, const MutationPath& path2)
{
    // both paths must apply to the same branch and thus have the same length
    assert(path1.getTotalTime() == path2.getTotalTime());
    double totalTime = path1.getTotalTime();
    
    // get paths info
    vector<size_t> path1States = path1.getStates();
    vector<double> path1Times = path1.getTimes();
    vector<size_t> path2States = path2.getStates();
    vector<double> path2Times = path2.getTimes();

    // init
    size_t path1CurrEvent = 0, path2CurrEvent = 0;
    size_t path1CurrState = path1.getInitialState(), path2CurrState = path2.getInitialState();
    double path1CurrTime = path1Times[path1CurrEvent], path2CurrTime = path2Times[path2CurrEvent];
    double time = 0, distance = 0;


    while (time < totalTime)
    {
        // increment by the lower time value
        if (path1CurrTime < path2CurrTime)
        {
            time += path1CurrTime;
            path1CurrEvent += 1;
            path2CurrTime -= time;
        }
        else
        {
            time += path2CurrTime;
            path2CurrEvent += 1;
            path1CurrTime -= time;
        }

        // add disagreement to distance
        if (path1CurrState != path2CurrState)
            distance += time;

        // update times and states  
        path1CurrTime = path1CurrEvent < path1.getNumberOfEvents() ? path1Times[path1CurrEvent] : (path1.getTotalTime() - accumulate(path1Times.begin(), path1Times.end(), 0.0));
        path1CurrState = path1CurrEvent == 0 ? path1.getInitialState() : path1States[path1CurrEvent];
        path2CurrTime = path2CurrEvent < path2.getNumberOfEvents() ? path2Times[path2CurrEvent] : (path2.getTotalTime() - accumulate(path2Times.begin(), path2Times.end(), 0.0));
        path2CurrState = path2CurrEvent == 0 ? path2.getInitialState() : path2States[path2CurrEvent];
    }
    distance /= totalTime;
    assert(distance <= 1);
    
    return distance;
}

/******************************************************************************/

vector<TreeMapping> StochasticMapping::kMeansClustering(vector<TreeMapping>& mappings, unsigned int k, unsigned int epochs)
{
    // init
    vector<TreeMapping> centroids;
    
    for (unsigned int i = 0; i < k; ++i) {
        centroids.push_back(mappings[RandomTools::giveIntRandomNumberBetweenZeroAndEntry(mappings.size()-1)]);
    }
    map <unsigned int, pair<unsigned int, double>> mappingToCentroid; // maps mapping index (in mappings) to its centroid index (in centroids) and its distance from it

    unsigned it = 0;
    while (it < epochs)
    {
        // assign mappings centroids
        for (unsigned int m=0; m<mappings.size(); ++m)
        {
            unsigned int clusterId = 0;
            unsigned int closestCentroid = clusterId;
            double minDist = getDistance(mappings[m], centroids[clusterId]);
            while (clusterId < k)
            {
                double distance = getDistance(mappings[m], centroids[clusterId]);
                if (distance < minDist)
                {
                    closestCentroid = clusterId;
                    minDist = distance;
                }
            }
            mappingToCentroid[m] = make_pair(closestCentroid, minDist);
        }

        // compute the mean of each cluster based in its members
        for (size_t clusterId=0; clusterId<k; ++clusterId)
        {
            vector<TreeMapping> clusterMembers;
            for (unsigned int m=0; m<mappings.size(); ++m)
            {
                if (mappingToCentroid.at(m).first == clusterId)
                    clusterMembers.push_back(mappings[m]);
            }
            centroids[clusterId] = generateExpectedMapping(clusterMembers);

        }
        mappingToCentroid.clear();
    }

    return centroids;
}

/******************************************************************************/

vector<TreeMapping> StochasticMapping::generatedClusteredMappings(unsigned int numOfClusters)
{
  // init
  vector<TreeMapping> mappingsCentroids;
  mappingsCentroids.clear();
  mappingsCentroids.resize(numOfClusters);

  // compute pairwise distances - can I use additivity here?
  vector<TreeMapping> mappings = generateStochasticMappings();

  // apply clustering (for now, use k-means clustering)
  kMeansClustering(mappings, numOfClusters);

  return mappingsCentroids;
}

/******************************************************************************/

void StochasticMapping::setExpectedAncestrals(TreeMapping& expectedMapping, VVDouble& ancestralStatesFrequencies)
{
  auto nodes= tree_.get()->getAllNodes();
  for (const auto& node:nodes)
  {
    unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
    auto d = distance(ancestralStatesFrequencies[nodeIndex].begin(), max_element(ancestralStatesFrequencies[nodeIndex].begin(), ancestralStatesFrequencies[nodeIndex].end()));
    unsigned int state = static_cast<unsigned int>(d); 
    setNodeState(node, state, expectedMapping); // in the case of a leaf, the assigned state must be sampled
  }
}

/******************************************************************************/

TreeMapping StochasticMapping::generateExpectedMapping(const vector<TreeMapping>& mappings)
{
    // initialize the expected history
    TreeMapping expectedMapping;

    // compute a vector of the posterior asssignment probabilities for each inner node
    VVDouble ancestralStatesFrequencies;
    // initialize the vector
    ancestralStatesFrequencies.clear();
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();
    ancestralStatesFrequencies.resize(tree_.get()->getNumberOfNodes(), VDouble(modelStatesNum));

    computeStatesFrequencies(ancestralStatesFrequencies, mappings);

    // set the ancestral states accrdonig to the maximal posterior (i.e, conditional) probability
    setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);

    // update the expected history with the dwelling times
    PostOrderTreeIterator treeIt(*tree_.get());
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next()) // traverse the tree in post-order
    {
        if (tree_.get()->hasFather(node)) // for any node except to the root
        {
            double branchLength = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(node)).get()->getLength();
            unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(node));        
            // initialize vector of average dwelling times for the branch stemming from node
            VDouble AverageDwellingTimes;
            AverageDwellingTimes.clear();
            AverageDwellingTimes.resize(modelStatesNum, 0);
            // compute the average dwelling times of all the states
            for (size_t h = 0; h < mappings.size(); ++h)
            {
                MutationPath branchMapping = mappings[h].at(branchId);   
                vector<size_t> states = branchMapping.getStates();
                vector<double> times = branchMapping.getTimes();
                assert(branchMapping.getTotalTime() == branchLength);
                shared_ptr<PhyloNode> father = tree_.get()->getFatherOfNode(node);
                AverageDwellingTimes[branchMapping.getInitialState()] += times[0];
                for (size_t e=0; e<branchMapping.getNumberOfEvents()-1; ++e)
                {
                    AverageDwellingTimes[states[e]] += times[e+1]; // only model states, which are non-negative, are considered
                }
            }
            for (size_t state = 0; state < modelStatesNum; ++state)
            {
                AverageDwellingTimes[state] /= static_cast<double>(mappings.size());
            }
            // break the branch according to average dwelling times
            updateBranchByDwellingTimes(expectedMapping, node, AverageDwellingTimes, ancestralStatesFrequencies);
        }
    }
    return expectedMapping;
}

/******************************************************************************/

VVDouble StochasticMapping::getPosteriorProbabilities()
{
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();
    VVDouble posteriorProbabilities;
    posteriorProbabilities.clear();
    posteriorProbabilities.resize(tree_.get()->getNumberOfNodes(), VDouble(modelStatesNum));
    double nodeDataProb = 0;
    PostOrderTreeIterator treeIt(*tree_.get());
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next()) // traverse the tree in post-order
    {
        unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
        nodeDataProb = 0;
        for (size_t s = 0; s < modelStatesNum; ++s)
        {
            nodeDataProb = nodeDataProb + fractionalProbabilities_[nodeIndex][s];
        }
        for (size_t nodeState = 0; nodeState < modelStatesNum; ++nodeState)
        {
            posteriorProbabilities[nodeIndex][nodeState] = fractionalProbabilities_[nodeIndex][nodeState] / nodeDataProb;
        }
    }
    return posteriorProbabilities;
}

/******************************************************************************/

VVDouble StochasticMapping::getAnalyticDwellingTimes()
{
    // auxiliary variables
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();
    const Alphabet* alphabet = likelihood_.get()->getData()->getAlphabet();
    const TransitionModel* model = dynamic_cast<const TransitionModel*>(likelihood_->getSubstitutionProcess().getModel(0, 0)); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    auto nodes = tree_.get()->getAllNodes();

    // Compute the reward per state per site - expect the entries per site to be equal to the number of model states
    // we denote the reward for staying at model state k as rk
    // first, create a numeric alphabet whose states correspond to the model states
    UserAlphabetIndex1 alpha(alphabet);
    map <size_t, int> alphabetStatesToModelStates;
    size_t alphabetStatesNum = alphabet->getNumberOfStates();
    vector<string> resolvedStates = alphabet->getResolvedChars(); // here, only treat resolved states of the alphabet
    
    // Compute the expected dwelling times per branch and state as follows:
    // For branch b of length t, the average dwelling time in state k is rk*t (based on Minin and Suchard paper).
    // The average dwelling times of all the model states should sum uo to t (make sure of it!)
    VVDouble expectedDwellingTimes;
    expectedDwellingTimes.clear();
    expectedDwellingTimes.resize(tree_.get()->getNumberOfNodes(), VDouble(modelStatesNum));
    
    for (size_t s = 0; s<resolvedStates.size(); ++s)
    {
        int character_state_a = alphabet->getState(resolvedStates[s]).getNum();
        if (character_state_a < 0) // special case for gaps - treat as unknown (i.e., the last state in the alphabet) - this section should not be visited if only resolved states are regarded
            character_state_a = alphabet->getStateAt(alphabetStatesNum-1).getNum();
        
        vector<int> alphabetCorrespondingStates_a = alphabet->getAlias(character_state_a);
        for (size_t cs=0; cs<alphabetCorrespondingStates_a.size(); ++cs)
        {
            int state_a = alphabetCorrespondingStates_a[cs];
            alpha.setIndex(state_a, 1); // set the reward of the state as 1 and the reward for the rest of the states as 0
            //for (size_t m = 0; m < alphabetStatesNum; ++m)
            for (size_t m=0; m<resolvedStates.size(); ++m)
            {
                int character_state_b = alphabet->getState(resolvedStates[m]).getNum();
                // special case for gaps - treat as unknown (i.e., the last state in the alphabet)
                if (character_state_b < 0)
                    character_state_b = alphabet->getStateAt(alphabetStatesNum-1).getNum();
                vector<int> alphabetCorrespondingStates_b = alphabet->getAlias(character_state_b);
                for (size_t cm=0; cm<alphabetCorrespondingStates_b.size(); ++cm)
                {
                    int state_b = alphabetCorrespondingStates_b[cm];
                    bool in_a = false;
                    for (size_t i=0; i<alphabetCorrespondingStates_a.size(); ++i)
                    {
                        if (state_b == alphabetCorrespondingStates_a[i])
                            in_a = true;
                    }
                    if (!in_a)
                    {
                        alpha.setIndex(state_b, 0);
                    }
                }
            }
            DecompositionReward reward(dynamic_cast<const SubstitutionModel*>(model), alpha.clone()); // TO FIX 20.6: this line attempts to delete alpha which doesn't belong to it. cloning it didn't help - get help from Itay / Anat
            //LikelihoodCalculationSingleProcess* rewardLik = likelihood_.get()->clone();
            shared_ptr<LikelihoodCalculationSingleProcess> rewardLik = make_shared<LikelihoodCalculationSingleProcess>(*likelihood_);
            ProbabilisticRewardMapping mapping(RewardMappingTools::computeRewardVectors(*rewardLik.get(), tree_.get()->getAllNodesIndexes(), reward, false));
            
            for (auto node: nodes) 
            {
                if (tree_.get()->hasFather(node)) // for any node except to the root
                {
                    unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
                    double dwellingTime =  mapping.getReward(nodeIndex, 0);
                    vector <size_t> correspondingModelStates = model->getModelStates(state_a);
                    for (size_t ms=0; ms<correspondingModelStates.size(); ++ms)
                        expectedDwellingTimes[nodeIndex][correspondingModelStates[ms]] = dwellingTime / static_cast<double>(correspondingModelStates.size());
                }
            }
        }
    }

    // standardize expected dwelling times, if needed
    double sumOfDwellingTimes;
    for (auto node: nodes) 
    {
        if (tree_.get()->getNodeIndex(node)) // for any node except to the root
        {
            unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
            double branchLength = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(node)).get()->getLength();
            sumOfDwellingTimes = 0;
            for (size_t s = 0; s < modelStatesNum; ++s)
            {
                sumOfDwellingTimes = sumOfDwellingTimes + expectedDwellingTimes[nodeIndex][s];
            }
            if (sumOfDwellingTimes != branchLength)
			{
				for (size_t s = 0; s < modelStatesNum; ++s)
				{
					expectedDwellingTimes[nodeIndex][s] =  branchLength * (expectedDwellingTimes[nodeIndex][s]) / sumOfDwellingTimes;
				}
			}
        }
    }

    return expectedDwellingTimes;
}

/******************************************************************************/

TreeMapping StochasticMapping::generateAnalyticExpectedMapping()
{
    // initialize
    TreeMapping expectedMapping;
    auto nodes = tree_.get()->getAllNodes();

    // Compute the posterior assignment probabilities to internal nodes, based on the fractional probablities computed earlier
    // because the sum of partial likelihoods (i.e, the fractional probabilities) is in fact the probablity of the data, 
    //it is sufficient to standardize the vector of fractional probabilires for each node to obtain the posterior probabilities
    VVDouble posteriorProbabilities = getPosteriorProbabilities();

    // Assign states to internal nodes based on the majority rule over the posterior probabilities
    setExpectedAncestrals(expectedMapping, posteriorProbabilities);

    // update the mapping with the analytically expected dwelling times
    VVDouble expectedDwellingTimes = getAnalyticDwellingTimes();
    for (auto node: nodes) 
    {
        if (tree_.get()->hasFather(node)) // for any node except to the root
        {
            updateBranchByDwellingTimes(expectedMapping, node, expectedDwellingTimes[tree_.get()->getNodeIndex(node)], posteriorProbabilities);
        }
    }

    return expectedMapping;
}

/******************************************************************************/

unsigned int StochasticMapping::getNodeState(shared_ptr<PhyloNode> node, const TreeMapping& mapping)
{
  if (tree_.get()->isLeaf(node))
  {
      return getLeafCharacterState(node.get()->getName());
  }
  vector<shared_ptr<PhyloNode>> sons = tree_.get()->getSons(node);
  unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(sons[0]));
  return static_cast<unsigned int>(mapping.at(branchId).getInitialState());
}

/******************************************************************************/

void StochasticMapping::setNodeState(shared_ptr<PhyloNode> node, unsigned int state, TreeMapping& mapping)
{
    vector<shared_ptr<PhyloNode>> sons = tree_.get()->getSons(node);
    for (auto son: sons)
    {
        unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(son));
        double branchLength = tree_.get()->getEdgeToFather(son).get()->getLength();
        MutationPath branchMapping(likelihood_->getData()->getAlphabet(), state, branchLength); // create a path that starts an node and ends at son of length branchLength
        mapping.insert(pair<int,MutationPath>(branchId, branchMapping));
    }
}

/******************************************************************************/

vector<unsigned int> StochasticMapping::getLeafModelStates(PhyloNode& node)
{
    unsigned int leafAlphabetState = getLeafCharacterState(node.getName());
    const TransitionModel* model = dynamic_cast<const TransitionModel*>(likelihood_->getSubstitutionProcess().getModel(0, 0)); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    vector<int> leafCharacterStates = likelihood_.get()->getData()->getAlphabet()->getAlias(leafAlphabetState);
    vector<unsigned int> leafModelStates;
    leafModelStates.clear();
    for (size_t s=0; s<leafCharacterStates.size(); ++s)
    {
        vector <size_t> correspondingModelStates = model->getModelStates(leafCharacterStates[s]);
        for (size_t m=0; m<correspondingModelStates.size(); ++m)
        {
            leafModelStates.push_back(static_cast<unsigned int>(correspondingModelStates[m]));
        }
    }
    return leafModelStates; 
}

/******************************************************************************/

unsigned int StochasticMapping::getLeafCharacterState(string leafName)
{
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();
    int leafState = -1;
    for (unsigned int state=0; state<modelStatesNum; ++state)
    {
      if (likelihood_->getData()->getStateValueAt(0, leafName, state) > 0)
      {
        leafState = state;
        break;
      }
    }
    if (leafState < 0) // if leaf is not mapped to a concrete model state (occurs if its a gap or an unknown character)
    {
        VDouble baseDistribution;
        baseDistribution.clear(); 
        baseDistribution.resize(modelStatesNum);
        for (size_t s=0; s<modelStatesNum; ++s)
            baseDistribution[s] = 1./static_cast<double>(modelStatesNum);
        leafState = sampleState(baseDistribution);
    }
    return static_cast<unsigned int>(leafState);
}

/******************************************************************************/

void StochasticMapping::computeFractionals()
{
    // some auxiliiary variables
    const TransitionModel* model = dynamic_cast<const TransitionModel*>(likelihood_->getSubstitutionProcess().getModel(0, 0)); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();

    // compute the fractional probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s]
    PostOrderTreeIterator treeIt(*tree_.get());
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next()) // traverse the tree in post-order
    {
        unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
        if (tree_.get()->isLeaf(node)) // if the node is a leaf, set the fractional probability of its state to 1, and the rest ot 0
        {
            vector<unsigned int> leafModelStates = getLeafModelStates(*node.get());
            for (unsigned int s = 0; s < leafModelStates.size(); ++s)
            {
                fractionalProbabilities_[nodeIndex][leafModelStates[s]] = 1;
            }
        }
    else                // if the node is internal, follow the Felesenstein computation rule to compute the fractional probability
    {
      std::vector<unsigned int> sonsIndices = tree_.get()->getSons(nodeIndex);
      for (unsigned int nodeState = 0; nodeState < modelStatesNum; ++nodeState)
      {
        double fullProb = 1;
        for (unsigned int j = 0; j < tree_.get()->getNumberOfSons(node); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
        {
          unsigned int sonIndex = sonsIndices[j];
          auto son = tree_.get()->getNode(sonIndex);
          double sonProb = 0;
          double bl = tree_.get()->getEdgeToFather(sonIndex).get()->getLength();
          for (unsigned int sonState = 0; sonState < modelStatesNum; ++sonState)
          {
            sonProb += model->Pij_t(nodeState, sonState, bl) * fractionalProbabilities_[sonIndex][sonState];
          }
          fullProb *= sonProb;
        }
        fractionalProbabilities_[nodeIndex][nodeState] = fullProb;
      }
    }
  }
}

/******************************************************************************/

void StochasticMapping::ComputeConditionals()
{
    // some auxiliiary variables
    const TransitionModel* model = dynamic_cast<const TransitionModel*>(likelihood_->getSubstitutionProcess().getModel(0, 0)); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();

    // compute the fractional probabilities for each node and state
    fractionalProbabilities_.clear();
    fractionalProbabilities_.resize(tree_.get()->getNumberOfNodes(), VDouble(modelStatesNum));
    computeFractionals();

    //  compute the conditional probabilities: for each combination of nodes son, father, compute Pr(son recieves sonState | father has fatherState)
    ConditionalProbabilities_.clear();
    ConditionalProbabilities_.resize(tree_.get()->getNumberOfNodes(), VVDouble(modelStatesNum, VDouble(modelStatesNum)));

    PostOrderTreeIterator treeIt(*tree_.get());
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next())
    {
        unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
        if (!tree_.get()->isLeaf(node) || !tree_.get()->hasFather(node))  // the second condition will catch the root even if it has a single child (in which case, isLeaf() returns true)
        {
            if (!tree_.get()->hasFather(node))  // if the node is the root -> set the conditional probability to be same for all "fatherStates"
            {
                double sum = 0.0;
                for (size_t sonState = 0; sonState < modelStatesNum; ++sonState)
                {
                    double stateConditionalNominator = fractionalProbabilities_[nodeIndex][sonState] * model->freq(sonState);
                    for (size_t fatherState = 0; fatherState < modelStatesNum; ++fatherState)
                    {
                        ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
                    }
                    sum += stateConditionalNominator;
                }
                for (size_t fatherState = 0; fatherState < modelStatesNum; ++fatherState)
                {
                    for (size_t sonState = 0; sonState < modelStatesNum; ++sonState)
                    {
                        ConditionalProbabilities_[nodeIndex][fatherState][sonState] /= sum;
                    }
                }
            }
            else                            // else -> follow equation (10) from the paper to compute the consitional assingment probabilities given the ones of his father
            {
                for (size_t fatherState = 0; fatherState < modelStatesNum; ++fatherState)
                {
                    double sum = 0.0;
                    for (size_t sonState = 0; sonState < modelStatesNum; ++sonState)
                    {
                        double bl = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(node)).get()->getLength();
                        double stateConditionalNominator = fractionalProbabilities_[nodeIndex][sonState] * model->Pij_t(fatherState, sonState, bl);
                        ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
                        sum += stateConditionalNominator;
                    }
                    for (size_t sonState = 0; sonState < modelStatesNum; ++sonState)
                    {
                        ConditionalProbabilities_[nodeIndex][fatherState][sonState] /= sum;
                    }
                }
            }
        }
        else // the node is a leaf, so the conditional probabilities should correspond to the possible leaf states assignments
        {
            vector<unsigned int> leafModelStates = getLeafModelStates(*node.get());  // gap integer is negative and cannot be cast to size_t instance
            for (size_t fatherState = 0; fatherState < modelStatesNum; ++fatherState)
            {
                double sum = 0.0;
                for (size_t s = 0; s < leafModelStates.size(); ++s)
                {
                    size_t sonState = leafModelStates[s];
                    double bl = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(node)).get()->getLength();
                    double stateConditionalNominator = fractionalProbabilities_[nodeIndex][sonState] * model->Pij_t(fatherState, sonState, bl);
                    ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
                    sum += stateConditionalNominator;
                }
                for (size_t s = 0; s < leafModelStates.size(); ++s)
                {
                    ConditionalProbabilities_[nodeIndex][fatherState][leafModelStates[s]] /= sum;
                }
            }
        }
    }
}

/******************************************************************************/

void StochasticMapping::computeStatesFrequencies(VVDouble& ancestralStatesFreuquencies, const vector<TreeMapping>& mappings)
{
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();
    PreOrderTreeIterator treeIt(*tree_.get());
    
    // compute the node assignment probabilities based on their frequency in the mappings
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next())
    {
        size_t nodeIndex = tree_.get()->getNodeIndex(node);
        // go over all the mappings and collect the number of states assignment per node (inclusing leafs which can have ambiguous states like N, in which case they will not be consistent across mappings)
        fill(ancestralStatesFreuquencies[nodeIndex].begin(), ancestralStatesFreuquencies[nodeIndex].end(), 0); // reset all the values to 0
        for (size_t h = 0; h < mappings.size(); ++h)
        {
            ancestralStatesFreuquencies[nodeIndex][getNodeState(node, mappings[h])]++; // the node index is non-negative value mapped to the node id, which could be negative
        }
        // now divide the vector entries by the number of mappings
        for (size_t nodeState = 0; nodeState < modelStatesNum; ++nodeState)
        {
            ancestralStatesFreuquencies[nodeIndex][nodeState] = ancestralStatesFreuquencies[nodeIndex][nodeState] / static_cast<int>(mappings.size());
        }
    }
}

/******************************************************************************/

unsigned int StochasticMapping::sampleState(const VDouble& distibution)
{
    unsigned int state = 0;        // the default state is 0
    double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

    for (unsigned int i = 0; i < distibution.size(); ++i)
    {
        prob -= distibution[i];
        if (prob < 0)  // if the the sampled probability is smaller than the probability to choose state i -> set state to be i
        {
            state = i;
            break;
        }
    }
    return state;
}

/******************************************************************************/

void StochasticMapping::sampleAncestrals(TreeMapping& mapping)
{
    PreOrderTreeIterator treeIt(*tree_.get());
    for (auto node = treeIt.begin(); node != treeIt.end(); node = treeIt.next())
    {
        size_t nodeIndex = tree_.get()->getNodeIndex(node);
        if (!tree_.get()->isLeaf(node))
        {
            if (!tree_.get()->hasFather(node))
            {
                unsigned int rootState = sampleState(ConditionalProbabilities_[nodeIndex][0]); // set father state to 0 (all the entries in the fatherState level are the same anyway)
                setNodeState(node, rootState, mapping);
            }
            else
            {
                unsigned int fatherState = getNodeState(tree_.get()->getFatherOfNode(node), mapping);
                unsigned int sonState = sampleState(ConditionalProbabilities_[nodeIndex][fatherState]);
                setNodeState(node, sonState, mapping);
            }
        }
        else // the node is a leaf
        {
            vector<unsigned int> leafStates = getLeafModelStates(*node.get()); // get the model states corresponding to the node
            if (leafStates.size() == 1)
            {
                setNodeState(node, leafStates[0], mapping);
            }
            else // in the case of a leaf wih multiple possible assignments, the assigned state must be sampled
            {
                unsigned int fatherState = getNodeState(tree_.get()->getFatherOfNode(node), mapping);
                unsigned int sonState = sampleState(ConditionalProbabilities_[nodeIndex][fatherState]);
                setNodeState(node, sonState, mapping);
            }
        }
    }
}

/******************************************************************************/

 void StochasticMapping::sampleMutationsGivenAncestrals(TreeMapping& mapping)
{
    auto nodes = tree_.get()->getAllNodes();
    for (auto node: nodes) 
    {
        if (tree_.get()->hasFather(node))
        {
            sampleMutationsGivenAncestralsPerBranch(mapping, node);
        }
    }
}

/******************************************************************************/

void StochasticMapping::updateBranchMapping(TreeMapping& mapping, shared_ptr<PhyloNode> son, const MutationPath& sampledMapping)
{
    // if mapping was successfull -> copy tryMapping to branchMapping
    unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(son));
    MutationPath& branchMapping = mapping.at(branchId);
    
    // assert that the sampled mapping is consistent with the original mapping
    assert(branchMapping.getInitialState() == sampledMapping.getInitialState());
    assert(branchMapping.getTotalTime() == sampledMapping.getTotalTime());

    // copy the content of the sampled mutations path to the mapping
    vector<size_t> states = sampledMapping.getStates();
    vector<double> times = sampledMapping.getTimes();
    for (size_t e=0; e<sampledMapping.getNumberOfEvents(); ++e)
    {
        branchMapping.addEvent(states[e], times[e]);
    }
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestralsPerBranch(TreeMapping& mapping, shared_ptr<PhyloNode> son, unsigned int maxIterNum)
{
    shared_ptr<PhyloNode> father = tree_.get()->getFatherOfNode(son);
    unsigned int fatherState = getNodeState(father, mapping);
    unsigned int sonState = getNodeState(son, mapping);

    double branchLength = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(son)).get()->getLength();

    const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(likelihood_->getSubstitutionProcess().getModel(0, 0)); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    SimpleMutationProcess mutationProcess(model);

    // simulate mapping on a branch until you manage to finish at the son's state
    for (size_t i = 0; i < maxIterNum; ++i)
    {
        double disFromNode = 0.0;
        MutationPath tryMapping(likelihood_->getData()->getAlphabet(), fatherState, branchLength);
        unsigned int curState = static_cast<unsigned int>(mutationProcess.mutate(fatherState));               // draw the state to transition to after from initial state curState based on the relative tranistion rates distribution (see MutationProcess.cpp line 50)
            
        double timeTillChange;
        // if the father's state is not the same as the son's state -> use the correction corresponding to equation (11) in the paper
        if (fatherState != sonState)
        {   // sample timeTillChange conditional on it being smaller than branchLength
            double u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
            double waitingTimeParam = -1 * model->Qij(fatherState, fatherState); // get the parameter for the exoponential distribution to draw the waiting time from
            double tmp = u * (1.0 - exp(branchLength * -waitingTimeParam));
            timeTillChange =  -log(1.0 - tmp) / waitingTimeParam;
            assert (timeTillChange < branchLength);
        }
        else
        {
            timeTillChange = mutationProcess.getTimeBeforeNextMutationEvent(fatherState); // draw the time until a transition from exponential distribution with the rate of leaving fatherState
        }

        while (disFromNode + timeTillChange < branchLength)  // a jump occured but not passed the whole branch ->
        {
            tryMapping.addEvent(curState, timeTillChange);                                        // add the current state and time to branch history
            disFromNode += timeTillChange;
            timeTillChange = mutationProcess.getTimeBeforeNextMutationEvent(curState);            // draw the time until a transition from exponential distribution with the rate of leaving curState
            curState = static_cast<unsigned int>(mutationProcess.mutate(curState));               // draw a new state to mutate to
        }
        // the last jump passed the length of the branch -> finish the simulation and check if it's sucessfull (i.e, mapping is finished at the son's state)
        if (curState != sonState) // if the simulation failed, try again
        {
            continue;
        }
        else                      // if the simulation was sucessfull, add it to the built mapping
        {
            if (disFromNode > 0) // only if there is turely an event to add -> add it to the mapping
                tryMapping.addEvent(sonState, branchLength - disFromNode);
            updateBranchMapping(mapping, son, tryMapping);     // add the successfull simulation to the build mapping
            return;
        }
    }

    // if all simulations failed -> throw an exception
    throw Exception("could not produce simulations with father = " + TextTools::toString(fatherState) + " son " + TextTools::toString(sonState) + " branch length = " + TextTools::toString(branchLength));
}

/******************************************************************************/

 void StochasticMapping::updateBranchByDwellingTimes(TreeMapping& mapping, shared_ptr<PhyloNode> node, VDouble& dwellingTimes, VVDouble& ancestralStatesFrequencies)
{

    // first, convert the dwelling times vector to a mutation path of the branch
    size_t modelStatesNum = likelihood_->getSubstitutionProcess().getNumberOfStates();

    shared_ptr<PhyloNode> father = tree_.get()->getFatherOfNode(node);
    unsigned int nodeIndex = tree_.get()->getNodeIndex(node);
    unsigned int fatherIndex = tree_.get()->getNodeIndex(father);
    unsigned int sonState = getNodeState(node, mapping);
    unsigned int fatherState = getNodeState(father, mapping);

    double Pf = 1, Ps = 1;
    double shareOfFather = 0, shareOfSon = 0;

    double branchLength = tree_.get()->getEdgeToFather(tree_.get()->getNodeIndex(node)).get()->getLength();
    unsigned int branchId = tree_.get()->getEdgeIndex(tree_.get()->getEdgeToFather(node));
    MutationPath& branchMapping = mapping.at(branchId);

    // set the first event with the dwelling time that matches the state of the father
    if (fatherState == sonState)
    {
        if (tree_.get()->hasFather(node))
        {
            Pf = ancestralStatesFrequencies[fatherIndex][fatherState];
        }
        Ps = ancestralStatesFrequencies[nodeIndex][sonState];
        shareOfFather = Pf / (Pf + Ps);
        branchMapping.addEvent(fatherState, dwellingTimes[fatherState] * shareOfFather);

    }
    else
    {
        branchMapping.addEvent(fatherState, dwellingTimes[fatherState]);
    }
    // set all events except for the one entering the son
    for (size_t state = 0; state < modelStatesNum; ++state)
    {
        if (state != fatherState && state != sonState && dwellingTimes[state] > 0) // if the state matches an event which is not the first or the last -> add it
        {
            branchMapping.addEvent(state, dwellingTimes[state]);
        }
    }
    // change the length of the branch whose bottom node is the son according to the dwelling time of the relevant state
    if (fatherState == sonState)
    {
        shareOfSon = 1 - shareOfFather;
        branchMapping.addEvent(sonState, dwellingTimes[sonState] * shareOfSon);
    }
    else
    {
        branchMapping.addEvent(sonState, dwellingTimes[sonState]);
    }
    assert(branchMapping.getTotalTime() == branchLength);
}
