// Tester for StochasticMapping implementation

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>

// From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-phyl
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Io/Nhx.h>
// #include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
// #include <Bpp/Phyl/Model/G2001.h>
#include <Bpp/Phyl/Model/Character/SingleRateModel.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
// #include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
// #include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>
// #include <Bpp/Phyl/Mapping/RewardMappingTools.h>
// #include <Bpp/Phyl/Mapping/Reward.h>
// #include <Bpp/Phyl/Mapping/DecompositionReward.h>
#include <Bpp/Phyl/Simulation/MutationProcess.h>

using namespace bpp;
using namespace std;

typedef map<unsigned int, bpp::MutationPath> TreeMapping; // maps branch index to its mutations path


/******************************************************************************/

void checkIfMappingLegal(StochasticMapping& stocMapping, TreeMapping& mapping, PhyloTree& baseTree, LikelihoodCalculationSingleProcess& tl)
{
    auto nodes = baseTree.getAllNodes();
    auto leaves = baseTree.getAllLeaves();
    size_t modelStatesNum = tl.getSubstitutionProcess().getNumberOfStates();
  
    // make sure that the states at the leaves are set properly in mapping
    for (auto leaf: leaves)
    {
        // extract the leaf state from the character data in tl
        int leafState = -1;
        for (unsigned int state=0; state<modelStatesNum; ++state)
        {
            if (tl.getData()->getStateValueAt(0, leaf.get()->getName(), state) > 0)
            {
                leafState = state;
                break;
            }
        }
        if (leafState < 0) // leaf is not mapped to a concrete state ->  this test is redumdant
            continue;
        
        // get the state of the respective leaf from the mapping using its respective branch id
        unsigned int branchId = baseTree.getEdgeIndex(baseTree.getEdgeToFather(leaf));
        MutationPath branchMapping = mapping.at(branchId);
        size_t leafStateInMapping = branchMapping.getInitialState(); // if no events have occured, the father state equals the son (leaf) state
        if (branchMapping.getNumberOfEvents() > 0)
        {
            vector<size_t> states = branchMapping.getStates();
            leafStateInMapping = states[states.size()-1];
        }
        if (static_cast<size_t>(leafState) != leafStateInMapping) // if the node's state corresaponds to a concrete character state (not unknown character or a combination of state and rate in the case of a markov modulated model)
        {
            throw Exception("Leaf state for " + leaf.get()->getName() + " was not properly maintained in mapping and its value is " + TextTools::toString(leafStateInMapping) + " instead of " + TextTools::toString(leafState));
        }
    }
            
    // make sure that branch lengths are maintained in the mapping
    for (auto node: nodes)
    {
        if (baseTree.hasFather(node))
        {
            double branchLength = baseTree.getEdgeToFather(baseTree.getNodeIndex(node)).get()->getLength();
            unsigned int branchId = baseTree.getEdgeIndex(baseTree.getEdgeToFather(node));
            MutationPath branchMapping = mapping.at(branchId);
            double branchLengthInMapping = branchMapping.getTotalTime();
            if (abs(branchLengthInMapping - branchLength) > 0.0001) // expected history fails in the root - but the branch length of the branch coming from the root is insignificant as it isn't included in the tree. what happened here?
            {
                throw Exception("branch lengths for branch (" + node.get()->getName() + ", " + baseTree.getFatherOfNode(node).get()->getName() + ") not maintained in mapping and its value is " + TextTools::toString(branchLengthInMapping) + " instead of " + TextTools::toString(branchLength));
            }
        }
    }
    
    // make sure there are no two nodes in a row such that both don't exist in the base tree and both recieve the same state (indicator of illegal transition)
    for (auto node: nodes)
    {
        if (baseTree.hasFather(node))
        {
            unsigned int branchId = baseTree.getEdgeIndex(baseTree.getEdgeToFather(node));
            MutationPath branchMapping = mapping.at(branchId);
            vector<size_t> states = branchMapping.getStates();
            for (size_t e=0; e<branchMapping.getNumberOfEvents()-1; ++e)
            {
                if (states[e] == states[e+1])
                    throw Exception("illegal transitions in the mapping between two identical states " + TextTools::toString(states[e]) + " at branch (" + node.get()->getName() + ", " + baseTree.getFatherOfNode(node).get()->getName() + ")");
            }
        }
    }
}


int main() 
{
  try
    {
        //fix seed for debugging purposes
        double seedUb = 10000000;
        double seed = RandomTools::giveRandomNumberBetweenZeroAndEntry(seedUb);
        RandomTools::setSeed(static_cast<long int>(seed));
        cout << "seed: " << seed << endl; // for debugging purposes in case the tester fails
        
        // process tree
        Newick reader;
        auto phyloTree = std::unique_ptr<PhyloTree>(reader.parenthesisToPhyloTree("(S15:0.85385,((S19:0.0854569,S16:0.139158):0.248594,(((((S12:0.0215813,S14:0.0122578):0.00733911,((S20:0.0133406,S18:0.02058):0.00622244,S8:0.0616991):0.00855007):0.0194517,S21:0.0361841):0.0260926,(S10:6.01257,S9:0.0572114):0.00432963):0.0582364,((((S11:0.00192042,((S13:0.00546429,S22:0.00413541):0.1,S7:0.00544892):0.00223313):0.0224013,S6:0.0147796):0.012621,(S24:0.1,S23:0.1):0.020303):0.0480321,((S2:0.0212492,((S1:0.029627,S3:0.322449):0.1,S17:0.0303775):0.1):0.0311297,(S5:0.00337913,S4:0.1):0.0451854):0.00880453):0.0445887):0.133367):0.85385);", false, "", false, false));
        ParametrizablePhyloTree* paramPhyloTree = new ParametrizablePhyloTree(*phyloTree);

        // process character data
        const IntegerAlphabet* alpha = new IntegerAlphabet(2);
        VectorSiteContainer sites(alpha);
		sites.addSequence(BasicSequence("S1", "1", alpha));
		sites.addSequence(BasicSequence("S2", "0", alpha));
		sites.addSequence(BasicSequence("S3", "1", alpha));
		sites.addSequence(BasicSequence("S4", "0", alpha));
		sites.addSequence(BasicSequence("S5", "0", alpha));
		sites.addSequence(BasicSequence("S6", "0", alpha));
		sites.addSequence(BasicSequence("S7", "0", alpha));
		sites.addSequence(BasicSequence("S8", "1", alpha));
		sites.addSequence(BasicSequence("S9", "0", alpha));
		sites.addSequence(BasicSequence("S10", "0", alpha));
		sites.addSequence(BasicSequence("S11", "0", alpha));
		sites.addSequence(BasicSequence("S12", "0", alpha));
		sites.addSequence(BasicSequence("S13", "0", alpha));
		sites.addSequence(BasicSequence("S14", "1", alpha));
		sites.addSequence(BasicSequence("S15", "1", alpha));
		sites.addSequence(BasicSequence("S16", "1", alpha));
		sites.addSequence(BasicSequence("S17", "0", alpha));
		sites.addSequence(BasicSequence("S18", "1", alpha));
		sites.addSequence(BasicSequence("S19", "0", alpha));
		sites.addSequence(BasicSequence("S20", "1", alpha));
		sites.addSequence(BasicSequence("S21", "1", alpha));
		sites.addSequence(BasicSequence("S22", "0", alpha));
		sites.addSequence(BasicSequence("S23", "0", alpha));
		sites.addSequence(BasicSequence("S24", "-", alpha));
        SiteContainerTools::changeGapsToUnknownCharacters(sites);

        // create a binary model
        vector<double> freqVals(alpha->getNumberOfStates());
        freqVals[0] = 0.5;
        freqVals[1] = 1-freqVals[0];
        shared_ptr<IntegerFrequencySet> freqs;
        freqs.reset(new FullIntegerFrequencySet(alpha, freqVals)); // frequencies set must be reset upon each usage to avoid concatanation of namespaces of prevous models
        auto model = std::make_shared<SingleRateModel>(alpha, freqs, false);
        model.get()->setParameterValue("global_rate", 42.);
        model.get()->setFrequencySet(*freqs);
        SimpleSubstitutionProcess* process = new SimpleSubstitutionProcess(model, paramPhyloTree);
        
        // create tree likelihood computation process
        Context context;
        auto characterTreeLikelihood = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, *process);
        auto copy = characterTreeLikelihood.get()->clone();
        
        // generate 1000 sotchastic mappings
        cout << "**** Generating 1000 stochastic mappings ****" << endl;
        unsigned int mappingsNum = 1000;
        auto stocMapping = std::make_shared<StochasticMapping>(characterTreeLikelihood, mappingsNum); // it is a good general practice to use "explicit" keyword on constructors with a single argument: https://stackoverflow.com/questions/121162/what-does-the-explicit-keyword-mean
        vector<TreeMapping> mappings = stocMapping->generateStochasticMappings();
        
        // make sure all the mappings are legal
        for (size_t i=0; i<mappingsNum; ++i)
        {
            checkIfMappingLegal(*stocMapping.get(), mappings[i], *phyloTree.get(), *characterTreeLikelihood.get());
        }
        cout << "**** All the generated mappings are legal ****" << endl;

		// compute ancestral frequencies over the stochastic mappings
        VVDouble ancestralFrequencies;
        ancestralFrequencies.clear();
        size_t statesNum = characterTreeLikelihood.get()->getSubstitutionProcess().getNumberOfStates();
        ancestralFrequencies.resize(phyloTree.get()->getNumberOfNodes(), VDouble(statesNum));
        stocMapping->computeStatesFrequencies(ancestralFrequencies, mappings);
        
		// generate an expected history
        cout << "**** Generating mappings-based expected mapping ****" << endl;
        TreeMapping expectedHistory = stocMapping->generateExpectedMapping(mappings);
        checkIfMappingLegal(*stocMapping.get(), expectedHistory, *phyloTree.get(), *characterTreeLikelihood.get());
        cout << "**** The generated expected mapping is legal ****" << endl;

        // make sure ancestral assignments correspond to frequencies in the mappings and the the posterior probabilities
        auto nodes = stocMapping->getBaseTree().get()->getAllNodes();
		for (auto node: nodes)
        {
           if (!stocMapping->getBaseTree().get()->isLeaf(node))
           {
                unsigned int state = stocMapping->getNodeState(node, expectedHistory);
                double stateFrequency = ancestralFrequencies[stocMapping->getBaseTree().get()->getNodeIndex(node)][state];
                if (stateFrequency < 0.5)
                {
                    cout << "Failed to assign ancestral state to node " << node->getName() << " according to the frequency: Assigned state is " << state << " while its frequency is " << stateFrequency << endl;
                    return 1;
                }
           }
        }
        
		// compute the average dwelling times at node S19
        shared_ptr<const PhyloTree> baseTree = stocMapping->getBaseTree();
        VDouble AverageDwellingTimes;
        AverageDwellingTimes.clear();
        AverageDwellingTimes.resize(statesNum, 0);
        unsigned int nodeIndex = 1;
        std::shared_ptr<PhyloNode> node = baseTree.get()->getNode(nodeIndex); // this node index correspond to S19
        unsigned int branchId = baseTree.get()->getEdgeIndex(baseTree.get()->getEdgeToFather(node));
    
		// compute the average dwelling times of all the states
        for (size_t i=0; i<mappings.size(); ++i)
        {
            TreeMapping& mapping =  mappings[i];
            MutationPath branchMapping = mapping.at(branchId);
            vector<double> times = branchMapping.getTimes();
            vector<size_t> states = branchMapping.getStates();
            AverageDwellingTimes[branchMapping.getInitialState()] += times[0];
            for (size_t e=0; e<branchMapping.getNumberOfEvents()-1; ++e)
            {
                AverageDwellingTimes[states[e]] += times[e+1];
            }
        }
        for (size_t s=0; s<statesNum; ++s)
        {
            AverageDwellingTimes[s] = AverageDwellingTimes[s] / static_cast<double>(mappings.size());
        }
        MutationPath expectedBranchMapping = expectedHistory.at(branchId);
        vector<double> expectedTimes = expectedBranchMapping.getTimes();
        double splitToFather = expectedTimes[0];
        double splitFromSon = expectedTimes[expectedBranchMapping.getNumberOfEvents()-1];

		// check state of father of S19: if the state of the father is 0 -> also make sure the division of dwelling time under state 0 corresponds to the frequency of 0 in the father
        unsigned int nodeState = stocMapping->getNodeState(node, expectedHistory);
        std::shared_ptr<PhyloNode> father = baseTree.get()->getFatherOfNode(node);
        unsigned int fatherIndex = baseTree.get()->getNodeIndex(father);
        unsigned int fatherState = stocMapping->getNodeState(father, expectedHistory);
        if (fatherState == nodeState)
        {
            // compute division of dwelling time according to the states freuqncies at the father
            double fatherFrequency = ancestralFrequencies[fatherIndex][nodeState];
			double sonFrequency = ancestralFrequencies[nodeIndex][nodeState];
            double fatherShare = fatherFrequency / (sonFrequency+fatherFrequency) * AverageDwellingTimes[nodeState];
            double sonShare = AverageDwellingTimes[nodeState] - fatherShare;

            if (abs(splitFromSon - sonShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch of son is of length " << splitFromSon << " instead of " << sonShare << endl;
                return 1;
            }
            if (abs(splitToFather - fatherShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch beneath father is of length " << splitToFather << " instead of " << fatherShare << endl;
                return 1;               
            }
        }
        else
        {
            if (abs(splitFromSon - AverageDwellingTimes[nodeState]) > 0.0001)
            {
                cout << "Error in dwelling time from son. Branch of son is of length " << splitFromSon << " instead of " << AverageDwellingTimes[nodeState] << endl;
                return 1;
            }
            if (abs(splitToFather - AverageDwellingTimes[fatherState]) > 0.0001)
            {
                cout << "Error in dwelling time to father. Branch of son is of length " << splitToFather << " instead of " << AverageDwellingTimes[fatherState] << endl;
                return 1;
            }
        }

		// repeat the same tests for the analytic expected history
        cout << "**** Generating analytic expected mapping ****" << endl;
        TreeMapping analyticExpectedHistory = stocMapping->generateAnalyticExpectedMapping();
        checkIfMappingLegal(*stocMapping.get(), analyticExpectedHistory, *phyloTree.get(), *characterTreeLikelihood.get());
        cout << "**** The generated analytic mapping is legal ****" << endl;

        // // create another analytic expected mapping and make sure its equal to the former one (reconstruction is deterministic)
        // TreeMapping analyticExpectedHistory2 = stocMapping->generateAnalyticExpectedMapping();
        
		// // compare for each node in post-order traversal: name, state and distance to father
        // vector<Node*> hist1Nodes = dynamic_cast<TreeTemplate<Node>*>(analyticExpectedHistory)->getNodes();
        // vector<Node*> hist2Nodes = dynamic_cast<TreeTemplate<Node>*>(analyticExpectedHistory2)->getNodes();
        // if (hist1Nodes.size() != hist2Nodes.size())
        // {
        //     cerr << "Error! in repeated reconstruction of analytic expected history the number of nodes varies" << endl;
        //     return 1;
        // }
        // for (size_t n=0; n<hist1Nodes.size(); ++n)
        // {
        //     if (hist1Nodes[n]->getId() != hist2Nodes[n]->getId())
        //     {
        //         cerr << "Error! the nodes IDs in the two reconstrcued analyitc histories don't match for index " << n << endl;
        //         return 1;
        //     }
        //     if (StochasticMapping::getNodeState(hist1Nodes[n]) != StochasticMapping::getNodeState(hist2Nodes[n]))
        //     {
        //         cerr << "Error! the nodes states in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
        //         return 1;
        //     }
        //     if (hist1Nodes[n]->getId() != analyticExpectedHistory->getRootId())
        //     {
        //         if (abs(hist1Nodes[n]->getDistanceToFather() - hist2Nodes[n]->getDistanceToFather()) > 0.0001)
        //         {
        //             cerr << "Error! the branch lengths in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
        //             return 1;
        //         }
        //     }
        // } */
        
        delete process;
        delete alpha;

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}
