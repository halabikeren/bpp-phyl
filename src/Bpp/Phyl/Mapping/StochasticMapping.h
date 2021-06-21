//
// File: StochasticMapping.h
// Created by: Keren Halabi
// Created on: June 2018
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef ___STOCHASTIC_MAPPING_H
#define ___STOCHASTIC_MAPPING_H

#include "../NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../Simulation/MutationProcess.h"
#include "../NewLikelihood/DataFlow/DataFlowCWise.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

typedef vector<vector<vector<double> > > VVVDouble;
typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;

typedef map<unsigned int, bpp::MutationPath> TreeMapping; // maps branch index to its mutations path top to bottom (for ancestor to descendant)


/* class for reprenting the framework of Stochastic mapping
 *
 *   A StochasticMapping instance can be used to sample histories of
 *   state transitions along a tree, given a substitution model and the
 *   states at the tip taxa. For more information, see: Nielsen, Rasmus.
 *   "Mapping mutations on phylogenies." Systematic biology 51.5 (2002):
 *   729-739.‏
 */

namespace bpp
{
    class StochasticMapping
    {
    protected:

        /*
        * @brief The tree likelihood instance is used for computing the
        * the conditional sampling probabilities of the ancestral states as
        * well as the root assignment probabilities.
        *
        */
    
        shared_ptr<LikelihoodCalculationSingleProcess> likelihood_;  
        shared_ptr<const PhyloTree> tree_;

        /*
        * @briefvector that holds the fractional probabilities per state
        * per node in the tree, based on which the conditional and
        * posterior probabilities are computed
        *
        */

        VVDouble fractionalProbabilities_;               // vector that holds the fractional probabilities of states assignments to the root node
        VVVDouble ConditionalProbabilities_;             // vector that holds the conditional probabilities of states assignment probabilities of the nodes in the tree (node*father_states*son_states)
        unsigned int numOfMappings_;                     // the number of stochastic mappings to generate

    public:
        /* constructors and destructors */

        explicit StochasticMapping(shared_ptr<LikelihoodCalculationSingleProcess> drl, unsigned int numOfMappings = 10000); // it is a good general practice to use "explicit" keyword on constructors with a single argument: https://stackoverflow.com/questions/121162/what-does-the-explicit-keyword-mean

        ~StochasticMapping() {}

        StochasticMapping(const StochasticMapping& sm) : 
        likelihood_(sm.likelihood_),
        tree_(sm.tree_),
        fractionalProbabilities_(sm.fractionalProbabilities_),
        ConditionalProbabilities_(sm.ConditionalProbabilities_)
        {}

        /**
         * @brief cloning function used by the copy constructor of JointPhyloLikelihood
         *
         */
        StochasticMapping* clone() const { return new StochasticMapping(*this); }

        /**
         * @brief Assignment operator
         */
        StochasticMapping& operator=(const StochasticMapping& sm)
        {
            likelihood_ = sm.likelihood_;
            tree_ = sm.tree_;
            fractionalProbabilities_ = sm.fractionalProbabilities_;
            ConditionalProbabilities_ = sm.ConditionalProbabilities_;
            numOfMappings_ = sm.numOfMappings_;
            return *this;
        }

        shared_ptr<const PhyloTree> getBaseTree() { return tree_; }

        /*
        *
        * @brief generates a stochastic mappings based on the sampling
        * 
        */
        vector<TreeMapping> generateStochasticMappings();

        /*
        * @brief compute the distance between two mappings as the complement of their overlap precentage
        * parameters
        * 
        * @param  mapping1  first mapping
        * @param  mapping2  second mapping
        * 
        */        
        double getDistance(const TreeMapping& mapping1, const TreeMapping& mapping2);

        /*
        * @brief compute the distance between two branch mappings in the form of mutation paths as the complement of their overlap precentage
        * parameters
        * 
        * @param  path1  first branch mapping
        * @param  path2  second branch mapping
        * 
        */    
        double getPathsDistance(const MutationPath& path1, const MutationPath& path2);

        /*
        * @brief clusters stochastic mappings using k-means based on their distannce
        * parameters
        * 
        * @param mappings   the mappings to cluster
        * @param k          number of clusters to generate
        * 
        */        
        vector<TreeMapping> kMeansClustering(vector<TreeMapping>& mappings, unsigned int k, unsigned int epochs = 20);

        /*
        * @brief generates a stochastic mappings based on the sampling and then clusters them based on similarity
        * parameters
        * 
        * @param            number of clusters to generate
        * 
        */
        vector<TreeMapping> generatedClusteredMappings(unsigned int numOfClusters = 10);

        /**
         *@brief Creates a single expected (i.e, average) history based on
        * a given set of mappings steps. Correspond to Nielsen, Rasmus.
        * "Mapping mutations on phylogenies." Systematic biology 51.5
        * (2002): 729-739.‏
        *
        * the function assumes that there is only one site to simulate history
        *
        * @param mappings         the mappings to average over
        **/
        TreeMapping generateExpectedMapping(const vector<TreeMapping>& mappings);

        /**
         *@brief returns the analytically estimated dwelling times under each model state in each branch
         * the function assumes that there is only one site to simulate history for
        */
        VVDouble getAnalyticDwellingTimes();

        /**
         *@brief Creates a single expected (i.e, average) history based the rewards provided by te algorithm of Minin and Suchard (2008)
        * the function assumes that there is only one site to simulate history for
        */
        TreeMapping generateAnalyticExpectedMapping();

        /* extracts the state of a node in a mapping
        * @param node              The node to get the state of
        * @param mapping           The mapping from which the node state should be extracted
        * @return                  Node state is int
        */
        unsigned int getNodeState(shared_ptr<PhyloNode> node, const TreeMapping& mapping);

        /* compute the ancestral frequenceis of character states of all the nodes based on the mappings
        * @param                     A vector of the posterior probabilities probabilities to fill in (node**state combinaion in each entry)
        * @param                     A vector of mappings to base the frequencies on
        */
        void computeStatesFrequencies(VVDouble& ancestralStatesFreuquencies, const vector<TreeMapping>& mappings);

    private:

        /* returns the possible model states assigned to a leaf based on its character state
        * @param node - the leaf of interest
        */
        vector<unsigned int> getLeafModelStates(PhyloNode& node);

        /* returns the the character state of a leaf based on the character data
        * @param leafName - the name of the leaf of interest
        */
        unsigned int getLeafCharacterState(string leafName);

        /* sets the state of a node in a mapping
        * @param node               The node to get the state of
        * @param state              The state that needs to be assigned to the node
        * @param mapping            The mapping to apply the node state in
        */
        void setNodeState(shared_ptr<PhyloNode> node, unsigned int state, TreeMapping& mapping);

        /* return the posterior assignment probabilities of states to all the nodes, based on the computed fractional probabilities
        */
        VVDouble getPosteriorProbabilities();

        /* compute the fractional probabilities of all the nodes assignements
        * @param                     A vector of the fractional probabilities probabilities to fill in (node**state combinaion in each entry)
        */
        void computeFractionals();

        /* compute the conditional probabilities of all the nodes assignments
        * @param rootProbabilities  The root frequencies
        */
        void ComputeConditionals();

        /* auxiliary function that samples a state based on a given discrete distribution
        * @param distibution       The distribution to sample states based on
        */
        unsigned int sampleState(const VDouble& distibution); // k: best by ref

        /* samples ancestral states based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
        * @param mapping               The tree whose nodes names should be updated according to their assigned states.
        */
        void sampleAncestrals(TreeMapping& mapping);

        /* set ancestral states in the expected history based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
        * @param mapping                   The expected mapping instance whose nodes names should be updated according to their assigned states.
        * @param posteriorProbabilities    Vector of posterior assignment proabilities to inner node to decide on assignments
        */
        void setExpectedAncestrals(TreeMapping& expectedMapping, VVDouble& posteriorProbabilities);

        /* simulates mutations on phylogeny based the sampled ancestrals, tips data, and the simulation parameters
        * @param mapping               The tree that should be edited according to the sampled history
        */
        void sampleMutationsGivenAncestrals(TreeMapping& mapping);

        /* adds a branch mapping to the mapping in a tree format by repeatedly braking branches and adding internal nodes with single children
        * @param mapping               The mapping that contains node "son"
        * @param son                   The node at the bottom of the branch
        * @param sampledMapping        The successfully sampled MutationPath of transitions in the branch
        */
        void updateBranchMapping(TreeMapping& mapping, shared_ptr<PhyloNode> son, const MutationPath& sampledMapping);
        
        /* sample mutations based on the stochastic mapping parameters, the source and destination state, and the branch length, and updates the simulated history along the branch in the input tree, on the fly
        * @param mapping               The mapping that contains node "son"
        * @param son                   Node of interest
        * @param maxIterNum            Maximal number of imulation trials
        */
        void sampleMutationsGivenAncestralsPerBranch(TreeMapping& mapping, shared_ptr<PhyloNode> son, unsigned int maxIterNum = 10000);

        /* converts a vector of dwelling times to a mutation path and then updates the bracnh stemming from the given node */
        /* @param node                      The node at the bottom of the branch
        * @param dwellingTimes             A vector of dwelling times where the value at each entry i corresponds to the dwelling time under the i'th state
        *                                  Note that this function generates a new tree instance, that must be deleted by the calling function.
        * @param posteriorProbabilities    Posterior probaibitlies to divide the time spent in a shared state between father and son into two transitions
        */
        void updateBranchByDwellingTimes(TreeMapping& mapping, shared_ptr<PhyloNode> node, VDouble& dwellingTimes, VVDouble& ancestralStatesFrequencies);
    };
}

#endif// ___STOCHASTIC_MAPPING_H