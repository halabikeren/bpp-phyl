//
// File: ProcessComputationTree.h
// Created by: Laurent Guéguen
// Created on: Sat Dec 30 12:48 2006
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _PROCESS_COMPUTATION_TREE_H_
#define _PROCESS_COMPUTATION_TREE_H_

#include "../Tree/PhyloTree.h"
#include "../Tree/PhyloBranchParam.h"

#include "ParametrizablePhyloTree.h"
#include "SubstitutionProcess.h"

namespace bpp
{
/**
 * @brief Tree Organization of Computing Nodes
 *
 * Stores computation tools for all nodes, for all classes.
 *
 * This object has the parameters of the Tree and the rates
 * distribution, since it manages the SpeciationComputingNodes.
 *
 */

  // Node specific DataFlow objects
  class ProcessComputationNode:
    public PhyloNode
  {
    /*
     * @brief the index of the species in the phyloTree matching this node.
     *
     */
      
    const uint speciesIndex_;

  public:

    /*
     * @brief Build from a node in the phylo tree, with a specific
     * speciesIndex (because indexes are not the same as in the
     * ParametrizablePhyloTree.
     *
     */
      
    ProcessComputationNode(const PhyloNode& node, uint speciesIndex) :
      PhyloNode(node),
      speciesIndex_(speciesIndex) {}

    ProcessComputationNode(const ProcessComputationNode& node) :
      PhyloNode(node),
      speciesIndex_(node.speciesIndex_) {}

    uint getSpeciesIndex() const
    {
      return speciesIndex_;
    }

  };

  using ProcessComputationNodeRef=std::shared_ptr<ProcessComputationNode>;
  
  // Class for the edges
  class ProcessComputationEdge
  {
  private:
    /*
     * Model carried by the branch
     *
     */
     
    const TransitionModel* model_;

    /*
     * @brief numbers of the submodels used, if any.
     *
     * In practice, only vectors of size <=1 are implemented (ie
     * individual submodels), to be fixed later.
     */
    
    std::vector<uint> vSubNb_;

    /*
     * @brief use the probability associated to the edge in case of
     * mixture. If true, the probability is used and not the
     * transition probabilities.
     *
     */

    bool useProb_;
    
    /*
     * @brief the index of the species in the phyloTree matching this node.
     *
     */
      
    const uint speciesIndex_;


  public:
    ProcessComputationEdge(const TransitionModel* model, uint speciesIndex, bool useProb=false, const std::vector<uint>& vNb=Vuint(0)) :
      model_(model),
      vSubNb_(vNb),
      useProb_(useProb),
      speciesIndex_(speciesIndex)
    {
    }

    const TransitionModel* getModel() const
    {
      return model_;
    }

    uint getSpeciesIndex() const
    {
      return speciesIndex_;
    }

    const std::vector<uint>& subModelNumbers() const
    {
      return vSubNb_;
    }

    bool useProb() const
    {
      return useProb_;
    }
    
  };

  using BaseTree = AssociationTreeGlobalGraphObserver<ProcessComputationNode, ProcessComputationEdge>;

  class ProcessComputationTree :
    public BaseTree
  {
  private:
    
    const SubstitutionProcess& process_;

    /*
     * Build rest of the tree under given father (event on father will
     * be set in this method)
     *
     */
    
    void _build_following_path(std::shared_ptr<ProcessComputationNode> father, const ModelPath& path);

    void _build_following_scenario(std::shared_ptr<ProcessComputationNode> father, const ModelScenario& scenario,  std::map<std::shared_ptr<MixedTransitionModel>, uint>& mMrca);


  public:
    /*
     * @brief construction of a ProcessComputationTree from a SubstitutionProcess
     *
     * @param process the SubstitutionProcess
     *
     */
     
    ProcessComputationTree(const SubstitutionProcess& process);
     

    // ProcessComputationTree(const ProcessComputationTree& tree) :
    //   BaseTree(tree)
 
    // ProcessComputationTree& operator=(const ProcessComputationTree& tree)
    // {
    //   BaseTree::operator=(tree);
    //   return *this;
    // }
    // ProcessComputationTree* clone() const { return new ProcessComputationTree(*this);}
      

  //   /*
  //    * @brief construction of a complete ProcessComputationTree.
  //    *
  //    * @param pSubMod a  pointer of TransitionModel.
  //    * @param vBr a vector of attribution of the model on the
  //    * branches of the tree that need a model.
  //    *
  //    */
     
  //   void addModel(const TransitionModel* pSubMod, std::vector<unsigned int> vBr);

  //   /*
  //    * @brief construction of an homogeneous ProcessComputationTree.
  //    *
  //    * @param pSubMod a  pointer of TransitionModel.
  //    *
  //    */
     
  //   void addModel(const TransitionModel* pSubMod);
    
  //   size_t getNumberOfClasses() const { return vTree_.size();}

  // private:
    
  //   void clearAllModels_();

  // public:

  //   /*
  //    * @brief Checks if every SpeciationComputingNode has a Model
  //    *
  //    */
    
  //   void checkModelOnEachNode();

  //   /*
  //    * @brief operator to get numbered TreeTemplate<SpeciationComputingNode>*
  //    *
  //    */
    
  //   std::shared_ptr<CompTree> operator[](size_t ntree) { return vTree_[ntree];}

  //   const std::shared_ptr<CompTree> operator[](size_t ntree) const { return vTree_[ntree];}

  //   /*
  //    *@brief update Distribution parameters and says to the
  //    * SpeciationComputingNodes to be ready to update if the Branch lengths are
  //    * changed.
  //    *
  //    */
    
  //   void fireParameterChanged(const ParameterList& pl);
        
  //   /*
  //    * @brief Says to specific nodes to be ready for update
  //    *
  //    * If flag = true (default), node has to be updated (false otherwise).
  //    */

  //   void update(std::vector<unsigned int>& vId, bool flag = true);

  //   void update(unsigned int id, bool flag = true);

  //   /*
  //    * @brief Says to all nodes to be ready for update
  //    *
  //    */
    
  //   void updateAll();

  //   /*
  //    * @brief Returns the list of the updated Nodes ids, excluding the root.
  //    *
  //    */
    
  //   Vuint toBeUpdatedNodes() const;

  // private:

  //   /*
  //    * @brief Initializes the Computing tree in homogeneous constructors.
  //    *
  //    */
    
  //   void init_();
    
  };
  
} //end of namespace bpp.

#endif //_COMPUTINGTREE_H_
