/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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
encoutaged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/


#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/Codon/RELAX.h>
#include <Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/MixtureSequenceEvolution.h>


#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MixtureOfAlignedPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/FormulaOfPhyloLikelihood.h>

#include <iostream>

using namespace bpp;
using namespace std;


int main() 
{
    try
    {
        // process tree
        Newick reader;
        unique_ptr<PhyloTree> tree(reader.parenthesisToPhyloTree("(((A:0.01, B:0.01):0.02,C:0.03):0.01,D:0.04);"));
        auto parTree = make_shared<ParametrizablePhyloTree>(*tree);

        // process sequence data
        const CodonAlphabet* alphabet = &AlphabetTools::DNA_CODON_ALPHABET;
        VectorSiteContainer sites(alphabet);
        sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTCT", alphabet));
        sites.addSequence(BasicSequence("B", "AACTGGATCTGCATGTCT", alphabet));
        sites.addSequence(BasicSequence("C", "ATCTGGACGTGCACGTGT", alphabet));
        sites.addSequence(BasicSequence("D", "CAACGGGAGTGCGCCTAT", alphabet));

        //// create a branch-site model with two site-categories and three site categories
        auto gc = make_shared<StandardGeneticCode>(AlphabetTools::DNA_ALPHABET);
        auto rdist = make_shared<ConstantRateDistribution>();
        shared_ptr<FrequencySet> bgRootFreqs = CodonFrequencySet::getFrequencySetForCodons(CodonFrequencySet::F3X4, gc.get());
        auto bgModel = make_shared<RELAX>(gc.get(), bgRootFreqs);
        shared_ptr<FrequencySet> fgRootFreqs = CodonFrequencySet::getFrequencySetForCodons(CodonFrequencySet::F3X4, gc.get());
        auto fgModel = make_shared<RELAX>(gc.get(), fgRootFreqs);

        NonHomogeneousSubstitutionProcess* subProc=new NonHomogeneousSubstitutionProcess(rdist->clone(), parTree->clone());
        Vuint vP1m1{0}; // branches of bgModel
        Vuint vP1m2{1, 2, 3, 4, 5}; // branches of fgModel
        subProc->addModel(std::shared_ptr<RELAX>(bgModel->clone()),vP1m1);
        subProc->addModel(std::shared_ptr<RELAX>(fgModel->clone()),vP1m2);

        auto modelColl=make_shared<SubstitutionProcessCollection>();
        modelColl->addModel(bgModel, 1);
        modelColl->addModel(fgModel, 2);
        modelColl->addFrequencies(bgRootFreqs, 1);
        modelColl->addFrequencies(fgRootFreqs, 2);
        modelColl->addDistribution(rdist, 1);
        modelColl->addTree(parTree, 1);

        map<size_t, Vuint> mModBr;
        mModBr[1]=vP1m1;
        mModBr[2]=vP1m2;
        modelColl->addSubstitutionProcess(1, mModBr, 1, 1, 1);

        vector<size_t> vp(1);
        vp[0]=1;
        MixtureSequenceEvolution mse(modelColl.get(), vp);

        // alias parameters
        map<string,string> sharedParamsMap = {
            {"RELAX.kappa_1", "RELAX.kappa_2"},
            {"RELAX.p_1", "RELAX.p_2"},
            {"RELAX.omega1_1", "RELAX.omega1_2"},
            {"RELAX.omega2_1", "RELAX.omega2_2"},
            {"RELAX.theta1_1", "RELAX.theta1_2"},
            {"RELAX.theta2_1", "RELAX.theta2_2"},
            {"RELAX.1_Full.theta_1", "RELAX.1_Full.theta_2"}, // why do I need F3X4 paramters prefices to be both "RELAX" and "YN98"?
            {"RELAX.1_Full.theta1_1", "RELAX.1_Full.theta1_2"},
            {"RELAX.1_Full.theta2_1", "RELAX.1_Full.theta2_2"},
            {"RELAX.2_Full.theta_1", "RELAX.2_Full.theta_2"},
            {"RELAX.2_Full.theta1_1", "RELAX.2_Full.theta1_2"},
            {"RELAX.2_Full.theta2_1", "RELAX.2_Full.theta2_2"},
            {"RELAX.3_Full.theta_1", "RELAX.3_Full.theta_2"},
            {"RELAX.3_Full.theta1_1", "RELAX.3_Full.theta1_2"},
            {"RELAX.3_Full.theta2_1", "RELAX.3_Full.theta2_2"},
            {"YN98.1_Full.theta_1", "YN98.1_Full.theta_2"},
            {"YN98.1_Full.theta1_1", "YN98.1_Full.theta1_2"},
            {"YN98.1_Full.theta2_1", "YN98.1_Full.theta2_2"},
            {"YN98.2_Full.theta_1", "YN98.2_Full.theta_2"},
            {"YN98.2_Full.theta1_1", "YN98.2_Full.theta1_2"},
            {"YN98.2_Full.theta2_1", "YN98.2_Full.theta2_2"},
            {"YN98.3_Full.theta_1", "YN98.3_Full.theta_2"},
            {"YN98.3_Full.theta1_1", "YN98.3_Full.theta1_2"},
            {"YN98.3_Full.theta2_1", "YN98.3_Full.theta2_2"},
        };
        modelColl->aliasParameters(sharedParamsMap, true);
        ParameterList aliasedParams = modelColl->getIndependentParameters();
        cout << "aliased parameters: " << endl;
        for (size_t p=0; p<aliasedParams.size(); ++p)
        {
            cout << aliasedParams[p].getName() << endl;
        }

        // create a likelihood nodes for each branch category
        Context context;
        auto pc(std::make_shared<PhyloLikelihoodContainer>(context, *modelColl));
        SubstitutionProcess* sP=subProc->clone();
        auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, *sP);
        pc->addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(context, lik));
        auto collNodes = pc->getCollectionNodes();

        // create a mixture likelihood function
        MixtureProcessPhyloLikelihood mlc(*sites.clone(), mse, *collNodes);

        // print model parameter values
        cout << "\n\ninitial parameter values" << endl;
        mlc.getParameters().printParameters(std::cout);

        // report likelihood
        cout << "\n\ninitial likelihood: " << -mlc.getValue() << endl;

        
        // optimize likelihood function
        OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
        OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));
        OptimizationTools::optimizeNumericalParameters2(&mlc, mlc.getParameters(), 0, 0.01, 100, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

        // print optimized parameter values - TO DO: aliasing didn't work! 
        cout << "\n\nparameter values following optimization" << endl;
        mlc.getParameters().printParameters(std::cout);

        // report likelihood 
        cout << "\n\nlikelihood following optimization: " << -mlc.getValue() << endl;
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
  return 0;
}
