//
// File: mkModel.cpp
// Created by: Keren Halabi
// Created on: 2021
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "mkModel.h"

// From the STL:
#include <cmath>

using namespace bpp;


//From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

//From bpp-seq:
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

using namespace std;

/******************************************************************************/

mkModel::mkModel(const IntegerAlphabet* alpha):
  AbstractParameterAliasable("mk."),
  AbstractReversibleSubstitutionModel(alpha, shared_ptr<const StateMap>(new CanonicalStateMap(alpha, false)), "mk."),
  rates_(size_, size_), // the user can select to constrain the rates to be similar across all pairs of states to obtain the actual mk model
  freqSet_(0)
{
  freqSet_.reset(new FixedIntegerFrequencySet(alpha, freq_));
  for (size_t i=0; i<rates_.getNumberOfRows(); ++i)
  {
    for (size_t j=0; j<rates_.getNumberOfRows(); ++j)
    {
      rates_(i, j) = 1.;
      addParameter_(new Parameter("mk.rate_" + TextTools::toString(i) + "_" + TextTools::toString(j), rates_(i, j), make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
    }
  }
  updateMatrices();
}

/******************************************************************************/

mkModel::mkModel(const IntegerAlphabet* alpha, shared_ptr<IntegerFrequencySet> freqSet, bool initFreqs) :
  AbstractParameterAliasable("mk+F."),
  AbstractReversibleSubstitutionModel(alpha, shared_ptr<const StateMap>(new CanonicalStateMap(alpha, false)), "mk+F."),
  rates_(size_, size_),
  freqSet_(freqSet)
{
  freqSet_->setNamespace("mk+F."+freqSet_->getNamespace());
  if (initFreqs) freqSet_->setFrequencies(freq_);
  else freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  for (size_t i=0; i<rates_.getNumberOfRows(); ++i)
  {
    for (size_t j=0; j<rates_.getNumberOfRows(); ++j)
    {
      rates_(i, j) = 1.;
      addParameter_(new Parameter("mk+F.rate_" + TextTools::toString(i) + "_" + TextTools::toString(j), rates_(i, j), make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
    }
  }
  updateMatrices();  
}

/******************************************************************************/


void mkModel::updateMatrices()
{

  // update parameters
  for (size_t i=0; i<rates_.getNumberOfRows(); ++i)
  {
    for (size_t j=0; j<rates_.getNumberOfColumns(); ++j)
    {
      rates_(i, j) = getParameterValue("rate_" + TextTools::toString(i) + "_" + TextTools::toString(j));
    }
  }
  freq_ = freqSet_->getFrequencies();

  // set the exachangability matrix 
  for (size_t i=0; i<size_; ++i)
  {
    for (size_t j=0; j<size_; ++j)
    {
      if (i == j)
      {
        exchangeability_(i,j) = -1 * rates_(i, j);
      }
      else
      {
        exchangeability_(i,j) = rates_(i, j);
      }
    }
  }

  isScalable_ = false;
  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void mkModel::fireParameterChanged(const ParameterList& parameters)
{
  for (size_t i=0; i<rates_.getNumberOfRows(); ++i)
  {
    for (size_t j=0; j<rates_.getNumberOfColumns(); ++j)
    {
      rates_(i, j)= parameters.getParameterValue(getName() + ".rate+" + TextTools::toString(i) + "_" +  TextTools::toString(j));
    }
  }
  freqSet_->matchParametersValues(parameters);
  freq_ = freqSet_->getFrequencies();
  AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
}

/******************************************************************************/

void mkModel::setFrequencySet(const IntegerFrequencySet& freqSet)
{
  freqSet_ = shared_ptr<IntegerFrequencySet>(dynamic_cast<IntegerFrequencySet*>(freqSet.clone()));
  resetParameters_();
  addParameters_(freqSet_->getParameters());
}

/******************************************************************************/

void mkModel::setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
{
  map<int, double> counts;
  SequenceContainerTools::getFrequencies(data, counts, pseudoCount);
  for (auto i : counts)
    freq_[(size_t)i.first] = i.second;
  
  freqSet_->setFrequencies(freq_);
  //Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/

void mkModel::setRateBounds(int srcState, int destState, double lb, double ub)
{
  shared_ptr<IntervalConstraint> bounds(new IntervalConstraint(lb, ub, true, true)); 
  getParameter_("rate_" + TextTools::toString(srcState) + "_" + TextTools::toString(destState)).setConstraint(bounds);
}
