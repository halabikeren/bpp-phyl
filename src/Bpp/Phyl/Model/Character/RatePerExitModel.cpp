//
// File: RatePerExitModel.cpp
// Created by: Keren Halabi
// Created on: April 2021
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


#include "RatePerExitModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

RatePerExitModel::RatePerExitModel(const IntegerAlphabet* alpha):
  AbstractParameterAliasable("RatePerExit."),
  CharacterSubstitutionModel("RatePerExit", alpha),
  exitRates_(size_)
{
  for (size_t i=0; i<exitRates_.size(); ++i)
  {
    exitRates_[i] = 1.;
    addParameter_(new Parameter(getNamespace() + "exit_rate_"+TextTools::toString(i), exitRates_[i], std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  }
  updateMatrices(); 
}

/******************************************************************************/

RatePerExitModel::RatePerExitModel(const IntegerAlphabet* alpha, std::shared_ptr<IntegerFrequencySet> freqSet, bool initFreqs):
  AbstractParameterAliasable("RatePerExit+F."), // ask Tiana why this works
  CharacterSubstitutionModel("RatePerExit", alpha, freqSet, initFreqs),
  exitRates_(size_)
{
  for (size_t i=0; i<exitRates_.size(); ++i)
  {
    exitRates_[i] = 1.;
    addParameter_(new Parameter(getNamespace() + "exit_rate_"+TextTools::toString(i), exitRates_[i], std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  }
  updateMatrices(); 
}

/******************************************************************************/

void RatePerExitModel::updateMatrices()
{

  // update parameters
  for (size_t i=0; i<exitRates_.size(); ++i)
    exitRates_[i] = getParameterValue("exit_rate_" + TextTools::toString(i)); 
  freq_ = freqSet_->getFrequencies();

  // set the exachangability matrix 
  for (size_t i=0; i<size_; ++i)
  {
    for (size_t j=0; j<size_; ++j)
    {
      exchangeability_(i,j) = ((i == j) ?  - exitRates_[i]: (exitRates_[i]/static_cast<double>(size_-1)));
    }
  }

  MatrixTools::hadamardMult(exchangeability_, freq_, generator_, false); // Diagonal elements of the exchangeability matrix will be ignored.
  setDiagonal();
  isScalable_ = false;
  AbstractSubstitutionModel::updateMatrices();
}


/******************************************************************************/

void RatePerExitModel::fireParameterChanged(const ParameterList& parameters)
{
  for (size_t i=0; i<exitRates_.size(); ++i)
  {
    if (parameters.hasParameter(getNamespace() + "exit_rate_" + TextTools::toString(i)))
      exitRates_[i] = parameters.getParameterValue(getNamespace() + "exit_rate_" + TextTools::toString(i));
  }
  CharacterSubstitutionModel::fireParameterChanged(parameters);
}