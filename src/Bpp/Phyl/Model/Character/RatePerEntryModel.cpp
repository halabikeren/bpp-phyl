//
// File: RatePerEntryModel.cpp
// Created by: Laurent Gueguen
// Created on: mardi 28 septembre 2010, à 14h 43
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


#include "RatePerEntryModel.h"
#include "../AbstractSubstitutionModel.h"

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

RatePerEntryModel::RatePerEntryModel(const IntegerAlphabet* alpha, const vector<double>& fixedFreqs, const vector<double>& entryRates):
  AbstractParameterAliasable("RatePerEntry."), // ask Tiana why this works
  CharacterSubstitutionModel("RatePerEntry", alpha, fixedFreqs),
  entryRates_(entryRates)
{
  for (size_t i=0; i<entryRates_.size(); ++i)
    addParameter_(new Parameter(getNamespace() + "entry_rate_"+TextTools::toString(i), entryRates_[i], std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  updateMatrices(); 
}

/******************************************************************************/

RatePerEntryModel::RatePerEntryModel(const IntegerAlphabet* alpha, std::shared_ptr<IntegerFrequencySet> freqSet, const vector<double>& entryRates, bool initFreqs):
  AbstractParameterAliasable("RatePerEntry+F."), // ask Tiana why this works
  CharacterSubstitutionModel("RatePerEntry", alpha, freqSet, initFreqs),
  entryRates_(entryRates)
{
  for (size_t i=0; i<entryRates_.size(); ++i)
    addParameter_(new Parameter(getNamespace() + "entry_rate_"+TextTools::toString(i), entryRates_[i], std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, false, false)));
  updateMatrices(); 
}

/******************************************************************************/

void RatePerEntryModel::updateMatrices()
{

  // update parameters
  for (size_t i=0; i<entryRates_.size(); ++i)
    entryRates_[i] = getParameterValue("entry_rate_" + TextTools::toString(i)); 
  freq_ = freqSet_->getFrequencies();

  vector<double> exitRates(size_);
  for (size_t i=0; i<size_; ++i)
  {
      double exitRate = 0;
      for (size_t j=0; j<size_; ++j)
      {
          if (j != i)
            exitRate += entryRates_[j];
      }
      exitRates[i] = exitRate;
  }

  // set the exachangability matrix 
  for (size_t i=0; i<size_; ++i)
  {
    for (size_t j=0; j<size_; ++j)
    {
      if (i == j)
      {
        exchangeability_(i,j) = - exitRates[i];
      }
      else
      {
        exchangeability_(i,j) = entryRates_[j];
      }
    }
  }

  MatrixTools::hadamardMult(exchangeability_, freq_, generator_, false); // Diagonal elements of the exchangeability matrix will be ignored.
  setDiagonal();
  isScalable_ = false;
  AbstractSubstitutionModel::updateMatrices();
}
