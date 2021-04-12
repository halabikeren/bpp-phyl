//
// File: mkModel.h
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

#ifndef _MKMODEL_H_
#define _MKMODEL_H_

#include "AbstractSubstitutionModel.h"
#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>
#include "./FrequencySet/IntegerFrequencySet.h"

using namespace std;

namespace bpp
{

class mkModel :
  public AbstractReversibleSubstitutionModel
{
private:
  mutable RowMatrix<double> rates_;
  shared_ptr<IntegerFrequencySet> freqSet_;

public:
  /**
   * @brief Build a simple mkModel model, with original equilibrium frequencies.
   *
   * @param alpha An Integer alphabet.
   */
  mkModel(const IntegerAlphabet* alpha);

  /**
    * @brief Build a mkModel with special equilibrium frequencies.
    *
    * @param alpha An Integer alphabet.
    * @param freqSet A pointer toward a integer frequencies set, which will be owned by this instance.
    * @param initFreqs Tell if the frequency set should be initialized with the original mk values
    * (equal frequencies across all states). Otherwise, the values of the set will be used.
   */
  mkModel(const IntegerAlphabet* alpha, shared_ptr<IntegerFrequencySet> freqSet, bool initFreqs = false);
  
  /** 
   * copy constructor
   */
  mkModel(const mkModel& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleSubstitutionModel(model),
    rates_(model.rates_),
    freqSet_(dynamic_cast<IntegerFrequencySet *>(model.freqSet_->clone()))
  {}

  /** 
   * copy operator
   */
  mkModel& operator=(const mkModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleSubstitutionModel::operator=(model);
    freqSet_.reset(dynamic_cast<IntegerFrequencySet *>(model.freqSet_->clone()));
    rates_ = model.rates_;
    return *this;
  }

  /**
   * destructor
   */
  virtual ~mkModel() {}

  /**
   * clone function
   */
  mkModel* clone() const { return new mkModel(*this); }

  
public:

  string getName() const 
  { 
    if (freqSet_->getNamespace().find("mk+F.")!=std::string::npos)
      return "mk+F"; 
    else 
      return "mk"; 
  }
  
  size_t getNumberOfStates() const { return size_; }

  void fireParameterChanged(const ParameterList& parameters);

  void setFrequencySet(const IntegerFrequencySet& freqSet);
  
  const std::shared_ptr<FrequencySet> getFrequencySet() const { return freqSet_; }

  void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0);

  /**
   * @brief enables customization of rate search space during optimization
   * 
   * @param srcState the source state
   * @param destState the destination state
   * @param lb a lower bound value
   * @param ub an upper bound value
   */
  void setRateBounds(int srcState, int destState, double lb, double ub);

protected:
  void updateMatrices();
};
} // end of namespace bpp.

#endif  // _MKMODEL_H_