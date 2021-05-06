//
// File: CharacterSubstitutionModel.h
// Created by: Keren Halabi
// Created on: April 2021
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

#ifndef _CHARACTERSUBSTITUTIONMODEL_H_
#define _CHARACTERSUBSTITUTIONMODEL_H_

#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/IntegerFrequencySet.h"

#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

using namespace std;

namespace bpp
{

class CharacterSubstitutionModel :
  public AbstractSubstitutionModel
{
protected:
  shared_ptr<IntegerFrequencySet> freqSet_; 

public:
  /**
   * @brief Build a simple CharacterSubstitutionModel model, with original equilibrium frequencies.
   *
   * @param prefix A stirng of the model name
   * @param alpha An Integer alphabet
   */
  CharacterSubstitutionModel(const std::string& prefix, const IntegerAlphabet* alpha);

  /**
    * @brief Build a CharacterSubstitutionModel with special equilibrium frequencies.
    *
    * @param prefix A stirng of the model name
    * @param alpha An Integer alphabet.
    * @param freqSet A pointer toward a integer frequencies set, which will be owned by this instance.
    * @param initFreqs Tell if the frequency set should be initialized with the original mk values
    * (equal frequencies across all states). Otherwise, the values of the set will be used.
   */
  CharacterSubstitutionModel(const std::string& prefix, const IntegerAlphabet* alpha, shared_ptr<IntegerFrequencySet> freqSet, bool initFreqs = false);
  
  /** 
   * copy constructor
   */
  CharacterSubstitutionModel(const CharacterSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    freqSet_(dynamic_cast<IntegerFrequencySet *>(model.freqSet_->clone()))
  {}

  /** 
   * copy operator
   */
  CharacterSubstitutionModel& operator=(const CharacterSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractSubstitutionModel::operator=(model);
    freqSet_.reset(dynamic_cast<IntegerFrequencySet *>(model.freqSet_->clone()));
    return *this;
  }

  /**
   * destructor
   */
  virtual ~CharacterSubstitutionModel() {}

  /**
   * clone function
   */
  CharacterSubstitutionModel* clone() const { return new CharacterSubstitutionModel(*this); }

public:

  string getName() const // TO DO: make sure namespace doesn't already include the +F
  { 
    if (freqSet_->getNamespace().find("+F.")!=std::string::npos)
      return getNamespace()+"+F"; 
    else 
      return getNamespace(); 
  }
  
  size_t getNumberOfStates() const { return size_; }

  void fireParameterChanged(const ParameterList& parameters);

  void setFrequencySet(const IntegerFrequencySet& freqSet);
  
  const std::shared_ptr<FrequencySet> getFrequencySet() const { return freqSet_; }

  void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount = 0);

  /**
    * @brief Enables setting of search space of a parameter
    *
    * @param name Parameter name
    * @param lb Lower bound
    * @param ub Upper bound
   */
  void setParamBounds(const std::string& name, double lb, double ub);

protected:
  void updateMatrices();
};
} // end of namespace bpp.

#endif  // _CHARACTERSUBSTITUTIONMODEL_H_