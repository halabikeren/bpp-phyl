//
// File: RatePerEntryModel.h
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

#ifndef _RATEPERENTRYMODEL_H_
#define _RATEPERENTRYMODEL_H_

#include "CharacterSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/IntegerFrequencySet.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

// From CoreLib:
#include <Bpp/Text/TextTools.h>

namespace bpp
{

  /**
   * @brief A substitution model for characters, allowing the rate of entry to each state to vary, 
   * leading to number of rate parameters which is equal to the number of states.
   */
  class RatePerEntryModel :
    public  CharacterSubstitutionModel
  {
  private:
    vector<double> entryRates_; // entryRates_[i] corresponds to the rate of entry to state i

  public:
    /**
     * @brief Build a simple RatePerEntryModel model, with original equilibrium frequencies.
     *
     * @param alpha A character alphabet
     */
    RatePerEntryModel(const IntegerAlphabet* alpha);

    /**
     * @brief Build a RatePerEntryModel model with special equilibrium frequencies.
     *
     * @param alpha A character alphabet.
     * @param freqSet A pointer toward a character frequencies set, which will be owned by this instance.
     * @param initFreqs Tell if the frequency set should be initialized with the original RatePerEntryModel values.
     * Otherwise, the values of the set will be used.
     */
    RatePerEntryModel(const IntegerAlphabet* alpha, std::shared_ptr<IntegerFrequencySet> freqSet, bool initFreqs=false);

    /** 
     * copy constructor
     */
    RatePerEntryModel(const RatePerEntryModel& model) :
      AbstractParameterAliasable(model),
      CharacterSubstitutionModel(model),
      entryRates_(model.entryRates_)
    {}

    /** 
     * copy operator
     */
    RatePerEntryModel& operator=(const RatePerEntryModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      CharacterSubstitutionModel::operator=(model);
      entryRates_ = model.entryRates_;
      return *this;
    }

    /**
     * destructor
     */
    virtual ~RatePerEntryModel() {}

    /**
     * clone function
     */
    RatePerEntryModel* clone() const { return new RatePerEntryModel(*this); }

  public:
  
    void fireParameterChanged(const ParameterList& parameters);

  protected:
    void updateMatrices();

  };

} //end of namespace bpp.

#endif	//_RATEPERENTRYMODEL_H_

