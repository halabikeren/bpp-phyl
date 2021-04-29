//
// File: SingleRateModel.h
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

#ifndef _SINGLERATEMODEL_H_
#define _SINGLERATEMODEL_H_

#include "CharacterSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/IntegerFrequencySet.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

namespace bpp
{

  /**
   * @brief The mk substitution model for characters.
   *
   * The original frequencies can be used, or alternatively a
   * parametrized version, corresponding to the so-called SingleRateModel+F
   * model. Eigen values and vectors are obtained numerically.
   * 
   * Reference:
   * Pagel, M. 1994. Detecting correlated evolution on phylogenies: 
   * A general method for the comparative analysis of discrete characters. 
   * Proc. R. Soc. Lond. B Biol. Sci. 255:37–45. The Royal Society.
   *
   */
  class SingleRateModel :
    public  CharacterSubstitutionModel
  {
  private:
    double globalRate_;

  public:
    /**
     * @brief Build a simple SingleRateModel model, with original equilibrium frequencies.
     *
     * @param alpha A character alphabet
     * @param fixedFreqs vector of fixed frequencies
     * @param rate Global substitution rate value
     */
    SingleRateModel(const IntegerAlphabet* alpha, const vector<double>& fixedFreqs, double rate = 1.);

    /**
     * @brief Build a SingleRateModel model with special equilibrium frequencies.
     *
     * @param alpha A character alphabet.
     * @param freqSet A pointer toward a character frequencies set, which will be owned by this instance.
     * @param rate Global substitution rate value
     * @param initFreqs Tell if the frequency set should be initialized with the original SingleRateModel values.
     * Otherwise, the values of the set will be used.

     */
    SingleRateModel(const IntegerAlphabet* alpha, std::shared_ptr<IntegerFrequencySet> freqSet, double rate = 1., bool initFreqs=false);

    /** 
     * copy constructor
     */
    SingleRateModel(const SingleRateModel& model) :
      AbstractParameterAliasable(model),
      CharacterSubstitutionModel(model),
      globalRate_(model.globalRate_)
    {}

    /** 
     * copy operator
     */
    SingleRateModel& operator=(const SingleRateModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      CharacterSubstitutionModel::operator=(model);
      globalRate_ = model.globalRate_;
      return *this;
    }

    /**
     * destructor
     */
    virtual ~SingleRateModel() {}

    /**
     * clone function
     */
    SingleRateModel* clone() const { return new SingleRateModel(*this); }

  public:
  
    void fireParameterChanged(const ParameterList& parameters)
    {
      globalRate_ = parameters.getParameterValue("rate");
      CharacterSubstitutionModel::fireParameterChanged(parameters);
    }
    
  protected:
    void updateMatrices();

  };

} //end of namespace bpp.

#endif	//_SINGLERATEMODEL_H_

