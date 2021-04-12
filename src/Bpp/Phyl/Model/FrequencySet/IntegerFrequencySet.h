//
// File: IntegerFrequencySet.h
// Created by: Keren Halabi
// Created on: April 2021
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _INTEGERFREQUENCYSET_H_
#define _INTEGERFREQUENCYSET_H_

#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>
#include "FrequencySet.h"

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies for integer characters.
 */
  class IntegerFrequencySet :
    public virtual FrequencySet
  {
  public:
  
    IntegerFrequencySet* clone() const = 0;

    const IntegerAlphabet* getAlphabet() const = 0;
  };

/**
 * @brief Integer FrequencySet using k-1 independent parameters to
 * model the frequencies of k states.
 *
 * The parameters are called @f$ \theta_{i \in 1..k-1} @f$, and are
 * initialized so that all frequencies are equal to 1/k. The
 * parametrization depends on the method used. Default
 * method is 1 (ie global ratio).
 *
 * @see Simplex
 */
  class FullIntegerFrequencySet :
    public virtual IntegerFrequencySet,
    public FullFrequencySet
  {
  public:
    FullIntegerFrequencySet(const IntegerAlphabet* alphabet, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
      FullFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), allowNullFreqs, method, name) {}
    FullIntegerFrequencySet(const IntegerAlphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false, unsigned short method = 1, const std::string& name = "Full") :
      FullFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), initFreqs, allowNullFreqs, method, name) {}

    FullIntegerFrequencySet* clone() const { return new FullIntegerFrequencySet(*this); }

  public:
    const IntegerAlphabet* getAlphabet() const
    {
      return dynamic_cast<const IntegerAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };

/**
 * @brief FrequencySet useful for homogeneous and stationary models
 *
 * This set contains no parameter.
 */
  class FixedIntegerFrequencySet :
    public virtual IntegerFrequencySet,
    public FixedFrequencySet
  {
  public:
    FixedIntegerFrequencySet(const IntegerAlphabet* alphabet, const std::vector<double>& initFreqs, const std::string& name = "Fixed") :
      FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), initFreqs, name) {}

    /**
     * @brief Construction with uniform frequencies on the letters of
     * the alphabet.
     */
    FixedIntegerFrequencySet(const IntegerAlphabet* alphabet, const std::string& name = "Fixed") :
      FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), name) {}

    FixedIntegerFrequencySet* clone() const { return new FixedIntegerFrequencySet(*this); }

    const IntegerAlphabet* getAlphabet() const
    {
      return dynamic_cast<const IntegerAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };

  /**
   * @brief FrequencySet from file
   *
   * This set contains no parameter.
   */
  
  class UserIntegerFrequencySet :
    public virtual IntegerFrequencySet,
    public UserFrequencySet
  {
  public:
    UserIntegerFrequencySet(const IntegerAlphabet* alphabet, const std::string& path, size_t nCol=1) :
      UserFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), path, nCol) {}
    
    UserIntegerFrequencySet* clone() const { return new UserIntegerFrequencySet(*this); }

    const IntegerAlphabet* getAlphabet() const
    {
      return dynamic_cast<const IntegerAlphabet*>(AbstractFrequencySet::getAlphabet());
    }
  };


} // end of namespace bpp.

#endif // _INTEGERFREQUENCYSET_H_


