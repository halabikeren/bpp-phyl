//
// File: FrequenciesSet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
   Copyright or (c) or Copr. CNRS, (November 16, 2004)

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

#ifndef _FREQUENCIESSET_H_
#define _FREQUENCIESSET_H_

// From NumCalc:
#include <NumCalc/Parametrizable.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AbstractParametrizable.h>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/NucleicAlphabet.h>
#include <Seq/ProteicAlphabet.h>
#include <Seq/WordAlphabet.h>
#include <Seq/CodonAlphabet.h>
#include <Seq/GeneticCode.h>

namespace bpp
{
/**
 * @brief Parametrize a set of state frequencies.
 */
class FrequenciesSet :
  public virtual Parametrizable
{
public:
#ifndef NO_VIRTUAL_COV
  FrequenciesSet* clone() const = 0;
#endif

public:
  /**
   * @return The alphabet associated to this set.
   */
  virtual const Alphabet* getAlphabet() const = 0;

  /**
   * @return The frequencies values of the set.
   */
  virtual const std::vector<double>& getFrequencies() const = 0;

  /**
   * @brief Set the parameters in order to match a given set of frequencies.
   *
   * @param frequencies The set of frequencies to match.
   * @throw DimensionException If the number of frequencies does not match the size of the alphabet.
   * @throw Exception If the frequencies do not sum to 1.
   */
  virtual void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception) = 0;

  /**
   * @brief Set the Frequencies from the one of the map which keys
   *  match with a letter of the Alphabet.
   *  The frequencies are normalized so that the matching values sum 1.
   *
   * @param frequencies The set of frequencies to match.
   */
  virtual void setFrequenciesFromMap(const std::map<int, double>& frequencies) = 0;

  virtual std::string getName() const = 0;

  /**
   * @return The number of frequencies in the set. In most cases this will correspond to the size of the alphabet,
   * but it needs not be.
   */
  virtual unsigned int getNumberOfFrequencies() const = 0;

public:
  /**
   * @brief A helper function that provide frequencies set for codon models
   * according to PAML option.
   *
   * @param option A code describing the option, one of F61, F1X4 or F3X4.
   * @param gc The genetic code to use.
   */
  static FrequenciesSet* getFrequenciesSetForCodons(short option, const GeneticCode& gc);
  static const short F0;
  static const short F1X4;
  static const short F3X4;
  static const short F61;
};

/**
 * @brief Parametrize a set of state frequencies for nucleotides.
 */
class NucleotideFrequenciesSet :
  public virtual FrequenciesSet
{
public:
#ifndef NO_VIRTUAL_COV
  NucleotideFrequenciesSet* clone() const = 0;

  const NucleicAlphabet* getAlphabet() const = 0;
#endif
};

/**
 * @brief Parametrize a set of state frequencies for proteins.
 */
class ProteinFrequenciesSet :
  public virtual FrequenciesSet
{
public:
#ifndef NO_VIRTUAL_COV
  ProteinFrequenciesSet* clone() const = 0;

  const ProteicAlphabet* getAlphabet() const = 0;
#endif
};

/**
 * @brief Parametrize a set of state frequencies for codons.
 */
class CodonFrequenciesSet :
  public virtual FrequenciesSet
{
public:
#ifndef NO_VIRTUAL_COV
  CodonFrequenciesSet* clone() const = 0;

  const CodonAlphabet* getAlphabet() const = 0;
#endif
};


/**
 * @brief Basic implementation of the FrequenciesSet interface.
 */

class AbstractFrequenciesSet :
  public virtual FrequenciesSet,
  public AbstractParametrizable
{
private:
  const Alphabet* alphabet_;
  std::vector<double> freq_;

public:
  AbstractFrequenciesSet(unsigned int n, const Alphabet* alphabet, const std::string& prefix) :
    AbstractParametrizable(prefix),
    alphabet_(alphabet),
    freq_(n) {}

#ifndef NO_VIRTUAL_COV
  AbstractFrequenciesSet*
#else
  Clonable*
#endif
  clone() const = 0;

  AbstractFrequenciesSet(const AbstractFrequenciesSet& af) :
    AbstractParametrizable(af),
    alphabet_(af.alphabet_),
    freq_(af.freq_) {}

  AbstractFrequenciesSet & operator=(const AbstractFrequenciesSet& af)
  {
    AbstractParametrizable::operator=(af);
    alphabet_ = af.alphabet_;
    freq_ = af.freq_;
    return *this;
  }

public:
  const Alphabet* getAlphabet() const { return alphabet_; }

  const std::vector<double>& getFrequencies() const { return freq_; }

  void setFrequenciesFromMap(const std::map<int, double>& frequencies);

  unsigned int getNumberOfFrequencies() const { return freq_.size(); }

protected:
  std::vector<double>& getFrequencies_() { return freq_; }
  double& getFreq_(unsigned int i) { return freq_[i]; }
  const double& getFreq_(unsigned int i) const { return freq_[i]; }
  void setFrequencies_(const std::vector<double>& frequencies) { freq_ = frequencies; }
};

/**
 * @brief A generic FrequenciesSet allowing for the estimation of all frequencies.
 *
 * The FrequenciesSet has hence n-1 parameters, where n is the size of the input alphabet.
 */
class FullFrequenciesSet :
  public AbstractFrequenciesSet
{
public:
  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. If the alphabet is a CodonAlphabet, the stop codon
   * frequencies are null.
   */
  FullFrequenciesSet(const Alphabet* alphabet, bool allowNullFreqs = false);
  FullFrequenciesSet(const Alphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false) throw (Exception);

  FullFrequenciesSet* clone() const { return new FullFrequenciesSet(*this); }

public:
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  std::string getName() const {return "Full"; }

protected:
  void fireParameterChanged(const ParameterList& parameters);
};

/**
 * @brief A generic FrequenciesSet for Codon alphabets.
 *
 * It is very similar to FullFrequencySet, but only the non-stop codon
 *   frequencies are parameterized.
 */
class FullCodonFrequenciesSet :
  public virtual CodonFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. If the alphabet is a CodonAlphabet, the stop codon
   * frequencies are null.
   */
  FullCodonFrequenciesSet(const CodonAlphabet* alphabet, bool allowNullFreqs = false);
  FullCodonFrequenciesSet(const CodonAlphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false) throw (Exception);

#ifndef NO_VIRTUAL_COV
  FullCodonFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FullCodonFrequenciesSet(*this); }

public:
  std::string getName() const { return "Full"; }

  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

#ifndef NO_VIRTUAL_COV
  const CodonAlphabet* getAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

protected:
  void fireParameterChanged(const ParameterList& parameters);
};


/**
 * @brief Nucleotide FrequenciesSet using only one parameter, the GC content.
 */
class GCFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  GCFrequenciesSet(const NucleicAlphabet* alphabet) :
    AbstractFrequenciesSet(4, alphabet, "GC.")
  {
    Parameter p("GC.theta", 0.5, &Parameter::PROP_CONSTRAINT_IN);
    addParameter_(p);
    getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
  }

  GCFrequenciesSet(const NucleicAlphabet* alphabet, double theta) :
    AbstractFrequenciesSet(4, alphabet, "GC.")
  {
    Parameter p("GC.theta", theta, &Parameter::PROP_CONSTRAINT_IN);
    addParameter_(p);
    getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
    getFreq_(1) = getFreq_(2) = theta / 2.;
  }

#ifndef NO_VIRTUAL_COV
  GCFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new GCFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  std::string getName() const { return "GC"; }

protected:
  void fireParameterChanged(const ParameterList& parameters);
};

/**
 * @brief Nucleotide FrequenciesSet using three independent parameters to modelize the four frequencies.
 */
class FullNucleotideFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, bool allowNullFreqs = false);

  FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, double theta, double theta1, double theta2, bool allowNullFreqs = false);

#ifndef NO_VIRTUAL_COV
  FullNucleotideFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FullNucleotideFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif

  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  std::string getName() const { return "FullNucleotide"; }

protected:
  void fireParameterChanged(const ParameterList& parameters);
};


/**
 * @brief Protein FrequenciesSet using 19 independent parameters to modelize the 20 frequencies.
 *
 * The parameters are called @f$ \theta_{i \in 1..19} @f$, and are initialized so that all frequencies are equal to  0.005, that is
 * @f[ \theta_i = \frac{0.05}{0.956{i-1}},\quad i = 1..19 @f] or according to a user-specified vector of initial values.
 */
class FullProteinFrequenciesSet :
  public virtual ProteinFrequenciesSet,
  public FullFrequenciesSet
{
public:
  FullProteinFrequenciesSet(const ProteicAlphabet* alphabet, bool allowNullFreqs = false) :
    FullFrequenciesSet(alphabet, allowNullFreqs) {}
  FullProteinFrequenciesSet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs, bool allowNullFreqs = false) throw (Exception) :
    FullFrequenciesSet(alphabet, initFreqs, allowNullFreqs) {}

#ifndef NO_VIRTUAL_COV
  FullProteinFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FullProteinFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};

/**
 * @brief FrequenciesSet to be used with a Markov-modulated substitution model.
 *
 * This implementation uses one parameter per character state frequency.
 * The rate states are assumed to be fixed and are passed as an argument to the constructor, together with a 'regular'
 * FrequenciesSet. The number of parameters hence do not depends on the number of rates used.
 */
class MarkovModulatedFrequenciesSet :
  public AbstractFrequenciesSet
{
private:
  FrequenciesSet* freqSet_;
  std::vector<double> rateFreqs_;

public:
  MarkovModulatedFrequenciesSet(FrequenciesSet* freqSet, const std::vector<double>& rateFreqs) :
    AbstractFrequenciesSet(getAlphabet()->getSize() * rateFreqs.size(), freqSet->getAlphabet(), "MarkovModulated."),
    freqSet_(freqSet),
    rateFreqs_(rateFreqs)
  {
   freqSet_->setNamespace("MarkovModulated." + freqSet_->getNamespace());
    addParameters_(freqSet_->getParameters());
    setFrequencies_(VectorTools::kroneckerMult(rateFreqs, freqSet_->getFrequencies()));
  }

  MarkovModulatedFrequenciesSet(const MarkovModulatedFrequenciesSet& mmfs) :
    AbstractFrequenciesSet(mmfs),
    freqSet_(mmfs.freqSet_->clone()),
    rateFreqs_(mmfs.rateFreqs_)
  {}

  MarkovModulatedFrequenciesSet & operator=(const MarkovModulatedFrequenciesSet& mmfs)
  {
    AbstractFrequenciesSet::operator=(mmfs);
    freqSet_ = mmfs.freqSet_->clone();
    rateFreqs_ = mmfs.rateFreqs_;
    return *this;
  }

  MarkovModulatedFrequenciesSet* clone() const { return new MarkovModulatedFrequenciesSet(*this); }

  virtual ~MarkovModulatedFrequenciesSet() { delete freqSet_; }

public:
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception)
  {
    // Just forward this method to the sequence state frequencies set. This may change in the future...
    freqSet_->setFrequencies(frequencies);
  }

  void fireParameterChanged(const ParameterList& pl)
  {
   freqSet_->matchParametersValues(pl);
    setFrequencies_(VectorTools::kroneckerMult(rateFreqs_, freqSet_->getFrequencies()));
  }

  const FrequenciesSet& getStatesFrequenciesSet() const { return *freqSet_; }

  void setNamespace(const std::string& prefix)
  {
   AbstractFrequenciesSet::setNamespace(prefix);
   freqSet_->setNamespace(prefix + freqSet_->getNamespace());
  }

  std::string getName() const { return "MarkovModulated." + freqSet_->getName(); }
};


/**
 * @brief FrequenciesSet useful for homogeneous and stationary models.
 *
 * This set contains no parameter.
 */
class FixedFrequenciesSet :
  public AbstractFrequenciesSet
{
public:
  FixedFrequenciesSet(const Alphabet* alphabet, const std::vector<double>& initFreqs);

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. If the alphabet is a CodonAlphabet, the stop codon
   * frequencies are null.
   */
  FixedFrequenciesSet(const Alphabet* alphabet);

  FixedFrequenciesSet* clone() const { return new FixedFrequenciesSet(*this); }

public:
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  std::string getName() const { return "Fixed"; }

protected:
  void fireParameterChanged(const ParameterList& parameters) {}
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models, nucleotide implementation
 *
 * This set contains no parameter.
 */
class FixedNucleotideFrequenciesSet :
  public virtual NucleotideFrequenciesSet,
  public FixedFrequenciesSet
{
public:
  FixedNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, const std::vector<double>& initFreqs) :
    FixedFrequenciesSet(alphabet, initFreqs) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedNucleotideFrequenciesSet(const NucleicAlphabet* alphabet) :
    FixedFrequenciesSet(alphabet) {}

#ifndef NO_VIRTUAL_COV
  FixedNucleotideFrequenciesSet*
#else
  NucleotideFrequenciesSet*
#endif
  clone() const { return new FixedNucleotideFrequenciesSet(*this); }

#ifndef NO_VIRTUAL_COV
  const NucleicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const NucleicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models, protein implementation
 *
 * This set contains no parameter.
 */
class FixedProteinFrequenciesSet :
  public virtual ProteinFrequenciesSet,
  public FixedFrequenciesSet
{
public:
  FixedProteinFrequenciesSet(const ProteicAlphabet* alphabet, const std::vector<double>& initFreqs) :
    FixedFrequenciesSet(alphabet, initFreqs) {}

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet.
   */
  FixedProteinFrequenciesSet(const ProteicAlphabet* alphabet) :
    FixedFrequenciesSet(alphabet) {}

#ifndef NO_VIRTUAL_COV
  FixedProteinFrequenciesSet*
#else
  FixedFrequenciesSet*
#endif
  clone() const { return new FixedProteinFrequenciesSet(*this); }

#ifndef NO_VIRTUAL_COV
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models, codon implementation
 *
 * This set contains no parameter.
 */
class FixedCodonFrequenciesSet :
  public virtual CodonFrequenciesSet,
  public AbstractFrequenciesSet
{
public:
  FixedCodonFrequenciesSet(const CodonAlphabet* alphabet, const std::vector<double>& initFreqs);

  /**
   * @brief Construction with uniform frequencies on the letters of
   * the alphabet. The stop codon frequencies are null.
   */
  FixedCodonFrequenciesSet(const CodonAlphabet* alphabet);

#ifndef NO_VIRTUAL_COV
  FixedCodonFrequenciesSet*
#else
  Clonable*
#endif
  clone() const { return new FixedCodonFrequenciesSet(*this); }

public:
#ifndef NO_VIRTUAL_COV
  const CodonAlphabet* getAlphabet() const
  {
    return dynamic_cast<const CodonAlphabet*>(AbstractFrequenciesSet::getAlphabet());
  }
#endif
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  std::string getName() const { return "Fixed"; }

protected:
  void fireParameterChanged(const ParameterList& parameters) {}
};


/*********************************************************************/
/****   Frequencies Set in Words *****/
/*********************************************************************/


/**
 * @brief Frequencies in words computed from the  frequencies on
 * letters. The parameters are the parameters of the Frequencies on
 * letters.
 * The WordFrequenciesSet owns the FrequenciesSet* it is built on.
 * Interface class.
 * @author Laurent Guéguen
 */
class WordFrequenciesSet :
  public AbstractFrequenciesSet
{
protected:
  unsigned int getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector);

public:
  WordFrequenciesSet(unsigned int size, const Alphabet* palph);

  virtual ~WordFrequenciesSet();

  /**
   *@ brief Return the n-th FrequenciesSet*
   **/

  const FrequenciesSet& getFrequenciesSetForLetter(unsigned int i) const;

  /**
   *@ brief Return the length of the words
   **/

  virtual unsigned int getLength() const;
};


/**
 * @brief the Frequencies in words are the product of Independent Frequencies in letters
 * @author Laurent Guéguen
 */

class WordFromIndependentFrequenciesSet :
  public WordFrequenciesSet
{
private:
  std::vector<FrequenciesSet*> vFreq_;
  std::vector<std::string> vNestedPrefix_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and a vector of different FrequenciesSet*.
   * Throws an Exception if their lengths do not match.
   */
  WordFromIndependentFrequenciesSet(const WordAlphabet*, const std::vector<FrequenciesSet*>&) throw (Exception);

  WordFromIndependentFrequenciesSet(const WordFromIndependentFrequenciesSet& iwfs);

  ~WordFromIndependentFrequenciesSet();

  WordFromIndependentFrequenciesSet& operator=(const WordFromIndependentFrequenciesSet& iwfs);

  WordFromIndependentFrequenciesSet* clone() const { return new WordFromIndependentFrequenciesSet(*this); }

public:
  void fireParameterChanged(const ParameterList& pl);

  void updateFrequencies();

  /**
   *@ brief Independent letter frequencies from given word frequencies.
   * The frequencies of a letter at a position is the sum of the
   *    frequencies of the words that have this letter at this
   *    position.
   */
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  /**
   *@ brief Return the n-th FrequenciesSet*
   **/
  const FrequenciesSet& getFrequenciesSetForLetter(unsigned int i) const { return *vFreq_[i]; }

  /**
   *@ brief Return the length of the words
   **/

  virtual unsigned int getLength() const;

  void setNamespace(const std::string&);

  std::string getName() const;
};

/**
 * @brief the Frequencies in words are the product of the frequencies for a unique FrequenciesSet in letters
 * @author Laurent Guéguen
 */

class WordFromUniqueFrequenciesSet :
  public WordFrequenciesSet
{
private:
  FrequenciesSet* pFreq_;
  std::string NestedPrefix_;
  unsigned int length_;

public:
  /**
   * @brief Constructor from a WordAlphabet* and an AbstractFrequenciesSet* repeated as
   *  many times as the length of the words.
   */
  WordFromUniqueFrequenciesSet(const WordAlphabet* pWA, FrequenciesSet* pabsfreq);

  WordFromUniqueFrequenciesSet(const WordFromUniqueFrequenciesSet& iwfs);

  WordFromUniqueFrequenciesSet& operator=(const WordFromUniqueFrequenciesSet& iwfs);

  ~WordFromUniqueFrequenciesSet();

  WordFromUniqueFrequenciesSet* clone() const { return new WordFromUniqueFrequenciesSet(*this); }

public:
  void fireParameterChanged(const ParameterList& pl);

  /**
   *@ brief letter frequencies from given word frequencies. The
   * frequencies of a letter at a position is the sum of the
   * frequencies of the words that have this letter at this position.
   * The frequencies of each letter is the average of the frequencies
   * of that letter at all positions.
   */
  void setFrequencies(const std::vector<double>& frequencies) throw (DimensionException, Exception);

  void updateFrequencies();

  /**
   *@ brief Return the n-th FrequenciesSet*
   **/
  const FrequenciesSet& getFrequenciesSetForLetter(unsigned int i){ return *pFreq_; }

  unsigned int getLength() const { return length_; }

  void setNamespace(const std::string& prefix);

  std::string getName() const;
};

} // end of namespace bpp.

#endif // _FREQUENCIESSET_H_

