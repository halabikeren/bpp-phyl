//
// File: FrequenciesSet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
Copyright or <A9> or Copr. CNRS, (November 16, 2004)

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

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/NucleicAlphabet.h>

/**
 * @brief Parametrize a set of state frequencies.
 */
class FrequenciesSet:
  public virtual Parametrizable
{
  public:
#ifndef NO_VIRTUAL_COV
    FrequenciesSet * clone() const = 0;
#endif
  public:
    /**
     * @return The alphabet associated to this set.
     */  
    virtual const Alphabet * getAlphabet() const = 0;

    /**
     * @return The frequencies values of the set.
     */ 
    virtual vector<double> getFrequencies() const = 0;

    /**
     * @return The nimer of parameters associated to these frequencies.
     */ 
    virtual unsigned int getNumberOfParameters() const = 0;
};

/**
 * @brief Basic implementation of the FrequenciesSet interface.
 */
class AbstractFrequenciesSet:
  public virtual FrequenciesSet,
  public AbstractParametrizable
{
  protected:
    const Alphabet * _alphabet;
    vector<double> _freq;

  public:
    AbstractFrequenciesSet(const Alphabet * alphabet): _alphabet(alphabet) {}
  
  public:
    const Alphabet * getAlphabet() const { return _alphabet; }
    vector<double> getFrequencies() const { return _freq; }
    unsigned int getNumberOfParameters() const { return _parameters.size(); }
};

/**
 * @brief FrequenciesSet using one parameter per frequency.
 */
class FullFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    FullFrequenciesSet(const Alphabet * alphabet, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      _freq.resize(alphabet->getSize());
      for(unsigned int i = 0; i < alphabet->getSize(); i++)
      {
        _parameters.addParameter(Parameter(prefix + alphabet->intToChar((int)i), 1. / alphabet->getSize(), &Parameter::PROP_CONSTRAINT_IN));
        _freq[i] = 1. / alphabet->getSize();
      }
    }
    FullFrequenciesSet(const Alphabet * alphabet, const vector<double> & initFreqs, const string & prefix = "") throw (Exception):
      AbstractFrequenciesSet(alphabet)
    {
      if(initFreqs.size() != alphabet->getSize())
        throw Exception("FullFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
      double sum = 0.0;
      for(unsigned int i = 0; i < initFreqs.size(); i++)
      {
        sum += initFreqs[i];
      }
      if(fabs(1-sum) > 0.00000000000001)
      {
        throw Exception("Root frequencies must equal 1.");
      }
      _freq.resize(alphabet->getSize());
      for(unsigned int i = 0; i < alphabet->getSize(); i++)
      {
        _parameters.addParameter(Parameter(prefix + alphabet->intToChar((int)i), initFreqs[i], &Parameter::PROP_CONSTRAINT_IN));
        _freq[i] = initFreqs[i];
      }
    }
#ifndef NO_VIRTUAL_COV
    FullFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FullFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      for(unsigned int i = 0; i < _alphabet->getSize(); i++)
      {
        _freq[i] = _parameters[i]->getValue();
      }
    }
};

/**
 * @brief Nucleotide FrequenciesSet using only one parameter, the GC content.
 */
class GCFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    GCFrequenciesSet(const NucleicAlphabet * alphabet, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      _freq.resize(4);
      _parameters.addParameter(Parameter(prefix + "theta", 0.5, &Parameter::PROP_CONSTRAINT_IN));
      _freq[0] = _freq[1] = _freq[2] = _freq[3] = 0.25;
    }
    GCFrequenciesSet(const NucleicAlphabet * alphabet, double theta, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      _freq.resize(4);
      _parameters.addParameter(Parameter(prefix + "theta", theta, &Parameter::PROP_CONSTRAINT_IN));
      _freq[0] = _freq[3] = (1. - theta) / 2.;
      _freq[1] = _freq[2] = theta / 2.;
    }
#ifndef NO_VIRTUAL_COV
    GCFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new GCFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      double theta = _parameters[0]->getValue();
      _freq[0] = _freq[3] = (1. - theta) / 2.;
      _freq[1] = _freq[2] = theta / 2.;
    }
};

/**
 * @brief Nucleotide FrequenciesSet using three parameters to modelize the four frequencies.
 */
class FullNAFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    FullNAFrequenciesSet(const NucleicAlphabet * alphabet, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      _freq.resize(4);
      _parameters.addParameter(Parameter(prefix + "theta" , 0.5, &Parameter::PROP_CONSTRAINT_EX));
      _parameters.addParameter(Parameter(prefix + "theta1", 0.5, &Parameter::PROP_CONSTRAINT_EX));
      _parameters.addParameter(Parameter(prefix + "theta2", 0.5, &Parameter::PROP_CONSTRAINT_EX));
      _freq[0] = _freq[1] = _freq[2] = _freq[3] = 0.25;
    }
    FullNAFrequenciesSet(const NucleicAlphabet * alphabet, double theta, double theta1, double theta2, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      _freq.resize(4);
      _parameters.addParameter(Parameter(prefix + "theta" , theta , &Parameter::PROP_CONSTRAINT_EX));
      _parameters.addParameter(Parameter(prefix + "theta1", theta1, &Parameter::PROP_CONSTRAINT_EX));
      _parameters.addParameter(Parameter(prefix + "theta2", theta2, &Parameter::PROP_CONSTRAINT_EX));
      _freq[0] = theta1 * (1. - theta);
      _freq[1] = (1 - theta2) * theta;
      _freq[2] = theta2 * theta;
      _freq[3] = (1 - theta1) * (1. - theta);
    }
#ifndef NO_VIRTUAL_COV
    FullNAFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FullNAFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      double theta  = _parameters[0]->getValue();
      double theta1 = _parameters[1]->getValue();
      double theta2 = _parameters[2]->getValue();
      _freq[0] = theta1 * (1. - theta);
      _freq[1] = (1 - theta2) * theta;
      _freq[2] = theta2 * theta;
      _freq[3] = (1 - theta1) * (1. - theta);
    }
};

/**
 * @brief FrequenciesSet to be used with a Markov-modulated substitution model.
 * 
 * This implementation uses one parameter per character state frequency.
 * The rate states are assumed to be fixed and are passed as an argument to the constructor, together with a 'regular'
 * FrequenciesSet. The number of parameters hence do not depends on the number of rates used.
 */
class MarkovModulatedFrequenciesSet:
  public AbstractFrequenciesSet
{
  protected:
    FrequenciesSet * _freqSet;
    vector<double> _rateFreqs;
  public:
    MarkovModulatedFrequenciesSet(const MarkovModulatedFrequenciesSet & mmfs):
      AbstractFrequenciesSet(mmfs)
    {
      _freqSet = dynamic_cast<FrequenciesSet *>(mmfs._freqSet->clone());
      _rateFreqs = mmfs._rateFreqs;
    }
    MarkovModulatedFrequenciesSet & operator=(const MarkovModulatedFrequenciesSet & mmfs)
    {
      AbstractFrequenciesSet::operator=(mmfs);
      _freqSet = dynamic_cast<FrequenciesSet *>(mmfs._freqSet->clone());
      _rateFreqs = mmfs._rateFreqs;
      return *this;
    }
    MarkovModulatedFrequenciesSet(FrequenciesSet * freqSet, const vector<double> & rateFreqs):
      AbstractFrequenciesSet(freqSet->getAlphabet()), _rateFreqs(rateFreqs)
    {
      _freq.resize(_alphabet->getSize() * rateFreqs.size());
      _parameters.addParameters(_freqSet->getParameters());
      _freq = VectorTools::kroneckerMult(rateFreqs, _freqSet->getFrequencies());
    }
#ifndef NO_VIRTUAL_COV
    MarkovModulatedFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new MarkovModulatedFrequenciesSet(*this); }

    virtual ~MarkovModulatedFrequenciesSet() { delete _freqSet; }

  public:
    void fireParameterChanged(const ParameterList & pl)
    {
      _freqSet->matchParametersValues(pl);
      _freq = VectorTools::kroneckerMult(_rateFreqs, _freqSet->getFrequencies());
    }
};

/**
 * @brief FrequenciesSet useful for homogeneous and stationary models.
 *
 * This set contains no parameter.
 */
class FixedFrequenciesSet:
  public AbstractFrequenciesSet
{
  public:
    FixedFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs, const string & prefix = ""):
      AbstractFrequenciesSet(alphabet)
    {
      double sum = 0.0;
      for(unsigned int i = 0; i < initFreqs.size(); i++)
      {
        sum += initFreqs[i];
      }
      if(fabs(1-sum) > 0.00000000000001)
      {
        throw Exception("Root frequencies must equal 1.");
      }
      _freq = initFreqs;
    }
#ifndef NO_VIRTUAL_COV
    FixedFrequenciesSet *
#else
    Clonable *
#endif
    clone() const { return new FixedFrequenciesSet(*this); }

  public:
    void fireParameterChanged(const ParameterList & pl) {}
};

#endif //_FREQUENCIESSET_H_

