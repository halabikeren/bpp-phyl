// From bpp-core:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/NumTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl
#include <Bpp/Phyl/Model/Character/CharacterSubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/SingleRateModel.h>
#include <Bpp/Phyl/Model/Character/RatePerEntryModel.h>
#include <Bpp/Phyl/Model/Character/RatePerExitModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairSymModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairModel.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

// from std
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

enum RatesParameterizationScope { noRate, singleRate, ratePerEntry, ratePerExit, ratePerPairSym, ratePerPair};
enum FreqParameterizationScope {fixed, parameterized};

double getTransitionProb(size_t i, size_t j, double d, RowMatrix<double>& rates, vector<double>& freqs)
{
    double e = exp(rates(0,0)*d);
    if ((i == 0) && (j == 0))
    {
        return freqs[0] + freqs[1] * e;
    }
    if ((i == 0) && (j == 1))
    {
        return freqs[1] * (1 - e);
    }
    if ((i == 1) && (j == 0))
    {
        return freqs[0] * (1 - e);
    }
    else
    {
        return freqs[1] + freqs[0] * e;
    }
}

void test_binary_model_settings(CharacterSubstitutionModel& model, RatesParameterizationScope ratesParameterizationScope, RowMatrix<double>& rates, vector<double>& freqs)
{
    // test rate matrix
    RowMatrix<double> rateMatrix = model.getGenerator();
    for (size_t i=0; i<rateMatrix.getNumberOfRows(); ++i)
    {
        for (size_t j=0; j<rateMatrix.getNumberOfColumns(); ++j)
        {
            double expectedRate =  rates(i, j) * ((i == j) ? (1-freqs[i]) : freqs[j]);
            if (abs(rateMatrix(i, j)-expectedRate) > 0.001)
            {
                cerr << "accepted rate value in position (" << i << ", " << j << ") is " << rateMatrix(i, j) << " instead of " << expectedRate << endl;
                exit(1);
            }
        }
    }

    // test eigen values
    vector<double> eigenvalues = model.getEigenValues();
    vector<double> expectedEigenValues(2);
    expectedEigenValues[0] = expectedEigenValues[0] = -rates(1,0) * freqs[0] - rates(0,1) * freqs[1];
    expectedEigenValues[1] = 0;
    if (!((abs(eigenvalues[0]-expectedEigenValues[0]<0.001) && (abs(eigenvalues[1]-expectedEigenValues[1]<0.001)))
        | (abs(eigenvalues[1]-expectedEigenValues[0]<0.001) && (abs(eigenvalues[0]-expectedEigenValues[1]<0.001)))))
    {
        cerr << "accepted eigvalues are not as expected. Accepted are: (" << eigenvalues[0] << ", " << eigenvalues[1] << ") while expected are (" << expectedEigenValues[0] << ", " << expectedEigenValues[1] << ")" << endl;
        exit(1);
    }
    
    // test transition matrix - TO DO: this currently only fits the noRate and singleRate parameterization scopes, and I need to add handling of the other cases here by manually computing the closed fomulaas for eigenvectos and computing the eigen decomposition of them
    if (ratesParameterizationScope == RatesParameterizationScope::noRate || ratesParameterizationScope == RatesParameterizationScope::singleRate)
    {
        double d = 1.;
        RowMatrix<double> transitionMatrix = model.getPij_t(d);
        for (size_t i=0; i<transitionMatrix.getNumberOfRows(); ++i)
        {
            for (size_t j=0; j<transitionMatrix.getNumberOfColumns(); ++j)
            {
                double expectedProb = getTransitionProb(i, j, d, rates, freqs);
                if (abs(transitionMatrix(i, j)-expectedProb) > 0.001)
                {
                    cerr << "accepted transition probability value in position (" << i << ", " << j << ") is " << transitionMatrix(i, j) << " instead of " << expectedProb << endl;
                    exit(1);
                }
            }
        }
    }

    // make sure you are able to change rate and frequency values
    const ParameterList& modelParamters = model.getParameters();
    for (size_t p=0; p<modelParamters.size(); ++p)
    {
        double prevParamValue = modelParamters[p].getValue();
        double currParamValue = prevParamValue+0.1;
        const string& parName = model.getParameterNameWithoutNamespace(modelParamters[p].getName());
        model.setParameterValue(parName, currParamValue);
        if (abs(modelParamters[p].getValue()-currParamValue)>0.0001)
        {
            cerr << "failed to modify parameter value for " << modelParamters[p].getName() << ". Instead of " << currParamValue << " the parameter value is " << modelParamters[p].getValue() << endl;
            exit(1);
        }
    }

    // make sure you are able to limit the range of rates
    for (size_t p=0; p<modelParamters.size(); ++p)
    {
        const string& parName = model.getParameterNameWithoutNamespace(modelParamters[p].getName());
        if ((parName.find("rate") != std::string::npos))
        {
            model.setParamBounds(parName, 0.1, 50);
            const shared_ptr<Constraint> constraint = modelParamters[p].getConstraint();
            if ((!constraint->includes(0.1, 50)) | (constraint->includes(0.05, 50)) | (constraint->includes(0.1, 51)))
            {
                cerr << "failed to modify boundaries of rate parameter " << parName << " to (0.1, 10)" << endl;
                exit(1);
            }
        }
    }


}

void test_parameterization_scope(CharacterSubstitutionModel& model, RatesParameterizationScope ratesParameterizationScope, FreqParameterizationScope freqParameterizationScope)
{
    const ParameterList modelParams = model.getParameters();
    size_t statesNum = model.getNumberOfStates();

    size_t expectedParametersNum = 0;
    if (freqParameterizationScope == FreqParameterizationScope::parameterized)
        expectedParametersNum += statesNum-1;
    if (ratesParameterizationScope == RatesParameterizationScope::singleRate)
        expectedParametersNum += 1;
    if ((ratesParameterizationScope == RatesParameterizationScope::ratePerEntry) | (ratesParameterizationScope == RatesParameterizationScope::ratePerExit))
        expectedParametersNum += statesNum;
    if (ratesParameterizationScope == RatesParameterizationScope::ratePerPairSym)
        expectedParametersNum += statesNum * (statesNum-1) / 2;
    if (ratesParameterizationScope == RatesParameterizationScope::ratePerPair)
        expectedParametersNum += statesNum * (statesNum-1);
    
    if (modelParams.size() != expectedParametersNum)
    {
        cerr << "Invalid number of model parameters " << modelParams.size() << " with parameters: " << endl;
        for (size_t p=0; p<modelParams.size(); ++p)
            cerr << modelParams[p].getName() << "=" << modelParams[p].getValue() << endl;
        exit(1);
    }
}

void test_binary_models(const IntegerAlphabet* alpha)
{
    cout << "case 1: binary models" << endl;
    RowMatrix<double> rateVals(2, 2);
    vector<double> freqVals(2);
    freqVals[0] = freqVals[1] = 0.5;
    shared_ptr<IntegerFrequencySet> freqs;
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));

    // no rate no freqs model
    cout << "case 1.1: rate = 1, equal frequencies: 0.5, 0.5" << endl;
    for (size_t i=0; i<2; ++i)
    {
        for (size_t j=0; j<2; ++j)
        {
            rateVals(i,j) = ((i == j) ?  -1: 1);
        }
    }
    CharacterSubstitutionModel noRatenoFreqsModel("mk", alpha);
    noRatenoFreqsModel.setFrequencySet(*freqs);
    test_binary_model_settings(noRatenoFreqsModel, RatesParameterizationScope::noRate, rateVals, freqVals);
    test_parameterization_scope(noRatenoFreqsModel, RatesParameterizationScope::noRate, FreqParameterizationScope::fixed);

    // no rate model
    cout << "case 1.2: rate = 1, unequal frequencies: 0.1, 0.9" << endl;
    freqVals[0] = 0.1;
    freqVals[1] = 1-freqVals[0];
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));
    CharacterSubstitutionModel noRateModel("mk", alpha, freqs, false);
    noRateModel.setFrequencySet(*freqs);
    test_binary_model_settings(noRateModel, RatesParameterizationScope::singleRate, rateVals, freqVals);
    test_parameterization_scope(noRateModel, RatesParameterizationScope::noRate, FreqParameterizationScope::parameterized);

    // single rate model
    cout << "case 1.3: rate = 42., unequal frequencies: 0.1, 0.9" << endl;
    double globalRate = 42.;
    for (size_t i=0; i<2; ++i)
    {
        for (size_t j=0; j<2; ++j)
        {
            rateVals(i,j) = ((i == j) ?  -globalRate: globalRate);
        }
    }
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals)); // frequencies set must be reset upon each usage to avoid concatanation of namespaces of prevous models
    SingleRateModel oneRateModel(alpha, freqs, false);
    oneRateModel.setParameterValue("global_rate", rateVals(0,1));
    oneRateModel.setFrequencySet(*freqs);
    test_binary_model_settings(dynamic_cast<CharacterSubstitutionModel&>(oneRateModel), RatesParameterizationScope::singleRate, rateVals, freqVals);
    test_parameterization_scope(oneRateModel, RatesParameterizationScope::singleRate, FreqParameterizationScope::parameterized);

    // rate per entry model
    cout << "case 1.4: entry rates = [2., 42.], unequal frequencies: 0.1, 0.9" << endl;
    rateVals(1, 0) = 2.;
    rateVals(0, 1) = 42.;
    rateVals(0, 0) = - rateVals(0, 1);
    rateVals(1, 1) = - rateVals(1, 0);
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));
    RatePerEntryModel ratePerEntryModel(alpha, freqs, false);
    ratePerEntryModel.setParameterValue("entry_rate_0", rateVals(1, 0));
    ratePerEntryModel.setParameterValue("entry_rate_1", rateVals(0, 1));
    ratePerEntryModel.setFrequencySet(*freqs);
    test_binary_model_settings(dynamic_cast<CharacterSubstitutionModel&>(ratePerEntryModel), RatesParameterizationScope::ratePerEntry, rateVals, freqVals);
    test_parameterization_scope(ratePerEntryModel, RatesParameterizationScope::ratePerEntry, FreqParameterizationScope::parameterized);


    // rate per exit model
    cout << "case 1.5: exit rates = [2., 42.], unequal frequencies: 0.1, 0.9" << endl;
    rateVals(1, 0) = 2.;
    rateVals(0, 1) = 42.;
    rateVals(0, 0) = - rateVals(0, 1);
    rateVals(1, 1) = - rateVals(1, 0);
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));
    RatePerExitModel ratePerExitModel(alpha, freqs, false);
    ratePerExitModel.setParameterValue("exit_rate_0", rateVals(0, 1));
    ratePerExitModel.setParameterValue("exit_rate_1", rateVals(1, 0));
    ratePerEntryModel.setFrequencySet(*freqs);
    test_binary_model_settings(dynamic_cast<CharacterSubstitutionModel&>(ratePerExitModel), RatesParameterizationScope::ratePerExit, rateVals, freqVals);
    test_parameterization_scope(ratePerExitModel, RatesParameterizationScope::ratePerExit, FreqParameterizationScope::parameterized);

    // rate per pair symmetric model
    cout << "case 1.6: rates(0,1)=rates(1,0)=2. unequal frequencies: 0.1, 0.9" << endl;
    rateVals(0, 1) = rateVals(1, 0) = 2.;
    rateVals(0, 0) = - rateVals(0, 1);
    rateVals(1, 1) = - rateVals(1, 0);
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));
    RatePerPairSymModel ratePerPairSymModel(alpha, freqs, false);
    ratePerPairSymModel.setParameterValue("rate_0_1", rateVals(0, 1));
    ratePerEntryModel.setFrequencySet(*freqs);
    test_binary_model_settings(dynamic_cast<CharacterSubstitutionModel&>(ratePerPairSymModel), RatesParameterizationScope::ratePerPairSym, rateVals, freqVals);
    test_parameterization_scope(ratePerPairSymModel, RatesParameterizationScope::ratePerPairSym, FreqParameterizationScope::parameterized);

    // rate per pair model
    cout << "case 1.7: rates(0,1)=2., rates(1,0)=42. unequal frequencies: 0.1, 0.9" << endl;
    rateVals(0, 1) = 2.,
    rateVals(1, 0) = 42.;
    rateVals(0, 0) = - rateVals(0, 1);
    rateVals(1, 1) = - rateVals(1, 0);
    freqs.reset(new FullIntegerFrequencySet(alpha, freqVals));
    RatePerPairModel ratePerPairModel(alpha, freqs, false);
    ratePerPairModel.setParameterValue("rate_0_1", rateVals(0, 1));
    ratePerPairModel.setParameterValue("rate_1_0", rateVals(1, 0));
    ratePerEntryModel.setFrequencySet(*freqs);
    test_binary_model_settings(dynamic_cast<CharacterSubstitutionModel&>(ratePerPairModel), RatesParameterizationScope::ratePerPair, rateVals, freqVals);
    test_parameterization_scope(ratePerPairModel, RatesParameterizationScope::ratePerPair, FreqParameterizationScope::parameterized);
}

int main() 
{
    try
    {
        
        const IntegerAlphabet* int_alphabet = new IntegerAlphabet(1);
        test_binary_models(int_alphabet);
    }
        catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}