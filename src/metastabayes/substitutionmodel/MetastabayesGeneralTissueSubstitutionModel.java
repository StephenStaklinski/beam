package metastabayes.substitutionmodel;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

@Description("Substitution model with structure indicators to setup shared rate paramaters between tissues." +
                "This code is loosely based on the FixedTreeAnalysis package by Andrew Bouckaert")

public class MetastabayesGeneralTissueSubstitutionModel extends SVSGeneralSubstitutionModel {
	
    public Input<IntegerParameter> structure = new Input<IntegerParameter>("matrixStructure",
            "integer indices to indicate structure of shared transition rate matrix parameters between diferent tissues", Validate.REQUIRED);

    private IntegerParameter matrixStructure;
    protected BooleanParameter rateIndicator;

	@Override
    public void initAndValidate(){
        
        updateMatrix = true;
        matrixStructure = structure.get();
        nrOfStates = (int) Math.sqrt(matrixStructure.getDimension());
		eigenSystem = createEigenSystem();
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];
        frequencies = frequenciesInput.get();
        rateIndicator = indicator.get();

    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        // new code to allow specification of the rate matrix structure for parameter sharing
        int n = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[matrixStructure.getValue(n)];
                n++;
            }
        }
        
        // retained from svsGeneralSubstitutionModel code
        double [] fFreqs = frequencies.getFreqs();

        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }

        // bring in frequencies
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] *= fFreqs[j];
                rateMatrix[j][i] *= fFreqs[i];
            }
        }
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }
        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * fFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
    }
    
}
