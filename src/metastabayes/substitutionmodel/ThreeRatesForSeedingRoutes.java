package metastabayes.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.Parameter;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

@Description("Substitution model with three rate parameter whre one is fixed for leaving the primary, one fixed for met to met transitions," +
			"and one fixed for returning to the primary")

public class ThreeRatesForSeedingRoutes extends SVSGeneralSubstitutionModel {
	
	@Override
    public void initAndValidate(){

		
        frequencies = frequenciesInput.get();
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
		eigenSystem = createEigenSystem();
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];
        // indicator.get();
    }
    
    @Override
    public void setupRelativeRates() {

        Function rates = this.ratesInput.get();
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i);
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        double [] fFreqs = frequencies.getFreqs();
        
        // custom asymmetric rate matrix with 3 parameters with rates fixed for one leaving primary, one met to met, and one met to primary.
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            // row 1 is all parameter0 for leaving the primary            
            // column 1 is all parameter2 for returning to the parameter
            // all others are parameter1 for met to met transitions
            for (int j = 0; j < i; j++) {
                if (i == 0) {
                	rateMatrix[i][j] = relativeRates[0];
                } else if (j == 0) {
                	rateMatrix[i][j] = relativeRates[2];
                } else {
                	rateMatrix[i][j] = relativeRates[1];
                }
            }
            for (int j = i + 1; j < nrOfStates; j++) {
                if (i == 0) {
                	rateMatrix[i][j] = relativeRates[0];
                } else if (j == 0) {
                	rateMatrix[i][j] = relativeRates[2];
                } else {
                	rateMatrix[i][j] = relativeRates[1];
                }
            }
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
    } // setupRateMatrix
    
    
    @Override
    protected boolean requiresRecalculation() {
    
    	Function v = ratesInput.get(); 
    	if (v instanceof Parameter<?>) {
    		Parameter.Base<?> p = (Parameter.Base<?>) v;
    		if (p.somethingIsDirty()) {
    			if (frequencies.isDirtyCalculation()) {
			    	return super.requiresRecalculation();
    			}
    			// no calculation is affected
    			return false;
    		}
    	}

    	return super.requiresRecalculation();
    }
    
}
