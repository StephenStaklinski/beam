package metastabayes.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.Parameter;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

@Description("Substitution model with one rate parameter fixed between all tissues")

public class OneRateAllTissues extends SVSGeneralSubstitutionModel {
	
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
        
        // custom rate matrix setup with all rates fixed to a single a parameter
        int count = 0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            for (int j = i+1; j <  nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[count];
               	rateMatrix[j][i] = relativeRates[count];
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
