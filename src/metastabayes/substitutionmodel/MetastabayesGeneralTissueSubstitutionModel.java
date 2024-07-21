package metastabayes.substitutionmodel;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

@Description("Substitution model with structure indicators to setup shared rate paramaters between tissues.")

public class MetastabayesGeneralSubstitutionModel extends SVSGeneralSubstitutionModel {
	
    public Input<IntegerParameter> structure = new Input<IntegerParameter>("matrixStructure",
            "integer indices to indicate structure of shared transition rate matrix parameters between diferent tissues", Validate.REQUIRED);

    private IntegerParameter matrixStructure;
    protected BooleanParameter rateIndicator;

	@Override
    public void initAndValidate(){
        
        updateMatrix = true;
        matrixStructure = structure.get();
        nrOfStates = (int) Math.sqrt(matrixStructure.getDimension());
        System.out.println("nr of states" + nrOfStates);
		eigenSystem = createEigenSystem();
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];
        frequencies = frequenciesInput.get();
        rateIndicator = indicator.get();

    }

    @Override
    public void setupRelativeRates() {

        Function rates = this.ratesInput.get();
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i) * (rateIndicator.getValue(i)?1.:0.);
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {
        
        int n = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[matrixStructure.getValue(n)];
                n++;
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
    }
    

    @Override 
    protected boolean requiresRecalculation() {
    	// if the rate is only dirty for a value that the indicators block out,
    	// no recalculation is required, so check this first.
    	Function v = ratesInput.get(); 
    	if (v instanceof Parameter<?>) {
    		Parameter.Base<?> p = (Parameter.Base<?>) v;
    		if (p.somethingIsDirty()) {
        		Parameter<Boolean> indicator2 = indicator.get(); 
    			for (int i = 0; i < p.getDimension(); i++) {
    				if (indicator2.getValue(i) && p.isDirty(i)) {
    			    	return super.requiresRecalculation();
    				}
    			}
    			// no calculation is affected
    			return false;
    		}
    	}
        updateMatrix = true;
        return true;
    }
    
}
