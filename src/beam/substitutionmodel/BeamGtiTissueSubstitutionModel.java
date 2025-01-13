package beam.substitutionmodel;

import java.util.HashMap; 

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;

/**
 * @author Stephen Staklinski
 **/

@Description("General GTI substitution model implementation.")

public class BeamGtiTissueSubstitutionModel extends ComplexSubstitutionModel {

	@Override
    public void initAndValidate(){
        super.initAndValidate();

        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];

        // Verify the number of input rates is correct
        if (relativeRates.length != nrOfStates * nrOfStates - nrOfStates) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to ((nrOfStates * nrOfStates) - nrOfStates) for all off diagonal rates, but it is " 
                + relativeRates.length 
                + ". Check the dimension of the input rate parameters."
            );
        }

        // Verify that the input rates are normalized to sum to one substitution per unit time in the input xml file
        double sumOfRates = 0.0;
        for (double rate : relativeRates) {
            sumOfRates += rate;
        }
        if (Math.abs(sumOfRates - 1.0) > 1e-6) {
            throw new IllegalArgumentException(
                "The sum of the input rates must be 1.0 to fix the rate matrix to one substitution per unit time, but it is " 
                + sumOfRates 
                + ". Check that the input rate parameters sum to 1.0 and that the operator maintains the sum during MCMC proposals."
            );
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {
        
        int count = 0;
        for (int i=0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                if (i == j) {
                    // diagonal is setup later
                    rateMatrix[i][j] = 0;
                } else {
                    rateMatrix[i][j] = relativeRates[count];
                    count += 1;
                }
            }
        }
        
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum;
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }

}
