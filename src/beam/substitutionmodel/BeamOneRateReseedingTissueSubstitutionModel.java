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

@Description("Substitution model that fixes reseeding rates to one free parameter.")

public class BeamOneRateReseedingTissueSubstitutionModel extends ComplexSubstitutionModel {

	@Override
    public void initAndValidate(){
        super.initAndValidate();

        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];

        // Verify the number of input rates is correct
        if (relativeRates.length != 3 ) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to 3 (one primary seeding rate, one met to met seeding rate, and one primary reseeding rate), but it is " 
                + relativeRates.length 
                + ". Check the dimension of the input rate parameters."
            );
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                if ( i == j) {
                    // diagonal elements will be calculated later
                    rateMatrix[i][j] = 0;
                } else if ( i == 0) {
                    // primary seeding rate
                    rateMatrix[i][j] = relativeRates[0];
                } else if ( j == 0) {
                    // reseeding rate
                    rateMatrix[i][j] = relativeRates[2];
                } else {
                    // met-to-met seeding rate
                    rateMatrix[i][j] = relativeRates[1];
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

        // Normalize the rate matrix to one subsitution per unit time, since doing so as xml input is challenging
        double total = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            total += -rateMatrix[i][i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / total;
            }
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
