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

@Description("Substitution model that fixes reseeding rates to zero.")

public class BeamNoReseedingTissueSubstitutionModel extends ComplexSubstitutionModel {

	@Override
    public void initAndValidate(){
        super.initAndValidate();

        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];

        // Verify the number of input rates is correct
        int nrInputRates = ratesInput.get().getDimension();
        if (nrInputRates != 2 ) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to 2 (one primary seeding rate and one met to met seeding rate), but it is " 
                + nrInputRates 
                + ". Check the dimension of the input rate parameters."
            );
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        // set up rate matrix
        double[] rowSums = new double[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
            if (i == j || j == 0) {
                rateMatrix[i][j] = 0;
            } else if (i == 0) {
                rateMatrix[i][j] = relativeRates[0];
            } else {
                rateMatrix[i][j] = relativeRates[1];
            }
            rowSums[i] += rateMatrix[i][j];
            }
        }

        // set up diagonal and normalize the rate matrix to one substitution per unit time
        double total = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = -rowSums[i];
            total += rowSums[i];
        }

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] /= total;
            }
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
