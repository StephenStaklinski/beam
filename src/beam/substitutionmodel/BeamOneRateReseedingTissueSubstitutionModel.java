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
        int nrInputRates = ratesInput.get().getDimension();
        if (nrInputRates != 3 ) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to 3 (one primary seeding rate, one met to met seeding rate, and one primary reseeding rate), but it is " 
                + nrInputRates
                + ". Check the dimension of the input rate parameters."
            );
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        // assume equal frequencies
        double freq = 1.0 / nrOfStates;

        double primarySeedingRate = relativeRates[0] * freq;
        double reseedingRate = relativeRates[2] * freq;
        double metToMetSeedingRate = relativeRates[1] * freq;

        // set up rate matrix
        double[] rowSums = new double[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
            if (i == j) {
                rateMatrix[i][j] = 0;
            } else if (i == 0) {
                rateMatrix[i][j] = primarySeedingRate;
            } else if (j == 0) {
                rateMatrix[i][j] = reseedingRate;
            } else {
                rateMatrix[i][j] = metToMetSeedingRate;
            }
            rowSums[i] += rateMatrix[i][j];
            }
        }

        // set up diagonal and normalize rate matrix
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = -rowSums[i];
            subst += -rateMatrix[i][i] * freq;
        }

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
            rateMatrix[i][j] /= subst;
            }
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
