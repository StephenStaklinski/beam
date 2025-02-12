package beam.substitutionmodel;

import java.util.HashMap; 

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;

/**
 * @author Stephen Staklinski
 **/

@Description("Substitution model that allows all seeding topologies, fixing each of primary only, met to met only, and reseeding to their own free parameter.")

public class BeamPrimaryMetReseedingAllTissueSubstitutionModel extends ComplexSubstitutionModel {

    public Input<RealParameter> piInput = new Input<>("pi", "Stationary frequency of the first state", Validate.REQUIRED);

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

        /*
         * Sets up the rate matrix to be expected 1 substitution per unit time 
         * with respect to the stationary frequencies that are paramaterized by
         * pi specifying the primary tissue frequency and 1-pi for the remaining
         * tissues in a uniform distribution.
         */

         double pi = piInput.get().getValue();
         double[] piFreqs = new double[nrOfStates];
         piFreqs[0] = pi;
         double piUniformOthers = (1 - pi) / (nrOfStates - 1);
         for (int i = 1; i < nrOfStates; i++) {
             piFreqs[i] = piUniformOthers;
         }

        // set up rate matrix
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                if (i == j) {
                    rateMatrix[i][j] = 0;
                } else if (i == 0) {
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
                rateMatrix[i][j] *= piFreqs[j];
                rateMatrix[j][i] *= piFreqs[i];
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

        // normalise rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            subst += -rateMatrix[i][i] * piFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / subst;
            }
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
