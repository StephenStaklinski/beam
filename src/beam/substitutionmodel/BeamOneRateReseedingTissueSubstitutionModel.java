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

        // bring in frequencies
        double [] freqs = frequencies.getFreqs();
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] *= freqs[j];
                rateMatrix[j][i] *= freqs[i];
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
            subst += -rateMatrix[i][i] * freqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / subst;
            }
        }

    }
}
