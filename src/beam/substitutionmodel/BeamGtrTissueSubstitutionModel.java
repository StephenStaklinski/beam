package beam.substitutionmodel;

import java.io.PipedInputStream;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashMap; 

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * @author Stephen Staklinski
 **/

@Description("General GTI substitution model implementation.")

public class BeamGtrTissueSubstitutionModel extends GeneralSubstitutionModel {

    public Input<RealParameter> piInput = new Input<>("pi", "Stationary frequency of the first state", Validate.REQUIRED);

	@Override
    public void initAndValidate(){
        updateMatrix = true;
        frequencies = frequenciesInput.get();
        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];

        // Verify the number of input rates is correct
        int nrInputRates = ratesInput.get().getDimension();
        if (nrInputRates != ((nrOfStates * nrOfStates) - nrOfStates) / 2) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to ((nrOfStates * nrOfStates) - nrOfStates) / 2 for "
                + "all off diagonal rates on one side of the matrix, but it is " 
                + nrInputRates 
                + ". Check the dimension of the input rate parameters."
            );
        }

        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e.getMessage());
		}

        relativeRates = new double[ratesInput.get().getDimension()];
    }


    /** sets up rate matrix**/
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

        // setup off diagonal rates
        int count = 0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[count];
                rateMatrix[j][i] = relativeRates[count];
                count++;
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

}
