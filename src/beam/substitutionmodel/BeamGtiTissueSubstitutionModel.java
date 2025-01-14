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
        int nrInputRates = ratesInput.get().getDimension();
        if (nrInputRates != nrOfStates * nrOfStates - nrOfStates) {
            throw new IllegalArgumentException(
                "The number of input rates must be equal to ((nrOfStates * nrOfStates) - nrOfStates) for all off diagonal rates, but it is " 
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

        int count = 0;
        for (int i=0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                if (i == j) {
                    // diagonal is setup later
                    rateMatrix[i][j] = 0;
                } else {
                    rateMatrix[i][j] = relativeRates[count] * freq;
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

        // normalise rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            subst += -rateMatrix[i][i] * freq;

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
