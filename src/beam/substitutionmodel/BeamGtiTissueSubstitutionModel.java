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

        // set up rate matrix and diagonal in one pass
        double[] rowSums = new double[nrOfStates];
        int count = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
            if (i != j) {
                rateMatrix[i][j] = relativeRates[count] * freq;
                rowSums[i] += rateMatrix[i][j];
                count++;
            }
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
