package metastabayes.substitutionmodel;

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

@Description("Substitution model with structure indicators to setup shared rate paramaters between tissues")

public class MetastabayesGeneralTissueSubstitutionModel extends ComplexSubstitutionModel {
	
    public Input<IntegerParameter> structure = new Input<IntegerParameter>("matrixStructure",
            "integer indices to indicate structure of shared transition rate matrix parameters between diferent tissues", Validate.OPTIONAL);

    private IntegerParameter matrixStructure;

	@Override
    public void initAndValidate(){

        super.initAndValidate();

        matrixStructure = structure.get();
        nrOfStates = (int) Math.sqrt(matrixStructure.getDimension());
        rateMatrix = new double[nrOfStates][nrOfStates];

        int numRelativeRates = relativeRates.length;

        HashMap<Integer,Integer> uniqueMatrixStructure = new HashMap<Integer,Integer>();
        for (int j = 0; j < matrixStructure.getDimension(); j++) {   
            uniqueMatrixStructure.put(matrixStructure.getValue(j), j);   
        }
        int countUniqueMatrixStructure = uniqueMatrixStructure.keySet().size();

        if (numRelativeRates != countUniqueMatrixStructure) {
            throw new IllegalArgumentException("Number of relative rates must match the number of unique rates specified in the matrix structure.");
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {

        // new code to allow specification of the rate matrix structure for parameter sharing
        int n = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[matrixStructure.getValue(n)];
                n++;
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
