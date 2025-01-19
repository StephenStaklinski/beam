package beam.substitutionmodel;


import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * @author Stephen Staklinski
 **/

@Description("Random substitution model implementation that simply returns uniform transition probabilities.")

public class BeamRandomTissueSubstitutionModel extends GeneralSubstitutionModel {

    /*
     * This class simply always returns uniform transition probabilities 
     * for an ancestral state random sampling control model for model
     * selection purposes between informative vs. uninformative data.
     */

     @Override
     public void initAndValidate(){
         updateMatrix = true;
         frequencies = frequenciesInput.get();
         nrOfStates = frequencies.getFreqs().length;
         rateMatrix = new double[nrOfStates][nrOfStates];
 
         // Verify the number of input rates is correct
         int nrInputRates = ratesInput.get().getDimension();
         relativeRates = new double[ratesInput.get().getDimension()];
     }


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // always returns equal transition probabilities for random sampling
        int nrOfStates = getStateCount();
        double prob = 1.0 / nrOfStates;
        Arrays.fill(matrix, prob);
    }


    /** sets up rate matrix**/
    @Override
    public void setupRateMatrix() {
        /*
         * This does nothing because the rate matrix is not
         * used in the likelihood calculations.
         */
    }
}

