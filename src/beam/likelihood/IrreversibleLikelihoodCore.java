package beam.likelihood;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;

/**
 * @author Stephen Staklinski
 **/
@Description("Contains methods to calculate the partial likelihoods by using a simplified pruning " +
                "algorithm to save on computations given the irreversible assumptions of the substitution model.")
public class IrreversibleLikelihoodCore extends BeerLikelihoodCore {

    
    public IrreversibleLikelihoodCore(int nrOfStates) {
        super(nrOfStates);
    }

    /**
     * Initializes the partials at a node with the known states.
     */
    public void setNodePartials(int nodeIndex, int[] states) {

        // allocate memory
        createNodePartials(nodeIndex);
        
        // Initialize all partials to 0
        Arrays.fill(partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0.0);

        // Set the partials as 1.0 for the known state
        for (int i = 0; i < nrOfPatterns; i++) {
            int state = states[i];
            partials[currentPartialsIndex[nodeIndex]][nodeIndex][i * nrOfStates + state] = 1.0;
        }
    }


    /**
     * Calculates partial likelihoods at a node while
     * first checking if the node has to be unedited 
     * to skip unnecessary calculations.
     */
    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex, int[] uneditedStatus) {
        
        final double[] partials1 = partials[currentPartialsIndex[childIndex1]][childIndex1];
        final double[] matrices1 = matrices[currentMatrixIndex[childIndex1]][childIndex1];
        final double[] partials2 = partials[currentPartialsIndex[childIndex2]][childIndex2];
        final double[] matrices2 = matrices[currentMatrixIndex[childIndex2]][childIndex2];
        double[] partials3 = partials[currentPartialsIndex[parentIndex]][parentIndex]; // pointer to update partials at parent node

        // fill partials with 0
        Arrays.fill(partials3, 0.0);

        // modified pruning algorithm
        double sum1, sum2;
        int u = 0;

        for (int k = 0; k < nrOfPatterns; k++) {
            if (uneditedStatus[k] == 1) {
                // Set partials to the unedited state
                sum1 = 0.0;
                sum2 = 0.0;
                for (int j = 0; j < nrOfStates; j++) {
                    sum1 += matrices1[j] * partials1[u + j];
                    sum2 += matrices2[j] * partials2[u + j];
                }
                partials3[u] = sum1 * sum2;
            } else {
                // Calculate partials for all states
                for (int i = 0; i < nrOfStates; i++) {
                    sum1 = 0.0;
                    sum2 = 0.0;
                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[i * nrOfStates + j] * partials1[u + j];
                        sum2 += matrices2[i * nrOfStates + j] * partials2[u + j];
                    }
                    partials3[u + i] = sum1 * sum2;
                }
            }
            u += nrOfStates;
        }


        if (useScaling) {
            scalePartials(parentIndex);
        }
    }


    /*
     * Calculates partial likelihoods at the cell division origin node with a single child.
     * Since this is the start of the experiment, the origin is known to be in the unedited
     * state, so we can only calculate the partials for that state and set the other to 0
     * since they will not be used by the frequencies anyways.
     */
    public void calculatePartials(int rootIndex, int originIndex) {

        final double[] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        final double[] matrices1 = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        double[] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];  // pointer to update partials at origin

        // Initialize all partials to 0
        Arrays.fill(partials3, 0.0);

        for (int k = 0, u = 0, v = 0; k < nrOfPatterns; k++, u += nrOfStates, v += nrOfStates) {
            // Calculate the partial for the first unedited state only, which is known at the origin
            double sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
            sum1 += matrices1[j] * partials1[v + j];
            }
            partials3[u] = sum1;
        }
        
        if (useScaling) {
            scalePartials(originIndex);
        }
    }


    /**
     * Calculates pattern log likelihoods at a node more efficiently by assuming the substitution model is irreversible
     * and the origin is known to be the first state as unedited.
     */
    public void calculateLogLikelihoods(int originIndex, double[] outLogLikelihoods) {
        double[] partials1 = partials[currentPartialsIndex[originIndex]][originIndex];
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {
            outLogLikelihoods[k] = Math.log(partials1[v]) + getLogScalingFactor(k);
            v += nrOfStates;
        }
    }

}
