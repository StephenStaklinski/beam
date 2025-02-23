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


    /*
     * Calculates partial likelihoods at the cell division origin node with a single child.
     * Since this is the start of the experiment, the origin is known to be in the unedited
     * state, so we can only calculate the partials for that state and set the other to 0
     * since they will not be used by the frequencies anyways.
     */
    public void calculatePartials(int rootIndex, int originIndex) {

        double[] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        double[] matrices1 = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        double[] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];

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
     * Calculates partial likelihoods at a node while
     * first checking if the node has to be unedited 
     * to skip unnecessary calculations.
     */
    public void calculatePartials(Node child1, Node child2, Node node, int[] uneditedStatus) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex2 = child2.getNr();
        int nodeIndex3 = node.getNr();
        
        double[] partials1 = partials[currentPartialsIndex[nodeIndex1]][nodeIndex1];
        double[] matrices1 = matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1];
        double[] partials2 = partials[currentPartialsIndex[nodeIndex2]][nodeIndex2];
        double[] matrices2 = matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2];
        double[] partials3 = partials[currentPartialsIndex[nodeIndex3]][nodeIndex3];


        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }


    /**
     * Calculates pattern log likelihoods at a node more efficiently by assuming the substitution model is irreversible
     * and the origin is known to be the first state as unedited.
     */
    public void calculateLogLikelihoods(int nodeIndex, double[] outLogLikelihoods) {
        double[] partials1 = partials[currentPartialsIndex[nodeIndex]][nodeIndex];
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {
            outLogLikelihoods[k] = Math.log(partials1[v]) + getLogScalingFactor(k);
            v += nrOfStates;
        }
    }

}
