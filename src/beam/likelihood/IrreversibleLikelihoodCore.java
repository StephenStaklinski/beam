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
     * Calculates partial likelihoods at a node with two children, like a normal node in a bifurcating tree including the root.
     */
    public void calculatePartials(Node child1, Node child2, Node node) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex2 = child2.getNr();
        int nodeIndex3 = node.getNr();

            if (states[nodeIndex1] != null) {
                if (states[nodeIndex2] != null) {
                    calculateStatesStatesPruning(
                            states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                } else {
                    calculateStatesPartialsPruning(states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                }
            } else {
                if (states[nodeIndex2] != null) {
                    calculateStatesPartialsPruning(states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                } else {
                    calculatePartialsPartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                }
            }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    /*
     * Calculates partial likelihoods at the cell division origin node with a single child.
     * Since this is the start of the experiment, the origin is known to be in the unedited
     * state, so we can only calculate the partials for that state and set the other to 0
     * since they will not be used by the frequencies anyways.
     */
    public void calculatePartials(Node child1, Node node) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex3 = node.getNr();

        calculatePartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
                                    matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], 
                                    partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);

        
        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

    // calculate partials for the origin
    protected double[] calculatePartialsPruning(double[] partials1, double[] matrices1, double[] partials3) {
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
        
        return partials3;
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
