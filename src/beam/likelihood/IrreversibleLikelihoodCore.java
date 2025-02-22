package beam.likelihood;

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
     * Calculates partial likelihoods at a node, specifically designed for nodes with one child like the origin in cell division.
     */
    public void calculatePartials(Node child1, Node node) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex3 = node.getNr();

        if (states[nodeIndex1] != null) {
            calculateStatesPruning(states[nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
        } else {
            calculatePartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
            matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }


    protected double[] calculateStatesPruning(int[] stateIndex1, double[] matrices1, double[] partials3) {
        int v = 0;
        for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {
                int state1 = stateIndex1[k];
                int w = l * matrixSize;
                if (state1 < nrOfStates) {
                    for (int i = 0; i < nrOfStates; i++) {
                        partials3[v] = matrices1[w + state1];
                        v++;
                        w += nrOfStates;
                    }
                } else {
                    // single child has a gap or unknown state so set partials to 1
                    for (int j = 0; j < nrOfStates; j++) {
                        partials3[v] = 1.0;
                        v++;
                    }
                }
            }
        }
        return partials3;
    }

    // calculate partials for single child nodes
    protected double[] calculatePartialsPruning(double[] partials1, double[] matrices1, double[] partials3) {
        double sum1;
        int u = 0;
        int v = 0;
        for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {
                int w = l * matrixSize;
                for (int i = 0; i < nrOfStates; i++) {
                    sum1 = 0.0;
                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[w] * partials1[v + j];
                        w++;
                    }
                    partials3[u] = sum1;
                    u++;
                }
                v += nrOfStates;
            }
        }
        return partials3;
    }
}
