package beam.likelihood;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Stephen Staklinski
 **/
@Description("Contains methods to calculate the partial likelihoods by using a simplified pruning " +
                "algorithm to save on computations given the irreversible assumptions of the substitution model.")
public class IrreversibleLikelihoodCore extends BeerLikelihoodCore {

    
    public IrreversibleLikelihoodCore(int nrOfStates) {
        super(nrOfStates);
    }


	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities, int numNodes) {
        
        super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);

        numNodesNoOrigin = numNodes;

        ancestralStates = new HashSet[numNodesNoOrigin * nrOfPatterns];
        storedAncestralStates = new HashSet[numNodesNoOrigin * nrOfPatterns];
        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i] = new HashSet<>();
            storedAncestralStates[i] = new HashSet<>();
        }
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


    public void setLeafStatesSet(int leafIndex, int[] states, int missingDataState) {

        for (int i = 0; i < nrOfPatterns; i++) {
            int index = (i * nrOfPatterns) + leafIndex;
            if (states[i] != missingDataState) {
                ancestralStates[index].add(states[i]);
            }
        }
    }

    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int nodeIndex) {
        
        for (int i = 0; i < nrOfPatterns; i++) {
            int patternStart = i * nrOfPatterns;
            ancestralStates[patternStart + nodeIndex] = new HashSet<>();

            ancestralStates[patternStart + nodeIndex].addAll(ancestralStates[patternStart + childIndex1]);
            ancestralStates[patternStart + nodeIndex].addAll(ancestralStates[patternStart + childIndex2]);

            if (ancestralStates[patternStart + nodeIndex].size() >= 1) {
                ancestralStates[patternStart + nodeIndex].clear();
                ancestralStates[patternStart + nodeIndex].add(0);
            }  
        }
    }


    /**
     * Calculates partial likelihoods at a node while
     * first checking if the node has to be unedited 
     * to skip unnecessary calculations.
     */
    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex) {
        
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
            Set<Integer> possibleStates = new HashSet<>(ancestralStates[(k * nrOfPatterns) + parentIndex]);
            // Set<Integer> possibleStates = new HashSet<>();

            // if the subtree only has missing data, then the ancestral state can be anything
            if (possibleStates.isEmpty()) {
                for (int state = 0; state < nrOfStates; state++) {
                    possibleStates.add(state);
                }
            } else {
                // always add the unedited state
                possibleStates.add(0);
            }

            // Calculate partials for all states
            for (int i : possibleStates) {
                sum1 = 0.0;
                sum2 = 0.0;
                for (int j = 0; j < nrOfStates; j++) {
                    sum1 += matrices1[i * nrOfStates + j] * partials1[u + j];
                    sum2 += matrices2[i * nrOfStates + j] * partials2[u + j];
                }
                partials3[u + i] = sum1 * sum2;
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


    /**
     * Store the stored state
     */
    @Override
    public void store() {

        super.store();

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            storedAncestralStates[i] = new HashSet<>(ancestralStates[i]);
        }

    }

    /**
     * Restore current state
     */
    @Override
    public void restore() {

        super.restore();

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i] = new HashSet<>(storedAncestralStates[i]);
        }

    }


    protected int numNodesNoOrigin;

    protected Set<Integer>[] ancestralStates;
    protected Set<Integer>[] storedAncestralStates;


}
