package beam.likelihood;


import beast.base.core.Description;
import beast.base.evolution.likelihood.LikelihoodCore;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Stephen Staklinski
 **/
@Description("Contains methods to calculate the partial likelihoods by using a simplified pruning " +
                "algorithm to save on computations given the irreversible assumptions of the substitution model.")
public class IrreversibleLikelihoodCore extends LikelihoodCore {


    public IrreversibleLikelihoodCore(int nodeCount, int nrOfStates, int patternCount, int missingData) {
        this.nrOfStates = nrOfStates;
        this.nrOfNodes = nodeCount;
        this.nrOfPatterns = patternCount;
        this.missingDataState = missingData;
        this.matrixSize = nrOfStates * nrOfStates;
        this.partialsSize = patternCount * nrOfStates;
        this.numNodesNoOrigin = nodeCount - 1;

        // Initialize arrays
        allStates = new int[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) allStates[i] = i;

        partials = new double[2][nodeCount][partialsSize];
        matrices = new double[2][nodeCount][matrixSize];
        scalingFactors = new double[2][nodeCount][patternCount];
        
        // Initialize all arrays to 0
        for (int i = 0; i < nodeCount; i++) {
            Arrays.fill(partials[0][i], 0.0);
            Arrays.fill(partials[1][i], 0.0);
            Arrays.fill(matrices[0][i], 0.0);
            Arrays.fill(matrices[1][i], 0.0);
            Arrays.fill(scalingFactors[0][i], 0.0);
            Arrays.fill(scalingFactors[1][i], 0.0);
        }

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];
        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        matrixSize = nrOfStates * nrOfStates;

        matrices = new double[2][nodeCount][matrixSize];

        numNodesNoOrigin = nodeCount - 1;

        ancestralStates = new int[numNodesNoOrigin * nrOfPatterns];
        storedAncestralStates = new int[numNodesNoOrigin * nrOfPatterns];
        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i] = -1;
            storedAncestralStates[i] = -1;
        }
    }

    /**
     * Initializes the partials at a node with the known states.
     */
    public void setNodePartials(int leafIndex, int[] states) {
        for (int i = 0; i < nrOfPatterns; i++) {
            // Set the known state to 1.0 for tips
            partials[currentPartialsIndex[leafIndex]][leafIndex][i * nrOfStates + states[i]] = 1.0;

            // Set the ancestral state for tips
            ancestralStates[leafIndex * nrOfPatterns + i] = (states[i] != missingDataState) ? states[i] : -2;
        }
    }


    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int parentIndex) {
        for (int i = 0; i < nrOfPatterns; i++) {
            final int parentIndexOffset = parentIndex * nrOfPatterns + i;
            final int child1State = ancestralStates[childIndex1 * nrOfPatterns + i];
            final int child2State = ancestralStates[childIndex2 * nrOfPatterns + i];

            if (child1State == 0 || child2State == 0 || (child1State > 0 && child2State > 0 && child1State != child2State)) {
                ancestralStates[parentIndexOffset] = 0;
                continue;
            } else if (child1State < 0 && child2State < 0) {
                ancestralStates[parentIndexOffset] = -1;
                continue;
            } else {
                ancestralStates[parentIndexOffset] = child1State >= 0 ? child1State : child2State;
            }
        }
    }


    /**
     * Calculates partial likelihoods at a node while
     * first checking which partials need to be calculated
     * to save on computations when 0 values should just
     * be propagated.
     */
    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex) {
        currentPartialsIndex[parentIndex] = 1 - currentPartialsIndex[parentIndex];
        
        final double[] partials1 = partials[currentPartialsIndex[childIndex1]][childIndex1];
        final double[] matrices1 = matrices[currentMatrixIndex[childIndex1]][childIndex1];
        final double[] partials2 = partials[currentPartialsIndex[childIndex2]][childIndex2];
        final double[] matrices2 = matrices[currentMatrixIndex[childIndex2]][childIndex2];
        final double[] partials3 = partials[currentPartialsIndex[parentIndex]][parentIndex];

        for (int k = 0; k < nrOfPatterns; k++) {
            final int u = k * nrOfStates;
            final int[] possibleStates = getPossibleStates(ancestralStates[parentIndex * nrOfPatterns + k]);
            final int[] child1States = getPossibleStates(ancestralStates[childIndex1 * nrOfPatterns + k]);
            final int[] child2States = getPossibleStates(ancestralStates[childIndex2 * nrOfPatterns + k]);

            // Calculate partials for all states
            for (int i : possibleStates) {
                double sum1 = 0.0;
                double sum2 = 0.0;
                final int iOffset = i * nrOfStates;
                
                for (int j : child1States) {
                    sum1 += matrices1[iOffset + j] * partials1[u + j];
                }
                for (int j : child2States) {
                    sum2 += matrices2[iOffset + j] * partials2[u + j];
                }
                partials3[u + i] = sum1 * sum2;
            }

            if (useScaling) {
                scalePartials(parentIndex, k, possibleStates);
            }
        }
    }

    
    private int[] getPossibleStates(int state) {
        return state == 0 ? uneditedState :
               state > 0 ? new int[]{0, state} :
               state == -2 ? new int[]{missingDataState} :
               allStates;
    }



    /**
     * Calculates partial likelihoods and pattern log likelihoods at the cell division origin node with a single child.
     * Since this is the start of the experiment, the origin is known to be in the unedited state, so we can only calculate
     * the partials for that state and set the other to 0 since they will not be used by the frequencies anyways.
     */
    public double calculateLogLikelihoods(int rootIndex, int originIndex) {
        double logP = 0.0;
        currentPartialsIndex[originIndex] = 1 - currentPartialsIndex[originIndex];

        final double[] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        final double[] matrices1 = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        final double[] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];

        for (int k = 0; k < nrOfPatterns; k++) {
            final int offset = k * nrOfStates;
            double sum1 = 0.0;
            final int[] possibleStates = getPossibleStates(ancestralStates[rootIndex * nrOfPatterns + k]);

            for (int j : possibleStates) {
                sum1 += matrices1[j] * partials1[offset + j];
            }
            partials3[offset] = sum1;

            if (useScaling) {
                scalePartials(originIndex, k, new int[]{0});
                logP += Math.log(partials3[offset]) + getLogScalingFactor(k);
            } else {
                logP += Math.log(partials3[offset]);
            }
        }

        return logP;
    }


    /**
     * Sets probability matrix for a node
     */
    @Override
	public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex], matrixIndex * matrixSize, matrixSize);
    }


    protected void scalePartials(int nodeIndex, int patternNum, int[] possibleStates) {
        final int v = nrOfStates * patternNum;
        final double[] nodePartials = partials[currentPartialsIndex[nodeIndex]][nodeIndex];
        final double[] nodeScalingFactors = scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex];
        
        // Find the maximum partial value for scaling
        double scaleFactor = 0.0;
        for (int j : possibleStates) {
            scaleFactor = Math.max(scaleFactor, nodePartials[v + j]);
        }

        // Scale partials if needed and set scaling factor
        if (scaleFactor < scalingThreshold) {
            for (int j : possibleStates) {
                nodePartials[v + j] /= scaleFactor;
            }
        }
        nodeScalingFactors[patternNum] = scaleFactor < scalingThreshold ? Math.log(scaleFactor) : 0.0;
    }

    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @return the log scaling factor
     */
    @Override
    public double getLogScalingFactor(int patternIndex_) {    
        return java.util.stream.IntStream.range(0, nrOfNodes)
            .mapToDouble(i -> scalingFactors[currentPartialsIndex[i]][i][patternIndex_])
            .sum();
    }


    public void setUseScaling(boolean status) {
        useScaling = status;
    }


    /**
     * Store the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            storedAncestralStates[i] = ancestralStates[i];
        }
    }


    /**
     * Restore current state
     */
    @Override
    public void restore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);

        if (storedAncestralStates == null) {
            storedAncestralStates = new int[numNodesNoOrigin * nrOfPatterns];
            for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
                storedAncestralStates[i] = ancestralStates[i];
            }
        }

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i] = storedAncestralStates[i];
        }
    }


    /**
     * cleans up and deallocates arrays.
     */
    @Override
	public void finalize() throws java.lang.Throwable {
        nrOfNodes = 0;
        nrOfPatterns = 0;
        nrOfStates = 0;
        numNodesNoOrigin = 0;
        partialsSize = 0;
        matrixSize = 0;
        partials = null;
        currentPartialsIndex = null;
        storedPartialsIndex = null;
        matrices = null;
        currentMatrixIndex = null;
        storedMatrixIndex = null;
        scalingFactors = null;
        ancestralStates = null;
        storedAncestralStates = null;
        allStates = null;
        uneditedState = null;
    }


    // Bunch of things that need to be implemented but are not used in this class

    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities) {}

    
    @Override
	public void createNodePartials(int nodeIndex) {}

    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {}
    
    @Override
    public void getNodePartials(int nodeIndex, double[] partials) {}

    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {}
    
    @Override
    public void setNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {}

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {}
    
    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {}
    
    @Override
    public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {}
    
    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {}
    
    @Override
    public void setUseScaling(double scale) {}

    @Override
    public void unstore() {}


    protected int nrOfStates;
    protected int nrOfNodes;
    protected int numNodesNoOrigin;
    protected int nrOfPatterns;
    protected int partialsSize;
    protected int matrixSize;

    protected double[][][] partials;
    protected double[][][] matrices;

    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;
    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected int[] ancestralStates;
    protected int[] storedAncestralStates;

    private int[] allStates;
    private int[] uneditedState = new int[]{0};

    protected boolean useScaling = false;
    protected double[][][] scalingFactors;
    private double scalingThreshold = 1.0E-100;

    private int missingDataState;

}