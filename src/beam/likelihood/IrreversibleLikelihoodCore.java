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

    
    public IrreversibleLikelihoodCore(int nodeCount, int nrOfStates, int patternCount) {
        this.nrOfStates = nrOfStates;
        this.nrOfNodes = nodeCount;
        this.nrOfPatterns = patternCount;

        scalingFactors = new double[2][nrOfNodes][nrOfPatterns];

        partialsSize = patternCount * nrOfStates;

        partials = new double[2][nodeCount][];

        for (int i = 0; i < nodeCount; i++) {
            this.partials[0][i] = new double[partialsSize];
            this.partials[1][i] = new double[partialsSize];
        }

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        matrixSize = nrOfStates * nrOfStates;

        matrices = new double[2][nodeCount][matrixSize];

        numNodesNoOrigin = nodeCount - 1;

        ancestralStates = new int[numNodesNoOrigin * nrOfPatterns][nrOfStates];
        storedAncestralStates = new int[numNodesNoOrigin * nrOfPatterns][nrOfStates];
        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            Arrays.fill(ancestralStates[i], -1);
            Arrays.fill(storedAncestralStates[i], -1);
        }
    }

    /**
     * Initializes the partials at a node with the known states.
     */
    public void setNodePartials(int leafIndex, int[] states, int missingDataState) {
        
        // Initialize unedited partials to 0
        for (int i = 0; i < nrOfPatterns; i++) {
            partials[currentPartialsIndex[leafIndex]][leafIndex][i * nrOfStates] = 0.0;
        }

        // Set the partials as 1.0 for the known state
        for (int i = 0; i < nrOfPatterns; i++) {
            int state = states[i];
            partials[currentPartialsIndex[leafIndex]][leafIndex][i * nrOfStates + state] = 1.0;
            if (state != missingDataState) {
                ancestralStates[leafIndex * nrOfPatterns + i][0] = states[i];
            }
        }
    }


    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int parentIndex) {
        for (int i = 0; i < nrOfPatterns; i++) {
            final int parentIndexOffset = parentIndex * nrOfPatterns + i;
            final int[] parentStates = ancestralStates[parentIndexOffset];
            final int[] child1States = ancestralStates[childIndex1 * nrOfPatterns + i];
            final int[] child2States = ancestralStates[childIndex2 * nrOfPatterns + i];
            
            // Clear parent states
            Arrays.fill(parentStates, -1);
            
            // Add states from both children
            int count = 0;
            for (int j = 0; j < nrOfStates && child1States[j] != -1; j++) {
                parentStates[count++] = child1States[j];
            }
            for (int j = 0; j < nrOfStates && child2States[j] != -1; j++) {
                boolean found = false;
                for (int k = 0; k < count; k++) {
                    if (parentStates[k] == child2States[j]) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    parentStates[count++] = child2States[j];
                }
            }

            // If more than one state is possible, only allow unedited state (0)
            if (count > 1) {
                Arrays.fill(parentStates, -1);
                parentStates[0] = 0;
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

        // Pre-calculate array indices to avoid repeated multiplication
        final int[] stateIndices = new int[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            stateIndices[i] = i * nrOfStates;
        }

        for (int k = 0; k < nrOfPatterns; k++) {
            final int u = k * nrOfStates;
            final int[] possibleStates = getPossibleStates(ancestralStates[parentIndex * nrOfPatterns + k]);
            final int[] child1States = getPossibleStates(ancestralStates[childIndex1 * nrOfPatterns + k]);
            final int[] child2States = getPossibleStates(ancestralStates[childIndex2 * nrOfPatterns + k]);

            // Calculate partials for all states
            for (int i : possibleStates) {
                if (i == -1) continue;
                double sum1 = 0.0;
                double sum2 = 0.0;
                final int iOffset = stateIndices[i];
                
                for (int j : child1States) {
                    if (j == -1) continue;
                    sum1 += matrices1[iOffset + j] * partials1[u + j];
                }
                for (int j : child2States) {
                    if (j == -1) continue;
                    sum2 += matrices2[iOffset + j] * partials2[u + j];
                }
                partials3[u + i] = sum1 * sum2;
            }

            if (useScaling) {
                scalePartials(parentIndex, k, possibleStates);
            }
        }
    }

    public int[] getPossibleStates(int[] ancestralStates) {
        int[] possibleStates = new int[nrOfStates];
        int count = 0;

        // Copy non-negative states from ancestralStates
        for (int i = 0; i < nrOfStates && ancestralStates[i] != -1; i++) {
            possibleStates[count++] = ancestralStates[i];
        }

        // If no states are set, add all states except 0
        if (count == 0) {
            for (int state = 1; state < nrOfStates; state++) {
                possibleStates[count++] = state;
            }
        }

        // Always add the unedited state if not already present
        boolean hasZero = false;
        for (int i = 0; i < count; i++) {
            if (possibleStates[i] == 0) {
                hasZero = true;
                break;
            }
        }
        if (!hasZero) {
            possibleStates[count++] = 0;
        }

        // Fill remaining slots with -1
        for (int i = count; i < nrOfStates; i++) {
            possibleStates[i] = -1;
        }

        return possibleStates;
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
            double sum1 = 0.0;
            final int[] possibleStates = getPossibleStates(ancestralStates[rootIndex * nrOfPatterns + k]);

            for (int j : possibleStates) {
                if (j == -1) continue;
                sum1 += matrices1[j] * partials1[k * nrOfStates + j];
            }
            partials3[k * nrOfStates] = sum1;

            if (useScaling) {
                scalePartials(originIndex, k, new int[]{0, -1});
            }

            if (useScaling) {
                logP += Math.log(partials3[k * nrOfStates]) + getLogScalingFactor(k);
            } else {
                logP += Math.log(partials3[k * nrOfStates]);
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
        
        // Find the maximum partial value for scaling
        double scaleFactor = 0.0;
        for (int j : possibleStates) {
            if (j == -1) continue;
            final double partial = nodePartials[v + j];
            if (partial > scaleFactor) {
                scaleFactor = partial;
            }
        }

        // Scale partials if the scaleFactor is below the threshold
        if (scaleFactor < scalingThreshold) {
            for (int j : possibleStates) {
                if (j == -1) continue;
                nodePartials[v + j] /= scaleFactor;
            }
            scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternNum] = Math.log(scaleFactor);
        } else {
            scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternNum] = 0.0;
        }
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
        
        double logScalingFactor = 0.0;
        for (int i = 0; i < nrOfNodes; i++) {
            logScalingFactor += scalingFactors[currentPartialsIndex[i]][i][patternIndex_];
        }

        return logScalingFactor;
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
            System.arraycopy(ancestralStates[i], 0, storedAncestralStates[i], 0, nrOfStates);
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
            storedAncestralStates = new int[numNodesNoOrigin * nrOfPatterns][nrOfStates];
            for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
                System.arraycopy(ancestralStates[i], 0, storedAncestralStates[i], 0, nrOfStates);
            }
        }

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            System.arraycopy(storedAncestralStates[i], 0, ancestralStates[i], 0, nrOfStates);
        }
    }


    /**
     * cleans up and deallocates arrays.
     */
    @Override
	public void finalize() throws java.lang.Throwable {
        nrOfNodes = 0;
        nrOfPatterns = 0;
        partials = null;
        currentPartialsIndex = null;
        storedPartialsIndex = null;
        matrices = null;
        currentMatrixIndex = null;
        storedMatrixIndex = null;
        scalingFactors = null;
        ancestralStates = null;
        storedAncestralStates = null;
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

    protected boolean integrateCategories;

    protected double[][][] partials;
    protected double[][][] matrices;

    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;
    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected int[][] ancestralStates;
    protected int[][] storedAncestralStates;

    protected boolean useScaling = false;
    protected double[][][] scalingFactors;
    private double scalingThreshold = 1.0E-100;

}