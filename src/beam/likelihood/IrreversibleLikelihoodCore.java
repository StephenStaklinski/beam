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

        ancestralStates = new HashSet[numNodesNoOrigin * nrOfPatterns];
        storedAncestralStates = new HashSet[numNodesNoOrigin * nrOfPatterns];
        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i] = new HashSet<>();
            storedAncestralStates[i] = new HashSet<>();
        }

        // Initialize state tracking arrays
        mustBeUnedited = new boolean[nodeCount][patternCount];
        singleState = new int[nodeCount][patternCount];
        // Initialize all states to -1 (all states possible)
        for (int i = 0; i < nodeCount; i++) {
            Arrays.fill(singleState[i], -1);
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
                ancestralStates[leafIndex * nrOfPatterns + i].add(states[i]);
            }
        }
    }


    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int parentIndex) {

        for (int i = 0; i < nrOfPatterns; i++) {
            ancestralStates[parentIndex * nrOfPatterns + i] = new HashSet<>();

            ancestralStates[parentIndex * nrOfPatterns + i].addAll(ancestralStates[childIndex1 * nrOfPatterns + i]);
            ancestralStates[parentIndex * nrOfPatterns + i].addAll(ancestralStates[childIndex2 * nrOfPatterns + i]);

            if (ancestralStates[parentIndex * nrOfPatterns + i].size() > 1) {
                ancestralStates[parentIndex * nrOfPatterns + i].clear();
                ancestralStates[parentIndex * nrOfPatterns + i].add(0);
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

        Arrays.fill(partials3, 0.0);

        for (int k = 0; k < nrOfPatterns; k++) {
            final int patternOffset = k * nrOfStates;
            
            // Parent states
            if (mustBeUnedited[parentIndex][k]) {
                // Only calculate for unedited state
                final int matrixOffset = 0;  // i = 0
                double sum1 = 0.0;
                double sum2 = 0.0;

                // Child 1 states
                if (mustBeUnedited[childIndex1][k]) {
                    sum1 = matrices1[matrixOffset] * partials1[patternOffset];
                } else if (singleState[childIndex1][k] == -1) {
                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[matrixOffset + j] * partials1[patternOffset + j];
                    }
                } else {
                    sum1 = matrices1[matrixOffset] * partials1[patternOffset] + 
                          matrices1[matrixOffset + singleState[childIndex1][k]] * partials1[patternOffset + singleState[childIndex1][k]];
                }

                // Child 2 states
                if (mustBeUnedited[childIndex2][k]) {
                    sum2 = matrices2[matrixOffset] * partials2[patternOffset];
                } else if (singleState[childIndex2][k] == -1) {
                    for (int j = 0; j < nrOfStates; j++) {
                        sum2 += matrices2[matrixOffset + j] * partials2[patternOffset + j];
                    }
                } else {
                    sum2 = matrices2[matrixOffset] * partials2[patternOffset] + 
                          matrices2[matrixOffset + singleState[childIndex2][k]] * partials2[patternOffset + singleState[childIndex2][k]];
                }

                partials3[patternOffset] = sum1 * sum2;
            } else if (singleState[parentIndex][k] == -1) {
                // Calculate for all states
                for (int i = 0; i < nrOfStates; i++) {
                    final int matrixOffset = i * nrOfStates;
                    double sum1 = 0.0;
                    double sum2 = 0.0;

                    // Child 1 states
                    if (mustBeUnedited[childIndex1][k]) {
                        sum1 = matrices1[matrixOffset] * partials1[patternOffset];
                    } else if (singleState[childIndex1][k] == -1) {
                        for (int j = 0; j < nrOfStates; j++) {
                            sum1 += matrices1[matrixOffset + j] * partials1[patternOffset + j];
                        }
                    } else {
                        sum1 = matrices1[matrixOffset] * partials1[patternOffset] + 
                              matrices1[matrixOffset + singleState[childIndex1][k]] * partials1[patternOffset + singleState[childIndex1][k]];
                    }

                    // Child 2 states
                    if (mustBeUnedited[childIndex2][k]) {
                        sum2 = matrices2[matrixOffset] * partials2[patternOffset];
                    } else if (singleState[childIndex2][k] == -1) {
                        for (int j = 0; j < nrOfStates; j++) {
                            sum2 += matrices2[matrixOffset + j] * partials2[patternOffset + j];
                        }
                    } else {
                        sum2 = matrices2[matrixOffset] * partials2[patternOffset] + 
                              matrices2[matrixOffset + singleState[childIndex2][k]] * partials2[patternOffset + singleState[childIndex2][k]];
                    }

                    partials3[patternOffset + i] = sum1 * sum2;
                }
            } else {
                // Calculate for unedited and single state
                for (int i : new int[]{0, singleState[parentIndex][k]}) {
                    final int matrixOffset = i * nrOfStates;
                    double sum1 = 0.0;
                    double sum2 = 0.0;

                    // Child 1 states
                    if (mustBeUnedited[childIndex1][k]) {
                        sum1 = matrices1[matrixOffset] * partials1[patternOffset];
                    } else if (singleState[childIndex1][k] == -1) {
                        for (int j = 0; j < nrOfStates; j++) {
                            sum1 += matrices1[matrixOffset + j] * partials1[patternOffset + j];
                        }
                    } else {
                        sum1 = matrices1[matrixOffset] * partials1[patternOffset] + 
                              matrices1[matrixOffset + singleState[childIndex1][k]] * partials1[patternOffset + singleState[childIndex1][k]];
                    }

                    // Child 2 states
                    if (mustBeUnedited[childIndex2][k]) {
                        sum2 = matrices2[matrixOffset] * partials2[patternOffset];
                    } else if (singleState[childIndex2][k] == -1) {
                        for (int j = 0; j < nrOfStates; j++) {
                            sum2 += matrices2[matrixOffset + j] * partials2[patternOffset + j];
                        }
                    } else {
                        sum2 = matrices2[matrixOffset] * partials2[patternOffset] + 
                              matrices2[matrixOffset + singleState[childIndex2][k]] * partials2[patternOffset + singleState[childIndex2][k]];
                    }

                    partials3[patternOffset + i] = sum1 * sum2;
                }
            }

            if (useScaling) {
                // Only scale the states that were actually calculated
                if (mustBeUnedited[parentIndex][k]) {
                    boolean[] states = new boolean[nrOfStates];
                    states[0] = true;
                    scalePartials(parentIndex, k, states);
                } else if (singleState[parentIndex][k] == -1) {
                    boolean[] states = new boolean[nrOfStates];
                    Arrays.fill(states, true);
                    scalePartials(parentIndex, k, states);
                } else {
                    boolean[] states = new boolean[nrOfStates];
                    states[0] = true;
                    states[singleState[parentIndex][k]] = true;
                    scalePartials(parentIndex, k, states);
                }
            }
        }
    }

    public Set<Integer> getPossibleStates(Set<Integer> ancestralStates) {

        Set<Integer> possibleStates = new HashSet<>(ancestralStates);

        // if the subtree only has missing data, then the ancestral state can be anything
        if (possibleStates.isEmpty()) {
            for (int state = 1; state < nrOfStates; state++) {
                possibleStates.add(state);
            }
        }

        // always add the unedited state
        possibleStates.add(0);

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
        double[] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];  // pointer to update partials at origin

        for (int k = 0; k < nrOfPatterns; k++) {
            // Calculate the partial for the first unedited state only, which is known at the origin
            double sum1 = 0.0;

            // Check root node states
            if (mustBeUnedited[rootIndex][k]) {
                sum1 = matrices1[0] * partials1[k * nrOfStates];
            } else if (singleState[rootIndex][k] == -1) {
                // All states possible
                for (int j = 0; j < nrOfStates; j++) {
                    sum1 += matrices1[j] * partials1[k * nrOfStates + j];
                }
            } else {
                // Unedited state and single state
                sum1 = matrices1[0] * partials1[k * nrOfStates] + 
                      matrices1[singleState[rootIndex][k]] * partials1[k * nrOfStates + singleState[rootIndex][k]];
            }
            partials3[k * nrOfStates] = sum1;

            if (useScaling) {
                // Only scale the states that were actually calculated
                if (mustBeUnedited[rootIndex][k]) {
                    boolean[] states = new boolean[nrOfStates];
                    states[0] = true;
                    scalePartials(originIndex, k, states);
                } else if (singleState[rootIndex][k] == -1) {
                    boolean[] states = new boolean[nrOfStates];
                    Arrays.fill(states, true);
                    scalePartials(originIndex, k, states);
                } else {
                    boolean[] states = new boolean[nrOfStates];
                    states[0] = true;
                    states[singleState[rootIndex][k]] = true;
                    scalePartials(originIndex, k, states);
                }
            }

            // Calculate log likelihoods
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


    protected void scalePartials(int nodeIndex, int patternNum, boolean[] possibleStates) {
        double scaleFactor = 0.0;
        int v = nrOfStates * patternNum;

        // Find the maximum partial value for scaling
        for (int j = 0; j < nrOfStates; j++) {
            if (possibleStates[j]) {
                scaleFactor = Math.max(scaleFactor, partials[currentPartialsIndex[nodeIndex]][nodeIndex][v + j]);
            }
        }

        // Scale partials if the scaleFactor is below the threshold
        if (scaleFactor < scalingThreshold) {
            v = nrOfStates * patternNum;
            for (int j = 0; j < nrOfStates; j++) {
                if (possibleStates[j]) {
                    partials[currentPartialsIndex[nodeIndex]][nodeIndex][v + j] /= scaleFactor;
                }
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
            storedAncestralStates[i].clear();
            storedAncestralStates[i].addAll(ancestralStates[i]);
        }
    }


    /**
     * Restore current state
     */
    @Override
    public void restore() {

        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);

        // this is mainly for when scaling is turned on immediately, otherwise the storedAncestralStates will already be setup
        if (storedAncestralStates == null) {
            storedAncestralStates = new HashSet[numNodesNoOrigin * nrOfPatterns];
            for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
                storedAncestralStates[i] = new HashSet<>();
                storedAncestralStates[i].addAll(ancestralStates[i]);
            }
        }

        for (int i = 0; i < numNodesNoOrigin * nrOfPatterns; i++) {
            ancestralStates[i].clear();
            ancestralStates[i].addAll(storedAncestralStates[i]);
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

    protected Set<Integer>[] ancestralStates;
    protected Set<Integer>[] storedAncestralStates;

    protected boolean useScaling = false;
    protected double[][][] scalingFactors;
    private double scalingThreshold = 1.0E-100;

    protected boolean[][] mustBeUnedited;
    protected int[][] singleState;
}