package beam.likelihood;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Stephen Staklinski
 **/
@Description("Contains methods to calculate the partial likelihoods by using a simplified pruning " +
                "algorithm to save on computations given the irreversible assumptions of the substitution model.")
public class IrreversibleLikelihoodCore extends LikelihoodCore {

    
    public IrreversibleLikelihoodCore() {}


	public void initialize(int nodeCount, int nrOfStates, int patternCount) {

        this.nrOfStates = nrOfStates;
        this.nrOfNodes = nodeCount;
        this.nrOfPatterns = patternCount;

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
    }


    /**
     * Initializes the partials at a node with the known states.
     */
    public void setNodePartials(int leafIndex, int[] states, int missingDataState) {
        
        // Initialize all partials to 0
        Arrays.fill(partials[currentPartialsIndex[leafIndex]][leafIndex], 0.0);

        // Set the partials as 1.0 for the known state
        for (int i = 0; i < nrOfPatterns; i++) {
            int state = states[i];
            partials[currentPartialsIndex[leafIndex]][leafIndex][i * nrOfStates + state] = 1.0;
            if (state != missingDataState) {
                ancestralStates[leafIndex * nrOfPatterns + i].add(states[i]);
            }
        }
    }


    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int nodeIndex) {
        
        for (int i = 0; i < nrOfPatterns; i++) {
            ancestralStates[nodeIndex * nrOfPatterns + i] = new HashSet<>();

            ancestralStates[nodeIndex * nrOfPatterns + i].addAll(ancestralStates[childIndex1 * nrOfPatterns + i]);
            ancestralStates[nodeIndex * nrOfPatterns + i].addAll(ancestralStates[childIndex2 * nrOfPatterns + i]);

            if (ancestralStates[nodeIndex * nrOfPatterns + i].size() > 1) {
                ancestralStates[nodeIndex * nrOfPatterns + i].clear();
                ancestralStates[nodeIndex * nrOfPatterns + i].add(0);
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
        
        final double[] partials1 = partials[currentPartialsIndex[childIndex1]][childIndex1];
        final double[] matrices1 = matrices[currentMatrixIndex[childIndex1]][childIndex1];
        final double[] partials2 = partials[currentPartialsIndex[childIndex2]][childIndex2];
        final double[] matrices2 = matrices[currentMatrixIndex[childIndex2]][childIndex2];
        double[] partials3 = partials[currentPartialsIndex[parentIndex]][parentIndex]; // pointer to update partials at parent node

        // fill partials with 0
        Arrays.fill(partials3, 0.0);

        // modified pruning algorithm
        double sum1, sum2;

        for (int k = 0; k < nrOfPatterns; k++) {
            Set<Integer> possibleStates = new HashSet<>(ancestralStates[parentIndex * nrOfPatterns + k]);

            // if the subtree only has missing data, then the ancestral state can be anything
            if (possibleStates.isEmpty()) {
                for (int state = 1; state < nrOfStates; state++) {
                    possibleStates.add(state);
                }
            }

            // always add the unedited state
            possibleStates.add(0);

            // Calculate partials for all states
            int u = k * nrOfStates;
            for (int i : possibleStates) {
                sum1 = 0.0;
                sum2 = 0.0;
                for (int j = 0; j < nrOfStates; j++) {
                    sum1 += matrices1[i * nrOfStates + j] * partials1[u + j];
                    sum2 += matrices2[i * nrOfStates + j] * partials2[u + j];
                }
                partials3[u + i] = sum1 * sum2;
            }

            if (useScaling) {
                scalePartials(parentIndex, k, possibleStates);
            }
        }
    }


    /**
     * Calculates partial likelihoods and pattern log likelihoods at the cell division origin node with a single child.
     * Since this is the start of the experiment, the origin is known to be in the unedited state, so we can only calculate
     * the partials for that state and set the other to 0 since they will not be used by the frequencies anyways.
     */
    public void calculateLogLikelihoods(int rootIndex, int originIndex, double[] outLogLikelihoods) {

        final double[] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        final double[] matrices1 = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        double[] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];  // pointer to update partials at origin

        for (int k = 0; k < nrOfPatterns; k++) {
            // Calculate the partial for the first unedited state only, which is known at the origin
            double sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                sum1 += matrices1[j] * partials1[k * nrOfStates + j];
            }
            partials3[k * nrOfStates] = sum1;

            if (useScaling) {
                scalePartials(originIndex, k, new HashSet<>(Arrays.asList(0)));
            }

            // Calculate log likelihoods
            outLogLikelihoods[k] = Math.log(partials3[k * nrOfStates]) + getLogScalingFactor(k);
        }
    }


    /**
     * Sets probability matrix for a node
     */
    @Override
	public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex], matrixIndex * matrixSize, matrixSize);
    }

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
    }


    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }


    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param nodeIndex
     */
    protected void scalePartials(int nodeIndex, int patternNum, Set<Integer> possibleStates) {
        
        double scaleFactor = 0.0;
        int v = nrOfStates * patternNum;

        // Find the maximum partial value for scaling
        for (int j : possibleStates) {
            scaleFactor = Math.max(scaleFactor, partials[currentPartialsIndex[nodeIndex]][nodeIndex][v++]);
        }

        // Scale partials if the scaleFactor is below the threshold
        if (scaleFactor < scalingThreshold) {
            v = nrOfStates * patternNum;
            for (int j : possibleStates) {
                partials[currentPartialsIndex[nodeIndex]][nodeIndex][v++] /= scaleFactor;
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
        if (useScaling) {
            for (int i = 0; i < nrOfNodes; i++) {
                logScalingFactor += scalingFactors[currentPartialsIndex[i]][i][patternIndex_];
            }
        }
        return logScalingFactor;
    }


    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfNodes][nrOfPatterns];
        }
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

        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int[] tmp2 = currentPartialsIndex;
        currentPartialsIndex = storedPartialsIndex;
        storedPartialsIndex = tmp2;

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
    public void setNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {}
    
    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {}
    
    @Override
    public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {}
    
    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {}
    
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

}
