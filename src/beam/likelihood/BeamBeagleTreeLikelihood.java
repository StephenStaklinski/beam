package beam.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;
import beagle.InstanceDetails;
import beagle.ResourceDetails;

/**
 *
 * @author Stephen Staklinski
 */

@Description("Uses the beagle library to calculate the tree likelihood by modifying the BeagleTreeLikelihood" +
            "code to include the origin node branch and output a log likelihood identical to TideTree" +
            "for the case when editing occurs the entire duration of the experiment.")
public class BeamBeagleTreeLikelihood extends GenericTreeLikelihood {


    public Input<RealParameter> originInput = new Input<>("origin", "Start of the cell division process, usually start of the experiment.", Input.Validate.OPTIONAL);
    final public Input<Frequencies> rootFrequenciesInput = new Input<>("rootFrequencies", "prior state frequencies at root, optional", Input.Validate.OPTIONAL);

    @Override
    public void initAndValidate() {
        initialize();
    }

    private void initialize() {

        // get the origin input, representing the start of the experiment as a node of degree 1 leading to the root
        if (originInput.get() != null){
            origin = originInput.get();
            useOrigin = true;
        }
    
        // get site model and validate it
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());

        // get substitution model
        substitutionModel = (SubstitutionModel.Base) m_siteModel.substModelInput.get();
        if (!substitutionModel.canReturnComplexDiagonalization()) {
            throw new IllegalArgumentException("Substitution model needs to be able to return complex diagonalization in the current implementation.");
        }

        // get the branch model with the default branch model as a strict clock if not specified in the input
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
        	throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");
        }

        // get the number of possible outcome states in the substitution model (think unique barcode edits or tissue locations)
        m_nStateCount = substitutionModel.getStateCount();

        if (debugInputData) {
            System.out.println("There are " + m_nStateCount + " unique states.");
        }

        // number of sites
        patternCount = dataInput.get().getPatternCount();
        Log.warning.println("There are " + patternCount + " unique site patterns.");
        
        if (debugInputData) {
            System.out.println("There are " + patternCount + " unique site patterns.");
        }

        // number of rates for the sites
        if (m_siteModel.getCategoryRates(null).length != 1) {
            throw new IllegalArgumentException("Site categories are not supported in the current implementation.");
        }

        currentFreqs = new double[m_nStateCount];

        // setup probability vector for the correct size of transition probability matrix filling later on
        matrixDimensions = m_nStateCount * m_nStateCount;
        probabilities = new double[matrixDimensions];
        matrices = new double[matrixDimensions];

        // get the number of nodes in the tree
        m_nNodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = m_nNodeCount - tipCount;

        // initialize branch length parameters
        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];

        if (useOrigin) {
            // define where the origin partials will be stored
            originPartials = new double[patternCount * m_nStateCount];
            storedOriginPartials = new double[patternCount * m_nStateCount];
            // define the transition matrix for the branch from the origin to the root
            rootTransitionMatrix = new double[matrixDimensions];
            storedRootTransitionMatrix = new double[matrixDimensions];
        }

        // Setup beagle with the defined specs of the models/data above
        setupBeagle();

        // Initialize the nodes array with the nodes from the input tree
        Node [] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = getTaxonIndex(nodes[i].getID(), dataInput.get());  

            if (debugInputData) {
                System.out.println("Setting states for taxon: " + nodes[i].getID());
            }

            setStates(beagle, i, taxon);
        }

        // Initialize the site pattern weights, if there are any
        double[] patternWeights = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            patternWeights[i] = dataInput.get().getPatternWeight(i);
        }
        beagle.setPatternWeights(patternWeights);
    }
    

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {

        if (transitionMatrixDebug || partialsDebug || storeRestoreDebug) {
            System.out.println("calculateLogP() called");
        }

        if(useOrigin) {
            Double originHeight = origin.getValue();
            if (treeInput.get().getRoot().getHeight() >= originHeight) {

                if (transitionMatrixDebug || partialsDebug || storeRestoreDebug) {
                    System.out.println("Tree height larger than origin.");
                }

                return Double.NEGATIVE_INFINITY;
            }
        }

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[1][m_nNodeCount];
        }

        if (operations == null) {
                operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
        }

        branchUpdateCount = 0;
        operationCount = 0;

        // Traverse the tree to update any necessary transition matrices
        final Node root = treeInput.get().getRoot();
        traverse(root, true);

        double logL;

        // Update partial likelihood up to the root in beagle
        beagle.updatePartials(operations[0], operationCount, Beagle.NONE);

        int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

        // Get the root frequencies for possible states
        double[] frequencies = rootFrequenciesInput.get() == null ? substitutionModel.getFrequencies() : rootFrequenciesInput.get().getFreqs();
        
        // make sure the root frequencies are initialized properly and then updated if the frequencies change
        if (frequencies != currentFreqs) {
            beagle.setStateFrequencies(0, frequencies);
        }
        System.arraycopy(frequencies, 0, currentFreqs, 0, frequencies.length);

        double[] sumLogLikelihoods = new double[1];

        /*
        * If the origin is specified, it is included here. The origin is handled outside of beagle 
        * because beagle requires all nodes to have a degree of 2, which the origin does not, as its 
        * only child is the root. The strategy involves obtaining the beagle-calculated root partials 
        * and calculating the transition probability matrix for the branch from the origin to the root,
        * The root partial likelihoods and the transition probability matrix are then used to calculate 
        * the origin partial likelihoods. These origin partials are re-inserted into beagle along with 
        * the origin frequencies of each state to compute the sum logL internally in beagle.
        */
        if (useOrigin && root.getHeight() != origin.getValue()) {

            if (partialsDebug) {
                System.out.println("Using origin, so starting root to origin partials calculation.");
            }

            // get the beagle calculated root partials
            double[] rootPartials = new double[patternCount * m_nStateCount];
            beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);

            // get the root node transition matrix, normally ignored but computed based on the height from root to origin
            int rootNodeNum = root.getNr();

            double br = branchRateModel.getRateForBranch(root);
            substitutionModel.getTransitionProbabilities(root, origin.getValue(), root.getHeight(), br, probabilities);
            System.arraycopy(probabilities, 0, rootTransitionMatrix,  0, matrixDimensions);

            if (transitionMatrixDebug) {
                System.out.println("Root to origin transition matrix: " + Arrays.toString(probabilities));
            }

            // calculate the origin partials
            calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

            if (partialsDebug) {
                System.out.println("Root partials: " + Arrays.toString(rootPartials));
                System.out.println("Root transition matrix: " + Arrays.toString(rootTransitionMatrix));
                System.out.println("Branch rate: " + br);
                double length = origin.getValue() - root.getHeight();
                System.out.println("Branch length: " + length);
                System.out.println("Origin partials: " + Arrays.toString(originPartials));
            }

            // replace the root partials with the origin partials in beagle to allow for the final likelihood calculation in in beagle
            beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), originPartials);

            // calculate the likelihood with the new partials
            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{Beagle.NONE}, 1, sumLogLikelihoods);
            logL = sumLogLikelihoods[0];

            if (partialsDebug) {
                System.out.println("sumLogLikelihoods: " + Arrays.toString(sumLogLikelihoods));
            }
        
            // restore the original root partials in case the step is rejected or rescaling is required
            // this is also necessary to get the correct partials for sampling the tissue state at the root node
            beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), rootPartials);
        }
        else {
            // For no origin input, the logL simply comes from the root partials and frequencies already set in beagle
            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{Beagle.NONE}, 1, sumLogLikelihoods);
            logL = sumLogLikelihoods[0];
        }

        logP = logL;

        if (partialsDebug) {
            System.out.println("Returning logL: " + logL);
        }

        if (partialsDebug) {
            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                System.out.println("Likelihood calculation has underflow (is NaN or -Infinity) and partials scaling is not used in the current implementation, so returning -Infinity.");
            }
        }

        return logL;
    }


    /**
     * Traverse the tree to update transition probability matrices and subsequently calculate partial likelihoods for each node.
     *
     * @param node           node
     * @param flip           flip
     * @return boolean
     */
    private int traverse(Node node, boolean flip) {

        int nodeNum = node.getNr();

        // Decide if this node needs to be updated
        int update = (node.isDirty() | hasDirt);

        if (transitionMatrixDebug || partialsDebug) {
            System.out.println("Traversing node: " + node.getNr() + " with update: " + update);
        }

        // Get the clock rate for the branch
        final double branchRate = branchRateModel.getRateForBranch(node);

        /* Calculate the branch length in number of substitutions to store it, where it is dependent on 
        the clock rate to convert realTimeLength * clockRate = numSubstitutionsLength */ 
        final double branchTime = node.getLength() * branchRate;
        if (branchTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchTime);
        }

        // Update if its not the root and the node is dirty or the branch length has changed
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {

            // Store the current node branch length in case it was changed, causing the update
            m_branchLengths[nodeNum] = branchTime;

            // Flip places a flag to calculate these values later in beagle updatePartials()
            if (flip) {
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // Set which matrix to update
            final int updateCount = branchUpdateCount;
            matrixUpdateIndices[0][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRate, probabilities);
            System.arraycopy(probabilities, 0, matrices,  0, matrixDimensions);
            int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            beagle.setTransitionMatrix(matrixIndex, matrices, 1);

            if (transitionMatrixDebug) {
                System.out.println("Node: " + node.getNr() + " Branch Rate: " + branchRate + " Branch Time: " + branchTime);
                System.out.println("Node: " + node.getNr() + " Transition Matrix: " + Arrays.toString(probabilities));
            }

            branchUpdateCount++;

            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes to enforce post-order traversal
            Node child1 = node.getLeft();
            final int update1 = traverse(child1, flip);

            Node child2 = node.getRight();
            final int update2 = traverse(child2, flip);

            // If either child node was dirty, then update the parent node
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int x = operationCount * Beagle.OPERATION_TUPLE_SIZE;

                // Flip places a flag to calculate these values later in beagle updatePartials()
                if (flip) {
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[0];
                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                // Not using scaleFactors
                operations[x + 1] = Beagle.NONE;
                operations[x + 2] = Beagle.NONE;

                // specify operations for beagle to perform later in updatePartials()
                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                if (partialsDebug) {
                    System.out.println("Updating partials for node: " + nodeNum);
                    System.out.println("Operations: " + Arrays.toString(Arrays.copyOfRange(operations, x, x + Beagle.OPERATION_TUPLE_SIZE)));
                }

                operationCount++;

                update |= (update1 | update2);
            }
        }
        return update;
    }

    /*
     * Calculate partials for the origin node of degree 1 that goes from the start of the experiment to the root
     */
    protected double[] calculateOriginPartials(double[] partials1, double[] matrices1, double[] partials3) {
        double sum1;
        int u = 0;
        int v = 0;

        for (int k = 0; k < patternCount; k++) {
            int w = 0;
            for (int i = 0; i < m_nStateCount; i++) {
                sum1 = 0.0;
                for (int j = 0; j < m_nStateCount; j++) {
                    sum1 += matrices1[w] * partials1[v + j];
                    w++;
                }
                partials3[u] = sum1;
                u++;
            }
            v += m_nStateCount;
        }

        return partials3;
    }


    @Override
    protected boolean requiresRecalculation() {

        hasDirt = Tree.IS_CLEAN;

        if (substitutionModel instanceof CalculationNode) {
            if (((CalculationNode) substitutionModel).isDirtyCalculation()) {
                hasDirt = Tree.IS_DIRTY;
                return true;
            }
        }
        
        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void store() {

        partialBufferHelper.storeState();
        matrixBufferHelper.storeState();

        if (useOrigin) {
            // store origin partials
            System.arraycopy(originPartials, 0, storedOriginPartials, 0, originPartials.length);

            // Store root to origin branch transition matrix
            System.arraycopy(rootTransitionMatrix, 0, storedRootTransitionMatrix, 0, rootTransitionMatrix.length);

        }
        // Store logP and reset isDirty to false
        super.store();

        // store branch lengths
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        if (storeRestoreDebug) {
            System.out.println("Original logP: " + logP);
            System.out.println("Original branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials);
            System.out.println("Original root partials: " + Arrays.toString(rootPartials));
            System.out.println("Original origin partials: " + Arrays.toString(originPartials));
            System.out.println("Original root matrix: " + Arrays.toString(rootTransitionMatrix));
        }
    }

    /**
     * Restores the state that was stored.
     */
    @Override
    public void restore() {

        if (storeRestoreDebug) {
            System.out.println("New logP: " + logP);
            System.out.println("New branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials);
            System.out.println("New root partials: " + Arrays.toString(rootPartials));
            System.out.println("New origin partials: " + Arrays.toString(originPartials));
            System.out.println("New root matrix: " + Arrays.toString(rootTransitionMatrix));
        }
        
        partialBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useOrigin) {
            // restore origin partials
            double[] tmp3 = storedOriginPartials;
            storedOriginPartials = originPartials;
            originPartials = tmp3;

            // Restore root to origin branch transition matrix
            System.arraycopy(storedRootTransitionMatrix, 0, rootTransitionMatrix, 0, rootTransitionMatrix.length);
        }


        // Restore logP and reset isDirty to false
        logP = storedLogP;

        // Reset isDirty to false
        super.restore(); 

        // restore branch lengths
        double[] tmp = storedBranchLengths;
        storedBranchLengths = m_branchLengths;
        m_branchLengths = tmp;

        if (storeRestoreDebug) {
            System.out.println("Restored logP: " + logP);
            System.out.println("Restored branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials2 = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials2);
            System.out.println("Restored root partials: " + Arrays.toString(rootPartials2));
            System.out.println("Restored origin partials: " + Arrays.toString(originPartials));
            System.out.println("Restored root matrix: " + Arrays.toString(rootTransitionMatrix));
        }

    }


    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon         the taxon
     */
    protected final void setStates(Beagle beagle,
                                   int nodeIndex, int taxon) {
        Alignment data = dataInput.get();
        int i;

        int[] states = new int[patternCount];

        if (debugInputData) {
            System.out.println("Number of site patterns to set: " + patternCount);
        }

        for (i = 0; i < patternCount; i++) {
            int code = data.getPattern(taxon, i);
            int[] statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states[i] = statesForCode[0];
            else
                states[i] = code; // Causes ambiguous states to be ignored.
        
            if (debugInputData) {
                System.out.println("For pattern " + i + " the code is " + code + " and the stateForCode is " + Arrays.toString(statesForCode) + " so the state is " + states[i]);
            }

            }

        if (debugInputData) {
            System.out.println("Setting states for node: " + nodeIndex + " with states: " + Arrays.toString(states));
        }

        beagle.setTipStates(nodeIndex, states);
    }


    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {    	
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}


    /*
     * Basic initialization of a beagle instance taken from BeagleTreeLikelihood
     */
    private void setupBeagle() {
        // one partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(m_nNodeCount, tipCount);
        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(m_nNodeCount, 0);
        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
        scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);


        // Choose CPU resource specs
        List<Integer> preferredOrder = new ArrayList<>();
        String r = System.getProperty("beagle.preferred.flags");
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    preferredOrder.add(n);
                } catch (NumberFormatException nfe) {
                    Log.warning.println("Invalid entry '" + part + "' in " + "beagle.preferred.flags");
                }
            }
        }
        long preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        long requirementFlags = BeagleFlag.EIGEN_COMPLEX.getMask();

        instanceCount++;

        try {
            beagle = BeagleFactory.loadBeagleInstance(
	                tipCount,
	                partialBufferHelper.getBufferCount(),
	                tipCount,
	                m_nStateCount,
	                patternCount,
	                2, // EigenBufferHelper not used, so this is fixed
	                matrixBufferHelper.getBufferCount(),
	                1, // One site category only in the current implementation
	                scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
	                null,
	                preferenceFlags,
	                requirementFlags
	        );
        } catch (Exception e) {
            Log.warning.println("beagle failed to be initialized, check installation.");
            System.exit(1);
        }

        // Set single category weight to prevent errors
        beagle.setCategoryWeights(0, new double[]{1.0});
    }


    public class BufferIndexHelper {
        /**
         * @param maxIndexValue the number of possible input values for the index
         * @param minIndexValue the minimum index value to have the mirrored buffers
         */
        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;

            offsetCount = maxIndexValue - minIndexValue;
            indexOffsets = new int[offsetCount];
            storedIndexOffsets = new int[offsetCount];
        }
        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }
        void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            } // else do nothing
        }
        public int getOffsetIndex(int i) {
            if (i < minIndexValue) {
                return i;
            }
            return indexOffsets[i - minIndexValue] + i;
        }
        void getIndices(int[] outIndices) {
            for (int i = 0; i < maxIndexValue; i++) {
                outIndices[i] = getOffsetIndex(i);
            }
        }
        void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

        }
        void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
        }
        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;
        private int[] indexOffsets;
        private int[] storedIndexOffsets;
    }

    // Initialize origin variable
    protected RealParameter origin;
    protected boolean useOrigin = false;

    private static int instanceCount = 0;

    // set value for scaling
    private double scalingThreshold = 1.0E-100;

    int m_nStateCount;
    int m_nNodeCount;
    int matrixDimensions;
    private double [] matrices;

    private double [] currentFreqs;

    private int[][] matrixUpdateIndices;

    private int branchUpdateCount;

    private int[][] operations;
    private int operationCount;

    protected BufferIndexHelper partialBufferHelper;
    protected BufferIndexHelper matrixBufferHelper;
    protected BufferIndexHelper scaleBufferHelper;

    protected int tipCount;
    protected int internalNodeCount;
    protected int patternCount;

    protected Beagle beagle;

    protected double[] originPartials;
    protected double[] storedOriginPartials;

    protected double[] rootTransitionMatrix;
    protected double[] storedRootTransitionMatrix;

    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;


    protected double[] patternLogLikelihoods;
    protected double[] probabilities;

    // Various debug flags
    private boolean storeRestoreDebug = false;
    private boolean transitionMatrixDebug = false;
    private boolean partialsDebug = false;
    private boolean debugInputData = false;

}