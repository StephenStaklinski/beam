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

        // number of sites
        patternCount = dataInput.get().getPatternCount();
        Log.warning.println("There are " + patternCount + " unique site patterns.");

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
        }

        // Setup beagle with the defined specs of the models/data above
        setupBeagle();

        // Initialize the nodes array with the nodes from the input tree
        Node [] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = getTaxonIndex(nodes[i].getID(), dataInput.get());  
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

        // DEBUGGING
        System.out.println("calculateLogP() called");

        if(useOrigin) {
            Double originHeight = origin.getValue();
            if (treeInput.get().getRoot().getHeight() >= originHeight) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[1][m_nNodeCount];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
                operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
        }

        // beagle scaling setup
        recomputeScaleFactors = false;
        if (everUnderflowed) {
            useScaleFactors = true;
            if (rescalingCountInner < RESCALE_TIMES) {
                recomputeScaleFactors = true;
                hasDirt = Tree.IS_FILTHY;
            }
            rescalingCountInner++;
            rescalingCount++;
            if (rescalingCount > RESCALE_FREQUENCY) {
                rescalingCount = 0;
                rescalingCountInner = 0;
            }
        }

        branchUpdateCount = 0;
        operationCount = 0;

        // Traverse the tree to update any necessary transition matrices
        final Node root = treeInput.get().getRoot();
        traverse(root, true);

        // Initialize variables to return logL and determing when to return if there is numerical instability requiring scaling to recompute a NaN or -Infinity logL
        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {
            // Update partial likelihood up to the root in beagle
            beagle.updatePartials(operations[0], operationCount, Beagle.NONE);

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            // Get the root frequencies for possible states
            double[] frequencies = rootFrequenciesInput.get() == null ? substitutionModel.getFrequencies() : rootFrequenciesInput.get().getFreqs();

            // More beagle scaling stuff
            int cumulateScaleBufferIndex = Beagle.NONE;
            if (useScaleFactors) {
                if (recomputeScaleFactors) {
                    scaleBufferHelper.flipOffset(internalNodeCount);
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                    beagle.resetScaleFactors(cumulateScaleBufferIndex);
                    beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
                } else {
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                }
            } else if (useAutoScaling) {
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, Beagle.NONE);
            }
            
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
             * 
             * If scaling is required, we take advantage of the fact that the log scale factors are simply accumulated 
             * across nodes and added to the logL. We obtain the existing sum of log scale factors at the root, add the 
             * log scale factor for the origin-to-root branch, and then combine the total log scale factors sum back into
             * the calculated logL.
             */
            if (useOrigin && root.getHeight() != origin.getValue()) {

                // get the beagle calculated root partials
                double[] rootPartials = new double[patternCount * m_nStateCount];
                beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);

                // get the root node transition matrix, normally ignored but computed based on the height from root to origin
                double[] rootTransitionMatrix = new double[matrixDimensions];
                int rootNodeNum = root.getNr();

                double br = branchRateModel.getRateForBranch(root);
                substitutionModel.getTransitionProbabilities(root, origin.getValue(), root.getHeight(), br, probabilities);
                System.arraycopy(probabilities, 0, rootTransitionMatrix,  0, matrixDimensions);

                // calculate the origin partials
                calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

                // scale origin partials if scaling is on
                double totalScaleFactorsSum = 0.0;
                if (useScaleFactors) {
                    double[] originScaleFactors = new double[patternCount];
                    int u = 0;
                    for (int i = 0; i < patternCount; i++) {
                        double scaleFactor = 0.0;
                        int v = u;

                        for (int j = 0; j < m_nStateCount; j++) {
                            if (originPartials[v] > scaleFactor) {
                                scaleFactor = originPartials[v];
                            }
                            v++;
                        }
                        v += (patternCount - 1) * m_nStateCount;

                        if (scaleFactor < scalingThreshold) {
                            v = u;

                            for (int j = 0; j < m_nStateCount; j++) {
                                originPartials[v] /= scaleFactor; // inplace modification of originPartials to results in scaled form only
                                v++;
                            }
                            v += (patternCount - 1) * m_nStateCount;

                            originScaleFactors[i] = Math.log(scaleFactor);

                        } else {
                            originScaleFactors[i] = 0.0;
                        }
                        u += m_nStateCount;
                    }

                    // get the root cumulative scale factors in a round about way singe beagle.getLogScaleFactors() is not working properly
                    double[] sumLogLikelihoodsNoScaling = new double[1];
                    beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{Beagle.NONE}, 1, sumLogLikelihoodsNoScaling);
                    double[] sumLogLikelihoodsWithScaling = new double[1];
                    beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoodsWithScaling);
                    
                    // subtrace scaling - no scaling to get the scaling factors since the algorithm sums the cumulative scale factors to the log likelihood
                    // this is not setup for multiple rate categories since the sumLogLikelihoods is by default only 1 value
                    double rootScaleFactors = sumLogLikelihoodsWithScaling[0] - sumLogLikelihoodsNoScaling[0];

                    // combine the root scale factors already summed across patterns with the origin scale factors sum across patterns
                    double originScaleFactorsSum = 0.0;
                    for (int i = 0; i < patternCount; i++) {
                        originScaleFactorsSum += originScaleFactors[i];
                    }
                    totalScaleFactorsSum = originScaleFactorsSum + rootScaleFactors;
                }

                // replace the root partials with the origin partials in beagle to allow for the final likelihood calculation in in beagle
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), originPartials);

                // calculate the likelihood with the new partials
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{Beagle.NONE}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0] + totalScaleFactorsSum;
            
                // restore the original root partials in case the step is rejected or rescaling is required
                // this is also necessary to get the correct partials for sampling the tissue state at the root node
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), rootPartials);

                // save the root transition matrix for sampling the tissue state at the root node
                beagle.setTransitionMatrix(rootNodeNum, rootTransitionMatrix, 1);
            }
            else {
                // For no origin input, the logL simply comes from the root partials and frequencies already set in beagle
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0];
            }

            // If the logL has numerical instability, then generally repeat calculation with scaling on
            if (Double.isNaN(logL) || Double.isInfinite(logL)) {

                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt) {
                    // we have had a potential under/over flow so attempt a rescaling                	
                	useScaleFactors = true;
                    recomputeScaleFactors = true;
                    branchUpdateCount = 0;
                    operationCount = 0;
                    // Traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(root, false);
                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling, or always rescaling, so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        logP = logL;

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
            matrixUpdateIndices[0][branchUpdateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRate, probabilities);
            System.arraycopy(probabilities, 0, matrices,  0, matrixDimensions);
            beagle.setTransitionMatrix(matrixBufferHelper.getOffsetIndex(nodeNum), matrices, 1);

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

                // Set up scaling for beagle if necessary
                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;
                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);
                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);
                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;
                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }
                } else {
                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                // specify operations for beagle to perform later in updatePartials()
                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

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
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
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

        // Only store scale factors when actually used
        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
        }
        // Always store scaling and underflow flags to ensure they are restored if a rejected state turns on scaling
        storedEverUnderflowed = everUnderflowed;
        storedUseScaleFactors = useScaleFactors;

        // store origin partials
        System.arraycopy(originPartials, 0, storedOriginPartials, 0, originPartials.length);

        // Store logP and reset isDirty to false
        super.store();

        // store branch lengths
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        // DEBUGGING
        System.out.println("Original logP: " + logP);
        System.out.println("Original branch lengths: " + Arrays.toString(m_branchLengths));
        System.out.println("Original origin partials: " + Arrays.toString(originPartials));
        System.out.println("Original scale buffer indices: " + Arrays.toString(scaleBufferIndices));
        System.out.println("Original ever underflowed: " + everUnderflowed);
        System.out.println("Original use scale factors: " + useScaleFactors);
    }

    /**
     * Restores the state that was stored.
     */
    @Override
    public void restore() {

        // DEBUGGING
        System.out.println("New logP: " + logP);
        System.out.println("New branch lengths: " + Arrays.toString(m_branchLengths));
        System.out.println("New origin partials: " + Arrays.toString(originPartials));
        System.out.println("New scale buffer indices: " + Arrays.toString(scaleBufferIndices));
        System.out.println("New ever underflowed: " + everUnderflowed);
        System.out.println("New use scale factors: " + useScaleFactors);
        
        partialBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
        }

        // Restore flags for scaling and underflow
        everUnderflowed = storedEverUnderflowed;
        useScaleFactors = storedUseScaleFactors;

        // restore origin partials
        double[] tmp3 = storedOriginPartials;
        storedOriginPartials = originPartials;
        originPartials = tmp3;

        // Restore logP and reset isDirty to false
        logP = storedLogP;

        // Reset isDirty to false
        super.restore(); 

        // restore branch lengths
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;

        // DEBUGGING
        System.out.println("Restored logP: " + logP);
        System.out.println("Restored branch lengths: " + Arrays.toString(m_branchLengths));
        System.out.println("Restored origin partials: " + Arrays.toString(originPartials));
        System.out.println("Restored scale buffer indices: " + Arrays.toString(scaleBufferIndices));
        System.out.println("Restored ever underflowed: " + everUnderflowed);
        System.out.println("Restored use scale factors: " + useScaleFactors);

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

        for (i = 0; i < patternCount; i++) {
            int code = data.getPattern(taxon, i);
            int[] statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states[i] = statesForCode[0];
            else
                states[i] = code; // Causes ambiguous states to be ignored.
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

        // Attempt to get the resource order from the System Property
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }

        // Use dynamic rescaling, meaning that calculations are rescaled only when they underflow
        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }

        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        }

        if (requiredOrder.size() > 0) {
            requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
        }
        

        String r = System.getProperty(RESCALE_FREQUENCY_PROPERTY);
        if (r != null) {
            rescalingFrequency = Integer.parseInt(r);
            if (rescalingFrequency < 1) {
                rescalingFrequency = RESCALE_FREQUENCY;
            }
        }

        if (preferenceFlags == 0 && resourceList == null) { // else determine dataset characteristics
            if (m_nStateCount == 4 && patternCount < 10000) // TODO determine good cut-off
                preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }

        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();

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
	                resourceList,
	                preferenceFlags,
	                requirementFlags
	        );
        } catch (Exception e) {
            Log.warning.println("beagle failed to be initialized, check installation.");
            System.exit(1);
        }

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;
        resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());

        if (resourceDetails != null) {
            StringBuilder sb = new StringBuilder("  Using beagle version: " + BeagleInfo.getVersion()
                    + " resource ");
            sb.append(resourceDetails.getNumber()).append(": ");
            sb.append(resourceDetails.getName()).append("\n");
            if (resourceDetails.getDescription() != null) {
                String[] description = resourceDetails.getDescription().split("\\|");
                for (String desc : description) {
                    if (desc.trim().length() > 0) {
                        sb.append("    ").append(desc.trim()).append("\n");
                    }
                }
            }
            sb.append("    with instance flags: ").append(instanceDetails.toString());
            Log.info.println(sb.toString());
        } else {
            Log.warning.println("  Error retrieving beagle resource for instance: " + instanceDetails.toString());
            System.exit(1);
        }


        Log.warning.println("  Using DYNAMIC rescaling scheme only when likelihood underflows.");

        // For dynamic rescaling, set it to only start scaling the likelihood after the first underflow/instability issue
        everUnderflowed = false;


        // Set single category weight to prevent errors
        beagle.setCategoryWeights(0, new double[]{1.0});
    }


    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid getEigenDecompositionentry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
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

    public enum PartialsRescalingScheme {
        DEFAULT("default"), // whatever our current favourite default is
        NONE("none"),       // no scaling
        DYNAMIC("dynamic"), // rescale when needed and reuse scaling factors
        ALWAYS("always"),   // rescale every node, every site, every time - slow but safe
        DELAYED("delayed"), // postpone until first underflow then switch to 'always'
        AUTO("auto");       // beagle automatic scaling - currently playing it safe with 'always'
        PartialsRescalingScheme(String text) {
            this.text = text;
        }
        public String getText() {
            return text;
        }
        private final String text;
        public static PartialsRescalingScheme parseFromString(String text) {
            for (PartialsRescalingScheme scheme : PartialsRescalingScheme.values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return DEFAULT;
        }
    }

    // Initialize origin variable
    protected RealParameter origin;
    protected boolean useOrigin = false;

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each beagle instance to. If less than the number of instances then
    // will wrap around.
    // note: to use a different device, say device 2, start beast with
    // java -Dbeagle.resource.order=2 beast.app.BeastMCMC
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";

    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;

    private static final int RESCALE_FREQUENCY = 10000;
    private static final int RESCALE_TIMES = 1;

    // set value for scaling
    private double scalingThreshold = 1.0E-100;

    int m_nStateCount;
    int m_nNodeCount;
    int matrixDimensions;
    private double [] matrices;

    
    private double [] currentCategoryRates;
    private double [] currentFreqs;
    private double [] currentCategoryWeights;

    private int[][] matrixUpdateIndices;
    private int branchUpdateCount;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationCount;

    protected BufferIndexHelper partialBufferHelper;
    public BufferIndexHelper getPartialBufferHelper() {return partialBufferHelper;}
    
    protected BufferIndexHelper matrixBufferHelper;
    public BufferIndexHelper getMatrixBufferHelper() {return matrixBufferHelper;}
    protected BufferIndexHelper scaleBufferHelper;

    protected /*final*/ int tipCount;
    protected /*final*/ int internalNodeCount;
    protected /*final*/ int patternCount;

    private int rescalingFrequency = RESCALE_FREQUENCY;
    protected boolean useScaleFactors = false;
    protected boolean storedUseScaleFactors;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private boolean storedEverUnderflowed;

    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    // an array used to transfer tip partials
    protected double[] tipPartials;

    // the beagle library instance
    protected Beagle beagle;
    
    public Beagle getBeagle() {return beagle;}

    // Declare originPartials
    protected double[] originPartials;
    protected double[] storedOriginPartials;

    // Declare branch lengths arrays
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

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

}