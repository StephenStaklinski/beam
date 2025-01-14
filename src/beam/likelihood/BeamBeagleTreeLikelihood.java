package beam.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;
import beagle.InstanceDetails;
import beagle.ResourceDetails;
import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.evolution.likelihood.BeagleTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;

/**
 *
 * @author Stephen Staklinski
 */

@Description("Uses the BEAGLE library to calculate the tree likelihood by modifying the BeagleTreeLikelihood" +
            "code to include the origin node branch and output a log likelihood identical to TideTree" +
            "for the case when editing occurs the entire duration of the experiment.")
public class BeamBeagleTreeLikelihood extends TreeLikelihood {


    public Input<RealParameter> originInput = new Input<>("origin",
            "Start of the cell division process, usually start of the experiment.",
            Input.Validate.OPTIONAL);


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
        substitutionModel = m_siteModel.substModelInput.get();

        // get the branch model with the default branch model as a strict clock if not specified in the input
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
        	branchRateModel = new StrictClockModel();
        }

        // get the number of possible outcome states in the substitution model (think unique barcode edits or tissue locations)
        m_nStateCount = substitutionModel.getStateCount();

        // number of sites
        patternCount = dataInput.get().getPatternCount();
        Log.warning.println("There are " + patternCount + " unique site patterns.");

        // number of rates for the sites
        double[] categoryRates = m_siteModel.getCategoryRates(null);
        this.categoryCount = m_siteModel.getCategoryCount();
        currentCategoryRates = categoryRates;
        currentFreqs = new double[m_nStateCount];
        currentCategoryWeights = new double[categoryRates.length];

        // setup probability vector for the correct size of transition probability matrix filling later on
        int matrixSize = (m_nStateCount + 1) * (m_nStateCount + 1);
        probabilities = new double[matrixSize];
        matrices = new double[m_nStateCount * m_nStateCount * categoryCount];

        // get the number of nodes in the tree
        m_nNodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = m_nNodeCount - tipCount;

        // initialize branch length parameters
        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];

        // Setup BEAGLE with the defined specs of the models/data above
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

        // Initialize site rates if any
        beagle.setCategoryRates(categoryRates);
    }
    

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {

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
            branchLengths = new double[1][m_nNodeCount];
            branchUpdateCount = new int[1];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
                operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
                operationCount = new int[1];
        }

        recomputeScaleFactors = false;

        if (this.rescalingScheme == PartialsRescalingScheme.ALWAYS) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
        } else if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC && everUnderflowed) {
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
        } else if (this.rescalingScheme == PartialsRescalingScheme.DELAYED && everUnderflowed) {
            useScaleFactors = true;
            recomputeScaleFactors = true;
            hasDirt = Tree.IS_FILTHY;
            rescalingCount++;
        }

        for (int i = 0; i < 1; i++) {
            branchUpdateCount[i] = 0;
        }
        operationListCount = 0;

        operationCount[0] = 0;

        final Node root = treeInput.get().getRoot();
        traverse(root, null, true);

        if (updateSiteModel) {
            double[] categoryRates = m_siteModel.getCategoryRates(null);
            if (constantPattern != null) {
	            double [] tmp = new double [categoryRates.length - 1];

	            for (int k = 0; k < categoryRates.length; k++) {
	            	tmp[k-1] = categoryRates[k];
	            }
	            categoryRates = tmp;
	        }
            for (int i = 0; i < categoryRates.length; i++) {
            	if (categoryRates[i] != currentCategoryRates[i]) {
                    beagle.setCategoryRates(categoryRates);
                    i = categoryRates.length;
            	}
            }
            currentCategoryRates = categoryRates;
        }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            double[] categoryWeights = m_siteModel.getCategoryProportions(null);
            if (constantPattern != null) {
	            double [] tmp = new double [categoryWeights.length - 1];

	            for (int k = 0; k < categoryWeights.length; k++) {
	            	tmp[k-1] = categoryWeights[k];
	            }
	            categoryWeights = tmp;
            }
            double[] frequencies = rootFrequenciesInput.get() == null ?
                    				substitutionModel.getFrequencies() :
                    				rootFrequenciesInput.get().getFreqs();

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

            // these could be set only when they change but store/restore would need to be considered
            for (int i = 0; i < categoryWeights.length; i++) {
            	if (categoryWeights[i] != currentCategoryWeights[i]) {
                    beagle.setCategoryWeights(0, categoryWeights);
            		i = categoryWeights.length;
            	}
            }
            currentCategoryWeights = categoryWeights;
            for (int i = 0; i < frequencies.length; i++) {
            	if (frequencies[i] != currentFreqs[i]) {
                    beagle.setStateFrequencies(0, frequencies);
            		i = frequencies.length;
            	}
            }
            currentFreqs = frequencies;

            double[] sumLogLikelihoods = new double[1];

            /// replace partials with new partials at the origin based on these calculations external to BEAGLE since BEAGLE cannot handle nodes with only a single child
            if (useOrigin && root.getHeight() != origin.getValue()) {

                // get the BEAGLE calculated root partials
                double[] rootPartials = new double[patternCount * m_nStateCount * categoryCount];
                beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);

                // get the root node transition matrix, normally ignored but computed based on the height from root to origin
                double[] rootTransitionMatrix = new double[m_nStateCount * m_nStateCount * categoryCount];
                int rootNodeNum = root.getNr();

                double br = branchRateModel.getRateForBranch(root);
                for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                    double jointbr = m_siteModel.getRateForCategory(i, root) * br;
                    substitutionModel.getTransitionProbabilities(root, origin.getValue(), root.getHeight(), jointbr, probabilities);
                    System.arraycopy(probabilities, 0, rootTransitionMatrix,  m_nStateCount * m_nStateCount * i, m_nStateCount * m_nStateCount);
                }

                // define where the origin partials will be stored
                double[] originPartials = new double[patternCount * m_nStateCount * categoryCount];

                // calculate the origin partials
                calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

                // scale origin partials if scaling is on
                double originScaleFactorsSum = 0.0;
                if (useScaleFactors) {
                    double[] originScaleFactors = new double[patternCount];
                    int u = 0;
                    for (int i = 0; i < patternCount; i++) {
                        double scaleFactor = 0.0;
                        int v = u;
                        for (int k = 0; k < categoryCount; k++) {
                            for (int j = 0; j < m_nStateCount; j++) {
                                if (originPartials[v] > scaleFactor) {
                                    scaleFactor = originPartials[v];
                                }
                                v++;
                            }
                            v += (patternCount - 1) * m_nStateCount;
                        }

                        if (scaleFactor < scalingThreshold) {
                            v = u;
                            for (int k = 0; k < categoryCount; k++) {
                                for (int j = 0; j < m_nStateCount; j++) {
                                    originPartials[v] /= scaleFactor;
                                    v++;
                                }
                                v += (patternCount - 1) * m_nStateCount;
                            }
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
                    originScaleFactorsSum = originScaleFactors[0] + rootScaleFactors;
                }

                // replace the root partials with the origin partials in BEAGLE to allow for the final likelihood calculation in in BEAGLE
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), originPartials);

                // calculate the likelihood with the new partials
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{Beagle.NONE}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0] + originScaleFactorsSum;
            
                // restore the original root partials in case the step is rejected or rescaling is required
                // this is also necessary to get the correct partials for sampling the tissue state at the root node
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), rootPartials);

                // save the root transition matrix for sampling the tissue state at the root node
                beagle.setTransitionMatrix(rootNodeNum, rootTransitionMatrix, 1);

                // save the origin partials globally for sampling tissue state at the origin if frequencies are not assuming the state is known
                originPartialsGlobal = originPartials;
            }
            else {
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0];
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {

                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (rescalingScheme == PartialsRescalingScheme.DYNAMIC || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling                	
                	useScaleFactors = true;
                    recomputeScaleFactors = true;

                    for (int i = 0; i < 1; i++) {
                        branchUpdateCount[i] = 0;
                    }

                    operationCount[0] = 0;

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(root, null, false);

                    done = false; // Run through do-while loop again
                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        updateSiteModel = false;

        logP = logL;
        
        return logP;
    }


    /**
     * Traverse the tree to update transition probability matrices and subsequently calculate partial likelihoods for each node.
     *
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private int traverse(Node node, int[] operatorNumber, boolean flip) {

        int nodeNum = node.getNr();

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the transition probability matrix(ices) for this branch
        int update = (node.isDirty() | hasDirt);

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {
            m_branchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            if (flip) {
                // first flip the matrixBufferHelper
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // then set which matrix to update
            final int updateCount = branchUpdateCount[0];
            matrixUpdateIndices[0][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            if (substitutionModel.canReturnComplexDiagonalization()) {
                for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                    substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), jointBranchRate, probabilities);
                    System.arraycopy(probabilities, 0, matrices,  m_nStateCount * m_nStateCount * i, m_nStateCount * m_nStateCount);
                }
            	int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            	beagle.setTransitionMatrix(matrixIndex, matrices, 1);
            }

            branchLengths[0][updateCount] = branchTime;
            branchUpdateCount[0]++;

            update |= Tree.IS_DIRTY;

        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            Node child1 = node.getLeft();
            final int[] op1 = {-1};
            final int update1 = traverse(child1, op1, flip);

            Node child2 = node.getRight();
            final int[] op2 = {-1};
            final int update2 = traverse(child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

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

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                operationCount[operationListCount]++;

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

        for (int l = 0; l < categoryCount; l++) {
            for (int k = 0; k < patternCount; k++) {
                int w = l * matrixSize;
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
        }
        return partials3;
    }


    @Override
    protected boolean requiresRecalculation() {
    
        hasDirt = Tree.IS_CLEAN;
        
        double[] categoryRates = m_siteModel.getCategoryRates(null);
        if (constantPattern != null) {
            double [] tmp = new double [categoryRates.length - 1];

            for (int k = 0; k < categoryRates.length; k++) {
            	tmp[k-1] = categoryRates[k];
            }
            categoryRates = tmp;
        }
        for (int i = 0; i < categoryRates.length; i++) {
        	if (categoryRates[i] != currentCategoryRates[i]) {
        		updateSiteModel = true;
        		break;
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
            //m_nHasDirt = Tree.IS_FILTHY;
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

        // store origin partials
        storedOriginPartialsGlobal = originPartialsGlobal;

        // store logP
        storedLogP = logP;

        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    /**
     * Restores the state that was stored.
     */
    @Override
    public void restore() {
  		updateSiteModel = true; // this is required to upload the categoryRates to BEAGLE after the restore
        
        partialBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
        }

        // restore origin partials
        originPartialsGlobal = storedOriginPartialsGlobal;

        // restore logP
        logP = storedLogP;

        super.restore();

        // double[] tmp = m_branchLengths;
        // m_branchLengths = storedBranchLengths;
        // storedBranchLengths = tmp;
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
     * Basic initialization of a BEAGLE instance taken from BeagleTreeLikelihood
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
        if (scalingOrder == null) {
            scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
        }

        // first set the rescaling scheme to use from the parser
        rescalingScheme = PartialsRescalingScheme.DEFAULT;// = rescalingScheme;
        rescalingScheme = DEFAULT_RESCALING_SCHEME;
        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (scalingOrder.size() > 0) {
            this.rescalingScheme = PartialsRescalingScheme.parseFromString(
                    scalingOrder.get(instanceCount % scalingOrder.size()));
        }

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

        if (scaling.get().equals(Scaling.always)) {
        	this.rescalingScheme = PartialsRescalingScheme.ALWAYS;
        }
        if (scaling.get().equals(Scaling.none)) {
        	this.rescalingScheme = PartialsRescalingScheme.NONE;
        }
        
        // Define default behaviour here
        if (this.rescalingScheme == PartialsRescalingScheme.DEFAULT) {
            //if GPU: the default is^H^Hwas dynamic scaling in BEAST, now NONE
            if (resourceList != null && resourceList[0] > 1) {
                //this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                this.rescalingScheme = PartialsRescalingScheme.NONE;
            } else { // if CPU: just run as fast as possible
                //this.rescalingScheme = PartialsRescalingScheme.NONE;
                // Dynamic should run as fast as none until first underflow
                this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            }
        }

        if (this.rescalingScheme == PartialsRescalingScheme.AUTO) {
            preferenceFlags |= BeagleFlag.SCALING_AUTO.getMask();
            useAutoScaling = true;
        } else {
               // preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
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

        if (substitutionModel.canReturnComplexDiagonalization()) {
            requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        }

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
	                categoryCount,
	                scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
	                resourceList,
	                preferenceFlags,
	                requirementFlags
	        );
        } catch (Exception e) {
            Log.warning.println("BEAGLE failed to be initialized, check installation.");
            System.exit(1);
        }

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;
        resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());

        if (resourceDetails != null) {
            StringBuilder sb = new StringBuilder("  Using BEAGLE version: " + BeagleInfo.getVersion()
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
            Log.warning.println("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
            System.exit(1);
        }

                // Setup rescaling scheme for when the likelihood has numerical instability issues
                if (this.rescalingScheme == PartialsRescalingScheme.AUTO &&
                resourceDetails != null &&
                (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
            // If auto scaling in BEAGLE is not supported then do it here
            this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            Log.warning.println("  Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText());
        } else {
        	Log.warning.println("  Using rescaling scheme : " + this.rescalingScheme.getText());
        }

        // For dynamic rescaling, set it to only start scaling the likelihood after the first underflow/instability issue
        if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
            everUnderflowed = false;
        }
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


    // Initialize origin variable
    protected RealParameter origin;
    protected boolean useOrigin = false;

    List<Integer> constantPattern = null;
	public List<Integer> getConstantPattern() {return constantPattern;}
    public void setConstantPattern(List<Integer> constantPattern) {this.constantPattern = constantPattern;}


    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    // note: to use a different device, say device 2, start beast with
    // java -Dbeagle.resource.order=2 beast.app.BeastMCMC
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";

    // Which scheme to use if choice not specified (or 'default' is selected):
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

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
    private double [] matrices;

    
    private double [] currentCategoryRates;
    private double [] currentFreqs;
    private double [] currentCategoryWeights;

    private int[][] matrixUpdateIndices;
    private double[][] branchLengths;
    private int[] branchUpdateCount;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;

    protected BufferIndexHelper partialBufferHelper;
    public BufferIndexHelper getPartialBufferHelper() {return partialBufferHelper;}
    
    protected BufferIndexHelper matrixBufferHelper;
    public BufferIndexHelper getMatrixBufferHelper() {return matrixBufferHelper;}
    protected BufferIndexHelper scaleBufferHelper;

    protected /*final*/ int tipCount;
    protected /*final*/ int internalNodeCount;
    protected /*final*/ int patternCount;

    private PartialsRescalingScheme rescalingScheme = DEFAULT_RESCALING_SCHEME;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    // the number of rate categories
    protected int categoryCount;

    // an array used to transfer tip partials
    protected double[] tipPartials;

    // the BEAGLE library instance
    protected Beagle beagle;
    
    public Beagle getBeagle() {return beagle;}

    // Declare originPartialsGlobal
    protected double[] originPartialsGlobal;
    protected double[] storedOriginPartialsGlobal;

    // Flag to specify that the site model has changed
    protected boolean updateSiteModel = true;

    // Class BufferIndexHelper
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
        AUTO("auto");       // BEAGLE automatic scaling - currently playing it safe with 'always'

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

}