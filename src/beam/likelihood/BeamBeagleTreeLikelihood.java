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
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.BeagleTreeLikelihood.PartialsRescalingScheme;
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

    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());


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

        // make sure the substitution model can return transition probabilities, since not currently setup to do within beagle
        if (!substitutionModel.canReturnComplexDiagonalization()) {
            throw new IllegalArgumentException("Substitution model must be able to return transition probabilities.");
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

        this.categoryCount = m_siteModel.getCategoryCount();

        // setup probability vector for the correct size of transition probability matrix filling later on
        matrixDimensions = m_nStateCount * m_nStateCount;
        probabilities = new double[matrixDimensions];
        matrices = new double[matrixDimensions * categoryCount];

        // get the number of nodes in the tree
        m_nNodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = m_nNodeCount - tipCount;

        // initialize branch length parameters
        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];

        eigenCount = 1;

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

        // updateSubstitutionModel = true;
        // // some subst models (e.g. WAG) never become dirty, so set up subst models right now
        // setUpSubstModel();
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
            matrixUpdateIndices = new int[eigenCount][m_nNodeCount];
            branchLengths = new double[eigenCount][m_nNodeCount];
            branchUpdateCount = new int[eigenCount];
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
                hasDirt = Tree.IS_FILTHY;// makeDirty();
//                System.err.println("Recomputing scale factors");
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


        for (int i = 0; i < eigenCount; i++) {
            branchUpdateCount[i] = 0;
        }

        operationListCount = 0;
        operationCount[0] = 0;

        // Traverse the tree to update any necessary transition matrices
        final Node root = treeInput.get().getRoot();
        traverse(root, true);

        // if (updateSubstitutionModel) {
        //     setUpSubstModel();
        // }

        // if (!substitutionModel.canReturnComplexDiagonalization()) {
        //     for (int i = 0; i < eigenCount; i++) {
        //         if (branchUpdateCount[i] > 0) {
        //             beagle.updateTransitionMatrices(
        //                     eigenBufferHelper.getOffsetIndex(i),
        //                     matrixUpdateIndices[i],
        //                     null,
        //                     null,
        //                     branchLengths[i],
        //                     branchUpdateCount[i]);
        //         }
        //     }
        // }

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            // Update partial likelihood up to the root in beagle
            beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);

            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            // Get the root frequencies for possible states
            double[] frequencies = rootFrequenciesInput.get() == null ? substitutionModel.getFrequencies() : rootFrequenciesInput.get().getFreqs();

            // Beagle scaling
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
            */
            if (useOrigin && root.getHeight() != origin.getValue()) {

                if (partialsDebug) {
                    System.out.println("Using origin, so starting root to origin partials calculation.");
                }

                // get the beagle calculated root partials
                double[] rootPartials = new double[patternCount * m_nStateCount];
                beagle.getPartials(rootIndex, Beagle.NONE, rootPartials); // DO NOT pass in cumulativeScaleBufferIndex here (it will un-scale the partials) and we need to retain the scaled partials to propagate to the origin

                // get the root node transition matrix, normally ignored but computed based on the height from root to origin
                int rootNodeNum = root.getNr();

                double br = branchRateModel.getRateForBranch(root);

                if (transitionMatrixDebug) {
                    // System.out.println("Root to origin rate matrix matrix: " + );
                    double len = origin.getValue() - root.getHeight();
                    System.out.println("Root to origin length: " + len);
                    System.out.println("Root to origin clockRate: " + br);
                    double dist = len * br;
                    System.out.println("Root to origin distance (length * clockRate): " + dist);
                }

                substitutionModel.getTransitionProbabilities(root, origin.getValue(), root.getHeight(), br, probabilities);
                System.arraycopy(probabilities, 0, rootTransitionMatrix,  0, matrixDimensions);

                if (transitionMatrixDebug) {
                    System.out.println();
                    System.out.println("Root to origin transition matrix returned to likelihood:");
                    for (int i = 0; i < m_nStateCount; i++) {
                        System.out.println(Arrays.toString(Arrays.copyOfRange(probabilities, i * m_nStateCount, (i + 1) * m_nStateCount)));
                    }
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

                // scale origin partials if scaling is on
                double originScaleFactorsSum = 0.0;
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

                    /*
                     * The scale factors are summed across all sites at the origin and added to the 
                     * log likelihood calculated by beagle, which already includes other scale factors.
                     * Since all are in log space and just addition steps, this works out.
                     */
                    for (int i = 0; i < patternCount; i++) {
                        originScaleFactorsSum += originScaleFactors[i];
                    }
                }
                
                // replace the root partials with the origin partials in beagle to allow for the final likelihood calculation in beagle
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), originPartials);

                // calculate the likelihood with the new partials
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0] + originScaleFactorsSum;


                if (partialsDebug) {
                    System.out.println("sumLogLikelihoods: " + Arrays.toString(sumLogLikelihoods));
                }
            
                // restore the original root partials in case the step is rejected or rescaling is required
                // this is also necessary to get the correct partials for sampling the tissue state at the root node
                beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), rootPartials);
            }
            else {
                // For no origin input, the logL simply comes from the root partials and frequencies already set in beagle
                beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
                logL = sumLogLikelihoods[0];
            }

            // Check for underflow in the likelihood due to numerical instability
            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (rescalingScheme == PartialsRescalingScheme.DYNAMIC || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling                	
                	useScaleFactors = true;
                    recomputeScaleFactors = true;

                    for (int i = 0; i < eigenCount; i++) {
                        branchUpdateCount[i] = 0;
                    }

                    operationCount[0] = 0;

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    traverse(root, false);

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

        // updateSubstitutionModel = false;
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
            System.out.println();
            System.out.println();
            System.out.println("Traversing node: " + node.getNr());
        }

        // Get the clock rate for the branch
        final double branchRate = branchRateModel.getRateForBranch(node);

        /* Calculate the branch length in number of substitutions to store it, where it is dependent on 
        the clock rate to convert realTimeLength * clockRate = numSubstitutionsLength */ 
        final double branchTime = node.getLength() * branchRate;

        // Update if its not the root and the node is dirty or the branch length has changed
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {

            // Store the current node branch length in case it was changed, causing the update
            m_branchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            // Flip places a flag to calculate these values later in beagle updatePartials()
            if (flip) {
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // Set which matrix to update
            final int eigenIndex = 0;
            final int updateCount = branchUpdateCount[eigenIndex];
            matrixUpdateIndices[eigenIndex][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            if (transitionMatrixDebug) {
                System.out.println("Updating transition matrix for node: " + nodeNum);
                double len = node.getParent().getHeight() - node.getHeight();
                System.out.println("Node: " + node.getNr() + " real time length: " + len);
                System.out.println("Node: " + node.getNr() + " Clock Rate: " + branchRate);
                double dist = len * branchRate;
                System.out.println("Node: " + node.getNr() + " Branch dist (clock rate * length): " + branchTime);
            }

            // if (substitutionModel.canReturnComplexDiagonalization()) {
            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRate, probabilities);
            System.arraycopy(probabilities, 0, matrices,  0, matrixDimensions);
            int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            beagle.setTransitionMatrix(matrixIndex, matrices, 1);

            if (transitionMatrixDebug) {
                System.out.println();
                System.out.println("Node: " + node.getNr() + " transition matrix returned to likelihood:");
                for (int i = 0; i < m_nStateCount; i++) {
                    System.out.println(Arrays.toString(Arrays.copyOfRange(probabilities, i * m_nStateCount, (i + 1) * m_nStateCount)));
                }
            }
            // }

            branchLengths[eigenIndex][updateCount] = branchTime;
            branchUpdateCount[eigenIndex]++;

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

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                // Flip places a flag to calculate these values later in beagle updatePartials()
                if (flip) {
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

                // specify operations for beagle to perform later in updatePartials()
                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                if (partialsDebug) {
                    System.out.println("Updating partials for node: " + nodeNum);
                    System.out.println("Operations: " + Arrays.toString(Arrays.copyOfRange(operations, x, x + Beagle.OPERATION_TUPLE_SIZE)));
                }

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
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }

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
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
//            rescalingCount = storedRescalingCount;
        }

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
    
    
    protected int getScaleBufferCount() {
        return internalNodeCount + 1;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon the taxon
     */
    protected final void setPartials(Beagle beagle,
                                     int nodeIndex, int taxon) {
        Alignment data = dataInput.get();

        double[] partials = new double[patternCount * m_nStateCount * categoryCount];

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

        	double[] tipProbabilities = data.getTipLikelihoods(taxon,i);
            if (tipProbabilities != null) {
            	for (int state = 0; state < m_nStateCount; state++) {
            		partials[v++] = tipProbabilities[state];
            	}
            }
            else {
            	int stateCount = data.getPattern(taxon, i);
                boolean[] stateSet = data.getStateSet(stateCount);
                for (int state = 0; state < m_nStateCount; state++) {
                	 partials[v++] = (stateSet[state] ? 1.0 : 0.0);                
                }
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * m_nStateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(nodeIndex, partials);
    }

    public int getPatternCount() {
        return patternCount;
    }

    // void setUpSubstModel() {
    //     // we are currently assuming a no-category model...
    //     // TODO More efficient to update only the substitution model that changed, instead of all
    // 	if (!substitutionModel.canReturnComplexDiagonalization()) {
	//         for (int i = 0; i < eigenCount; i++) {
	//             //EigenDecomposition ed = m_substitutionModel.getEigenDecomposition(i, 0);
	//             EigenDecomposition ed = substitutionModel.getEigenDecomposition(null);
	
	//             eigenBufferHelper.flipOffset(i);
	
	//             beagle.setEigenDecomposition(
	//                     eigenBufferHelper.getOffsetIndex(i),
	//                     ed.getEigenVectors(),
	//                     ed.getInverseEigenVectors(),
	//                     ed.getEigenValues());
	//         }
    // 	}
    // }


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
        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);
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
//                preferenceFlags |= BeagleFlag.SCALING_MANUAL.getMask();
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

        // if (substitutionModel.canReturnComplexDiagonalization()) {
        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        // }

        instanceCount++;

        try {
            beagle = BeagleFactory.loadBeagleInstance(
	                tipCount,
	                partialBufferHelper.getBufferCount(),
	                tipCount,
	                m_nStateCount,
	                patternCount,
	                eigenBufferHelper.getBufferCount(), // EigenBufferHelper not used, so this is fixed
	                matrixBufferHelper.getBufferCount(),
	                categoryCount,
	                scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
	                resourceList,
	                preferenceFlags,
	                requirementFlags
	        );
        } catch (Exception e) {
        	Log.warning.println("Error setting up BEAGLE. Check install.");
            System.exit(1);
        }

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
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
        } else {
        	Log.warning.println("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
            System.exit(1);
        }
        Log.warning.println("  With " + patternCount + " unique site patterns.");

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

        if (this.rescalingScheme == PartialsRescalingScheme.AUTO &&
                resourceDetails != null &&
                (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
            // If auto scaling in BEAGLE is not supported then do it here
            this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            Log.warning.println("  Auto rescaling not supported in BEAGLE, using : " + this.rescalingScheme.getText());
        } else {
        	Log.warning.println("  Using rescaling scheme : " + this.rescalingScheme.getText());
        }

        if (this.rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
            everUnderflowed = false; // If false, BEAST does not rescale until first under-/over-flow.
        }

        // NO site categories allowed in the current implementation
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


    public enum PartialsRescalingScheme {
        DEFAULT("default"), // whatever our current favourite default is
        NONE("none"),       // no scaling
        DYNAMIC("dynamic"), // rescale when needed and reuse scaling factors
        ALWAYS("always"),   // rescale every node, every site, every time - slow but safe
        DELAYED("delayed"), // postpone until first underflow then switch to 'always'
        AUTO("auto");       // BEAGLE automatic scaling - currently playing it safe with 'always'
//        KICK_ASS("kickAss"),// should be good, probably still to be discovered

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

    // Beagle related variables
    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";

    // set value for scaling
    private double scalingThreshold = 1.0E-100;

    int m_nStateCount;
    int m_nNodeCount;
    int matrixDimensions;
    private double [] matrices;

    private double [] currentFreqs;

    private int eigenCount;
    private int[][] matrixUpdateIndices;
    private double[][] branchLengths;
    private int[] branchUpdateCount;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;

    protected BufferIndexHelper partialBufferHelper;
    private BufferIndexHelper eigenBufferHelper;
    protected BufferIndexHelper matrixBufferHelper;
    protected BufferIndexHelper scaleBufferHelper;

    protected int tipCount;
    protected int internalNodeCount;
    protected int patternCount;
    protected int categoryCount;

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

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;
    protected boolean storedUpdateSubstitutionModel;


    protected double[] patternLogLikelihoods;
    protected double[] probabilities;

    // Scaling variables
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;
    private PartialsRescalingScheme rescalingScheme = DEFAULT_RESCALING_SCHEME;
    private static final int RESCALE_FREQUENCY = 10000;
    private int rescalingFrequency = RESCALE_FREQUENCY;
    private static final int RESCALE_TIMES = 1;
    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    // Various debug flags
    private boolean storeRestoreDebug = false;
    private boolean transitionMatrixDebug = false;
    private boolean partialsDebug = false;
    private boolean debugInputData = false;
}