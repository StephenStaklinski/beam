package beam.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

import beam.likelihood.IrreversibleLikelihoodCore;
import beam.substitutionmodel.BeamMutationSubstitutionModel;

import java.util.Arrays;

/**
 * @author Stephen Staklinski
 **/
@Description("Calculates the tree likelihood while considering the origin of the cell division process " +
                "and using a simplified pruning algorithm to save on computations given the irreversible " +
                "assumptions of the substitution model.")
public class BeamIrreversibleTreeLikelihood extends GenericTreeLikelihood {


    public Input<RealParameter> originInput = new Input<>("origin", "Start of the cell division process, usually start of the experiment.",Input.Validate.OPTIONAL);
    

    @Override
    public void initAndValidate() {
        
        data = dataInput.get();
        
        if (data.getTaxonCount() != treeInput.get().getLeafNodeCount() || originInput.get() == null || branchRateModelInput.get() == null || !(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("Invalid input: Ensure the number of nodes in the tree matches the number of sequences, the experiment origin time is specified, the branch rate model is specified, and the siteModel input is of type SiteModel.Base.");
        }

        originHeight = originInput.get().getValue();
        substitutionModel = (BeamMutationSubstitutionModel) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
        missingDataState = substitutionModel.getMissingState();
        int nrOfStates = substitutionModel.getStateCount();
        nrOfPatterns = dataInput.get().getPatternCount();
        probabilities = new double[(nrOfStates) * (nrOfStates)];

        likelihoodCore = new IrreversibleLikelihoodCore(treeInput.get().getNodeCount() + 1, nrOfStates, nrOfPatterns, missingDataState);

        setStates(treeInput.get().getRoot());
    }

        /**
     * set the data at the tips in the likelihood core
     */
    protected void setStates(Node node) {
        if (node.isLeaf()) {
            int taxonIndex = data.getTaxonIndex(node.getID());

            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            }

            int[] states = new int[nrOfPatterns];

            for (int i = 0; i < nrOfPatterns; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                states[i] = statesForCode[0];
            }

            // this just sets the known state partials and initializes the possible ancestral states feeder sets for the leaf
            likelihoodCore.setNodePartials(node.getNr(), states);
        } else {
            setStates(node.getLeft());
            setStates(node.getRight());
        }
    }


    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() {

        Node root = treeInput.get().getRoot();

        // return -infinity if the tree is larger than the experiment start time specified as the origin
        if (root.getHeight() >= originHeight) {
            return Double.NEGATIVE_INFINITY;
        }

        // calculate logP after updating relevant matrices
        traverse(root);

        // if there is numeric instability, turn on scaling and recalculate the likelihood
        if (logP == Double.NEGATIVE_INFINITY) {

            likelihoodCore.setUseScaling(true);

            // re-calculate logP with scaling
            traverse(root);

            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is still negative infinity after turning on scaling.");
            }
        }

        return logP;
    }


    /**
     * Computed the transition pre-order, then does a post-order traversal
     * calculating the partial likelihoods by a modified pruning algorithm
     * taking advantage of the irreversibility of the substitution model
     * to reduce the number of ancestral states to propagate partial
     * likelihoods for.
     */
    protected int traverse(Node node) {

        int update = (node.isDirty() | hasDirt);
        final int nodeIndex = node.getNr();

        // first update the transition probability matrix for this branch, for all nodes except the final root or origin
        if (update != Tree.IS_CLEAN) {

            double parentHeight = node.isRoot() ? originHeight : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), branchRateModelInput.get().getRateForBranch(node), probabilities);
            likelihoodCore.setNodeMatrix(nodeIndex, 0, probabilities);
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft();
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // calculate the partials at this node given it's children.
            // currently always does the calculation since children will always be dirty.
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int childIndex1 = child1.getNr();
                int childIndex2 = child2.getNr();

                // update possible ancestral states at the nodes for more efficient pruning algorithm in traverse to get the partial likelihoods
                likelihoodCore.setPossibleAncestralStates(childIndex1, childIndex2, nodeIndex);
                likelihoodCore.calculatePartials(childIndex1, childIndex2, nodeIndex);

                // if we are already back to the root in the post-order traversal, then propagate the partials to the origin
                if (node.isRoot()) {

                    // Calculate the origin partials and gets the logLikelihoods in an efficient way that assumes the origin frequencies are known as the unedited state
                    logP = likelihoodCore.calculateLogLikelihoods(nodeIndex, nodeIndex + 1);
                }

                update |= (update1 | update2);
            }
        }

        return update;
    }


    @Override
    public void store() {
        super.store();
        likelihoodCore.store();
    }


    @Override
    public void restore() {
        super.restore();
        likelihoodCore.restore();
    }


    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (((SiteModel.Base) siteModelInput.get()).isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }

        if (branchRateModelInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }

        return treeInput.get().somethingIsDirty();
    }


    protected Alignment data;
    protected Double originHeight;
    protected BeamMutationSubstitutionModel substitutionModel;
    protected int missingDataState;
    protected int nrOfPatterns;
    protected double[] probabilities;
    protected IrreversibleLikelihoodCore likelihoodCore;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;
    
}