package beam.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.BeagleTreeLikelihood.PartialsRescalingScheme;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
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
        
        // basic input checks
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

        if (originInput.get() == null){
            throw new IllegalArgumentException("The expriment origin time must be input.");
        }

        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }

        if (branchRateModelInput.get() == null) {
        	throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");
        }

        // initialize variables and the likelihood core
        int nodeCount = treeInput.get().getNodeCount() + 1;
        int nrOfStates = substitutionModel.getStateCount();
        data = dataInput.get();
        
        SiteModel.Base m_siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (BeamMutationSubstitutionModel) m_siteModel.substModelInput.get();
        branchRateModel = branchRateModelInput.get();
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        probabilities = new double[(nrOfStates + 1) * (nrOfStates + 1)];
        Arrays.fill(probabilities, 1.0);

        likelihoodCore = new IrreversibleLikelihoodCore(nodeCount, nrOfStates, data.getPatternCount());

        setStates(treeInput.get().getRoot());
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
        final double nodeHeight = node.getHeight();
        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // first update the transition probability matrix for this branch, for all nodes except the final root or origin
        if (update != Tree.IS_CLEAN  || branchTime != m_branchLengths[nodeIndex]) {
            m_branchLengths[nodeIndex] = branchTime;
            Node parent = node.getParent();
            if(node.isRoot()){
                parent= new Node();
                parent.setHeight(originInput.get().getValue());
            }

            substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probabilities);
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
                likelihoodCore.calculatePartials(childIndex1, childIndex2, nodeIndex);

                // if we are already back to the root in the post-order traversal, then propagate the partials to the origin
                if (node.isRoot()) {

                    // create the origin node
                    originNode = new Node();
                    originNode.setHeight(originInput.get().getValue());
                    originNode.setNr(node.getNr() + 1);

                    // Calculate the origin partials and gets the logLikelihoods in an efficient way that assumes the origin frequencies are known as the unedited state
                    int rootIndex = node.getNr();
                    int originIndex = originNode.getNr();
                    likelihoodCore.calculateLogLikelihoods(rootIndex, originIndex, logP);
                }

                update |= (update1 | update2);
            }
        }

        return update;
    }


    /**
     * set the data at the tips in the likelihood core
     */
    protected void setStates(Node node) {
        if (node.isLeaf()) {

            if (data.getTaxonIndex(node.getID()) == -1) {
                throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            }

            int[] states = new int[data.getPatternCount()];
            for (int i = 0; i < states.length; i++) {
                states[i] = data.getDataType().getStatesForCode(data.getPattern(data.getTaxonIndex(node.getID()), i))[0];
            }

            // this just sets the known state partials and initializes the possible ancestral states feeder sets for the leaf
            likelihoodCore.setNodePartials(node.getNr(), states, substitutionModel.getMissingState());
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

        final TreeInterface tree = treeInput.get();

        // return -infinity if the tree is larger than the experiment start time specified as the origin
        if (tree.getRoot().getHeight() >= originInput.get().getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        // update site likelihoods
        traverse(tree.getRoot());

        // if there is numeric instability, turn on scaling and recalculate the likelihood
        if (logP == Double.NEGATIVE_INFINITY) {
            Log.warning.println("Turning on scaling to prevent numeric instability.");
            likelihoodCore.setUseScaling(true);
            likelihoodCore.restore();

            traverse(tree.getRoot());

            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is negative infinity after turning on scaling.");
            }
        }

        return logP;
    }


    @Override
    public void store() {

        super.store();

        storedLogP = logP;
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        likelihoodCore.store();
    }


    @Override
    public void restore() {

        super.restore();

        logP = storedLogP;
        System.arraycopy(storedBranchLengths, 0, m_branchLengths, 0, m_branchLengths.length);
        
        likelihoodCore.restore();
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

        if (branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }


    protected Alignment data;
    protected BeamMutationSubstitutionModel substitutionModel;
    protected BranchRateModel.Base branchRateModel;
    protected Node originNode;
    protected double[] probabilities;
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;
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
