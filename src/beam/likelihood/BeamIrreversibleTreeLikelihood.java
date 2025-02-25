package beam.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
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
        substitutionModel = (BeamMutationSubstitutionModel) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
        m_branchLengths = new double[treeInput.get().getNodeCount() + 1];
        storedBranchLengths = new double[treeInput.get().getNodeCount() + 1];
        probabilities = new double[substitutionModel.getStateCount() * substitutionModel.getStateCount()];
        likelihoodCore = new IrreversibleLikelihoodCore(treeInput.get().getNodeCount() + 1, substitutionModel.getStateCount(), dataInput.get().getPatternCount());
        setStates(treeInput.get().getRoot());
    }

    
    /**
     * set the data at the tips in the likelihood core
     */
    protected void setStates(Node node) {
        if (node.isLeaf()) {

            if (dataInput.get().getTaxonIndex(node.getID()) == -1) {
                throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            }

            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                likelihoodCore.setNodePartials(node.getNr(), 
                                                i, 
                                                dataInput.get().getDataType().getStatesForCode(dataInput.get().getPattern(dataInput.get().getTaxonIndex(node.getID()), i))[0], 
                                                substitutionModel.getMissingState());
            }
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

        // return -infinity if the tree is larger than the experiment start time specified as the origin
        if (treeInput.get().getRoot().getHeight() >= originInput.get().getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        // do the pruning algorithm
        traverse(treeInput.get().getRoot());

        // if there is numeric instability, turn on scaling and recalculate the likelihood
        if (logP == Double.NEGATIVE_INFINITY) {
            System.out.println("Turning on scaling to prevent numeric instability.");
            
            likelihoodCore.restore();
            likelihoodCore.setUseScaling(true);

            // do the pruning algorithm with scaling this time
            traverse(treeInput.get().getRoot());

            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is negative infinity after turning on scaling.");
            }
        }

        return logP;
    }


    /**
     * Tree traversal to calculate the likelihood after updating relevant components
    **/
    protected int traverse(Node node) {

        int update = (node.isDirty() | hasDirt);

        final double branchTime = node.getLength() * branchRateModelInput.get().getRateForBranch(node);

        // first update the transition probability matrix for this branch, for all nodes except the final root or origin
        if (update != 0  || branchTime != m_branchLengths[node.getNr()]) {
            
            m_branchLengths[node.getNr()] = branchTime;

            double parentHeight = node.isRoot() ? originInput.get().getValue() : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), branchRateModelInput.get().getRateForBranch(node), probabilities);

            likelihoodCore.setNodeMatrix(node.getNr(), 0, probabilities);

            update |= 1;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final int update1 = traverse(node.getLeft());
            final int update2 = traverse(node.getRight());

            if (update1 != 0 || update2 != 0) {

                likelihoodCore.calculatePartials(node.getLeft().getNr(), node.getRight().getNr(), node.getNr());

                // if we are already back to the root in the post-order traversal, then propagate the partials to the origin
                if (node.isRoot()) {

                    // Calculate the origin partials and gets the logLikelihoods in an efficient way that assumes the origin frequencies are known as the unedited state
                    likelihoodCore.calculateLogLikelihoods(node.getNr(), node.getNr() + 1, logP);
                }

                update |= (update1 | update2);
            }
        }

        return update;
    }


    @Override
    public void store() {
        storedLogP = logP;
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        likelihoodCore.store();
    }


    @Override
    public void restore() {
        logP = storedLogP;
        System.arraycopy(storedBranchLengths, 0, m_branchLengths, 0, m_branchLengths.length);
        likelihoodCore.restore();
    }


    @Override
    protected boolean requiresRecalculation() {
        hasDirt = 0;

        if (treeInput.get().somethingIsDirty() || ((CalculationNode) substitutionModel).isDirtyCalculation() || branchRateModelInput.get().isDirtyCalculation()) {
            hasDirt = 1;
            return true;
        }

        return false;
    }


    protected BeamMutationSubstitutionModel substitutionModel;
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
