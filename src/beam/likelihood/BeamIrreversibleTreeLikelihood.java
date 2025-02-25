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
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount() ||
            originInput.get() == null ||
            !(siteModelInput.get() instanceof SiteModel.Base) ||
            branchRateModelInput.get() == null) {
            throw new IllegalArgumentException("Invalid input: either the number of nodes in the tree does not match the number of sequences, the experiment origin time is not specified, the siteModel input is not of type SiteModel.Base, or the branch rate model is not specified.");
        }

        // initialize variables and the likelihood core
        substitutionModel = (BeamMutationSubstitutionModel) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
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

        // first update the transition probability matrix for this branch, for all nodes except the final root or origin
        if (update != 0) {

            double parentHeight = node.isRoot() ? originInput.get().getValue() : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), branchRateModelInput.get().getRateForBranch(node), probabilities);
            likelihoodCore.setNodeMatrix(node.getNr(), 0, probabilities);
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            update = Math.max(update, traverse(node.getLeft()));
            update = Math.max(update, traverse(node.getRight()));

            if (update != 0) {

                if (update == 2) {
                    likelihoodCore.setPossibleAncestralStates(node.getLeft().getNr(), node.getRight().getNr(), node.getNr());
                }

                likelihoodCore.calculatePartials(node.getLeft().getNr(), node.getRight().getNr(), node.getNr());

                // if we are already back to the root in the post-order traversal, then propagate the partials to the origin to get the logP
                if (node.isRoot()) {
                    likelihoodCore.calculateLogLikelihoods(node.getNr(), node.getNr() + 1, logP);
                }
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
        super.store();
        likelihoodCore.restore();
    }


    @Override
    protected boolean requiresRecalculation() {
        hasDirt = 0;

        if (treeInput.get().somethingIsDirty()) {
            hasDirt = 2;
            return true;
        } else if (((CalculationNode) substitutionModel).isDirtyCalculation() || branchRateModelInput.get().isDirtyCalculation()) {
            hasDirt = 1;
            return true;
        }

        return false;
    }


    protected BeamMutationSubstitutionModel substitutionModel;
    protected double[] probabilities;
    protected IrreversibleLikelihoodCore likelihoodCore;

    /*
     * 0 = nothing to recalculate
     * 1 = recalculate transition probabilities and partials
     * 2 = recalculate possible ancestral states sets
     */
    protected int hasDirt;
    
}
