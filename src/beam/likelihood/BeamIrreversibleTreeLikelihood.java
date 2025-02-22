package beam.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
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
    
    final public Input<Frequencies> rootFrequenciesInput = new Input<>("rootFrequencies", "prior state frequencies at root, optional", Input.Validate.OPTIONAL);
    
    public LikelihoodCore getLikelihoodCore() {
        return likelihoodCore;
    }

    public BeamMutationSubstitutionModel getSubstitutionModel() {return substitutionModel;}

    @Override
    public void initAndValidate() {
        
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

        int nodeCount = treeInput.get().getNodeCount();
        if (originInput.get() != null){
            origin = originInput.get();
            useOrigin = true;
            nodeCount = nodeCount + 1;
        }

        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        Alignment alignment = dataInput.get();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = (BeamMutationSubstitutionModel) m_siteModel.substModelInput.get();

        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
        	throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");
        }

        int nrOfStates = substitutionModel.getStateCount();
        nrOfPatterns = dataInput.get().getPatternCount();
        nrOfMatrices = m_siteModel.getCategoryCount();

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        likelihoodCore = new IrreversibleLikelihoodCore(nrOfStates);
        patternLogLikelihoods = new double[nrOfPatterns];
        m_fRootPartials = new double[nrOfPatterns * nrOfStates];

        probabilities = new double[(nrOfStates + 1) * (nrOfStates + 1)];
        Arrays.fill(probabilities, 1.0);

        initCore();
    }

    protected void initCore() {

        int nodeCount = treeInput.get().getNodeCount();
        if (useOrigin){
            nodeCount = nodeCount + 1;
        }

        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, // current implementation does not allow for site model categories
                false); // current implementation does not use ambiguities

        int extNodeCount;
        int intNodeCount;
        
        if (useOrigin){
            nodeCount = nodeCount -1;
            extNodeCount = nodeCount / 2 + 1;
            intNodeCount = nodeCount / 2 + 1;
        }else{
            extNodeCount = nodeCount / 2 + 1;
            intNodeCount = nodeCount / 2;
        }

        setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /**
     * This method performs a recursion down the tree and computes the partial likelihood at each node
     * and stores them in the likelihood core
     * @param node
     * @return integer indicating whether the likelihood for the parent node has to be recalculated
     */
    protected int traverse(Node node) {

        int update = 2; //currently always recalculate for every node.

        final int nodeIndex = node.getNr();
        final double nodeHeight = node.getHeight();
        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if ( (!node.isRoot() || (node.isRoot() && useOrigin)) &&  (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex]) ) {
            m_branchLengths[nodeIndex] = branchTime;
            Node parent = node.getParent();
            if(node.isRoot()){
                parent= new Node();
                parent.setHeight(origin.getValue());
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);

            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            int update1 = traverse(child1);
            final double childHeight1 = child1.getHeight();

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);
            final double childHeight2 = child2.getHeight();

            update1 = 2;

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                likelihoodCore.calculatePartials(child1, child2, node);


                if (node.isRoot()) {

                    // frequencies are known at the root and fixed in the substitution model
                    double[] rootFrequencies = substitutionModel.getFrequencies();

                    if (useOrigin){
                        Double originHeight = origin.getValue();
                        originNode = new Node();
                        originNode.setHeight(originHeight);
                        originNode.setNr(node.getNr() + 1);

                        Double rootHeight = node.getHeight();


                        likelihoodCore.calculatePartials(node, originNode);

                        final double[] proportions = m_siteModel.getCategoryProportions(originNode);
                        likelihoodCore.integratePartials(originNode.getNr(), proportions, m_fRootPartials);

                        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);

                    }else{

                        //  calculate the pattern likelihoods at the root
                        final double[] proportions = m_siteModel.getCategoryProportions(node);
                        likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                    }
                }
            }

        }
        return 2;
    }


    /**
     * set the data at the tips in the likelihood core
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            String taxon = node.getID();
            int taxonIndex = data.getTaxonIndex(taxon);
            if (taxonIndex == -1) {
                if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                    taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
                }
                if (taxonIndex == -1) {
                    throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
                }
            }

            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                states[i] = statesForCode[0];
            }
            likelihoodCore.setNodeStates(node.getNr(), states);
        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
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
        if(useOrigin) {
            Double originHeight = origin.getValue();
            if (tree.getRoot().getHeight() >= originHeight) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // if the tree is not clean, the traverse call will update all pruning calculations
        if (traverse(tree.getRoot()) != Tree.IS_CLEAN) {
            // calculate the log likelihoods at the root by summing across site patterns given the site log likelihood already calculated
            logP = 0.0;
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }

        if (logP == Double.NEGATIVE_INFINITY) {
            Log.warning.println("Turning on scaling to prevent numeric instability.");
            likelihoodCore.setUseScaling(1.01);
            likelihoodCore.unstore();

            traverse(tree.getRoot());

            // calculate the log likelihoods at the root by summing across site patterns given the site log likelihood already calculated
            logP = 0.0;
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }

            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is negative infinity after turning on scaling.");
            }
        }

        return logP;
    }


    @Override
    public void store() {
        storedLogP = logP;

        if (likelihoodCore != null) {
            likelihoodCore.store();
        }

        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }


    @Override
    public void restore() {
        logP = storedLogP;

        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }

        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }


    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


    protected BeamMutationSubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;
    protected RealParameter origin;
    protected boolean useOrigin = false;
    protected Node originNode;
    protected int nrOfMatrices;
    protected int nrOfPatterns;
    protected double[] probabilities;
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;
    protected double[] patternLogLikelihoods;
    protected double[] m_fRootPartials;
    protected IrreversibleLikelihoodCore likelihoodCore;
}
