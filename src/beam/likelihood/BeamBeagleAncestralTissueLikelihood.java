package beam.likelihood;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beagle.Beagle;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import beastclassic.evolution.tree.TreeTrait;
import beastclassic.evolution.tree.TreeTraitProvider;
import beastclassic.evolution.likelihood.LeafTrait;

/**
 * Ancestral State Tree Likelihood that extends the custom Beagle likelihood with origin input
 * for the branch above the root.
 *
 * @author Stephen Staklinski
 */
@Description("Ancestral State Tree Likelihood that extends the custom Beagle likelihood with origin input" +
        " for the branch above the root.")
public class BeamBeagleAncestralTissueLikelihood extends BeamBeagleTreeLikelihood implements TreeTraitProvider {

    /** Key for states trait */
    public static final String STATES_KEY = "states";

    /** Input for tag label used in output tree posterior samples */
    public Input<String> tagInput = new Input<>("tag",
            "label used to report trait in the output tree posterior samples",
            Validate.REQUIRED);

    /** Input for leaf traits */
    public Input<List<LeafTrait>> leafTraitsInput = new Input<>("leaftrait",
            "list of leaf traits",
            new ArrayList<>());

    /** Data type for the alignment */
    protected DataType dataType;

    /** Reconstructed states for each node */
    private int[][] reconstructedStates;

    /** Stored reconstructed states for store/restore operations */
    private int[][] storedReconstructedStates;

    /** Tag label for trait reporting */
    private String tag;

    /** Flag indicating if states have been redrawn */
    private boolean areStatesRedrawn = false;

    /** Stored flag for states redrawn status */
    private boolean storedAreStatesRedrawn = false;

    /** Number of possible states */
    private int stateCount;

    /** Current state */
    private int state;

    /** States at tip nodes */
    private int[][] tipStates;

    /** Helper for tree traits */
    protected TreeTraitProvider.Helper treeTraits = new Helper();

    @Override
    public void initAndValidate() {
        validateInputData();
        super.initAndValidate();
        validateSiteModel();
        initializeDataStructures();
        setupTreeTraits();
    }

    private void validateInputData() {
        if (dataInput.get().getPatternCount() != 1) {
            throw new RuntimeException("More than one site in the input data not implemented.");
        }
    }

    private void validateSiteModel() {
        if (m_siteModel.getCategoryCount() > 1) {
            throw new RuntimeException("Multiple site rate categories not implemented.");
        }
    }

    private void initializeDataStructures() {
        tag = tagInput.get();
        Alignment data = dataInput.get();
        dataType = data.getDataType();
        stateCount = dataType.getStateCount();

        TreeInterface treeModel = treeInput.get();
        int nodeCount = treeModel.getNodeCount();

        initializeTipStates(treeModel, data);
        initializeReconstructedStates(nodeCount);
    }

    private void initializeTipStates(TreeInterface treeModel, Alignment data) {
        tipStates = new int[treeModel.getLeafNodeCount()][1];
        for (Node node : treeModel.getExternalNodes()) {
            String taxon = node.getID();
            int taxonIndex = findTaxonIndex(taxon, data);
            int[] states = tipStates[node.getNr()];
            int code = data.getPattern(taxonIndex, 0);
            states[0] = data.getDataType().getStatesForCode(code)[0];
        }
    }

    /**
     * Finds the taxon index in the alignment data.
     *
     * @param taxon The taxon name to find
     * @param data The alignment data
     * @return The index of the taxon in the alignment
     * @throws RuntimeException if the taxon is not found
     */
    private int findTaxonIndex(String taxon, Alignment data) {
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

    private void initializeReconstructedStates(int nodeCount) {
        reconstructedStates = new int[nodeCount][1];
        storedReconstructedStates = new int[nodeCount][1];
    }

    private void setupTreeTraits() {
        treeTraits.addTrait(STATES_KEY, new TreeTrait.IA() {
            @Override
            public String getTraitName() {
                return tag;
            }

            @Override
            public Intent getIntent() {
                return Intent.NODE;
            }

            @Override
            public int[] getTrait(TreeInterface tree, Node node) {
                return getStatesForNode(tree, node);
            }

            @Override
            public String getTraitString(TreeInterface tree, Node node) {
                return formattedState(getStatesForNode(tree, node), dataType);
            }
        });
    }

    @Override
    public void store() {
        super.store();
        storeReconstructedStates();
    }

    private void storeReconstructedStates() {
        System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
        storedAreStatesRedrawn = areStatesRedrawn;
    }

    @Override
    public void restore() {
        super.restore();
        restoreReconstructedStates();
    }

    private void restoreReconstructedStates() {
        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;
        areStatesRedrawn = storedAreStatesRedrawn;
    }

    /**
     * Gets the states for a given node in the tree.
     *
     * @param tree The tree to get states from
     * @param node The node to get states for
     * @return The states for the node
     */
    public int[] getStatesForNode(TreeInterface tree, Node node) {
        if (tree != treeInput.get()) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!areStatesRedrawn) {
            traverseSample(tree, tree.getRoot(), -1);
            areStatesRedrawn = true;
        }
        return reconstructedStates[node.getNr()];
    }

    @Override
    public double calculateLogP() {
        super.calculateLogP();
        areStatesRedrawn = false;
        return logP;
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    /**
     * Formats a state array into a string representation.
     *
     * @param state The state array to format
     * @param dataType The data type for state interpretation
     * @return Formatted string representation of the state
     */
    private static String formattedState(int[] state, DataType dataType) {
        StringBuilder sb = new StringBuilder();
        sb.append("\"");
        
        if (dataType instanceof UserDataType) {
            boolean first = true;
            for (int i : state) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }
                sb.append(dataType.getCode(i));
            }
        } else {
            for (int i : state) {
                sb.append(dataType.getChar(i));
            }
        }
        
        sb.append("\"");
        return sb.toString();
    }

    /**
     * Traverses the tree in pre-order to sample internal node states.
     *
     * @param tree The tree to traverse
     * @param node The current node
     * @param parentState The state of the parent node
     */
    public void traverseSample(TreeInterface tree, Node node, int parentState) {
        int nodeNum = node.getNr();
        Node parent = node.getParent();
        double[] conditionalProbabilities = new double[stateCount];

        if (!node.isLeaf()) {
            handleInternalNode(node, parent, parentState, conditionalProbabilities);
        } else {
            handleLeafNode(node);
        }
    }

    private void handleInternalNode(Node node, Node parent, int parentState, double[] conditionalProbabilities) {
        int nodeNum = node.getNr();

        if (parent == null && !useOrigin) {
            handleRootNode(node, conditionalProbabilities);
        } else {
            handleNonRootInternalNode(node, parent, parentState, conditionalProbabilities);
        }

        traverseChildren(node);
    }

    private void handleRootNode(Node node, double[] conditionalProbabilities) {
        beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNr()), Beagle.NONE, conditionalProbabilities);
        double[] rootFrequencies = getRootFrequencies();
        applyRootFrequencies(conditionalProbabilities, rootFrequencies);
        reconstructedStates[node.getNr()][0] = Randomizer.randomChoicePDF(conditionalProbabilities);
    }

    private void handleNonRootInternalNode(Node node, Node parent, int parentState, double[] conditionalProbabilities) {
        double[] partialLikelihood = new double[stateCount];
        beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNr()), Beagle.NONE, partialLikelihood);

        if (parent == null && useOrigin) {
            parentState = handleOriginNode(node);
        } else {
            beagle.getTransitionMatrix(matrixBufferHelper.getOffsetIndex(node.getNr()), probabilities);
        }

        calculateConditionalProbabilities(partialLikelihood, parentState, conditionalProbabilities);
        reconstructedStates[node.getNr()][0] = Randomizer.randomChoicePDF(conditionalProbabilities);
    }

    private double[] getRootFrequencies() {
        return rootFrequenciesInput.get() == null ? 
                substitutionModel.getFrequencies() : 
                rootFrequenciesInput.get().getFreqs();
    }

    private void applyRootFrequencies(double[] probabilities, double[] frequencies) {
        for (int i = 0; i < stateCount; i++) {
            probabilities[i] *= frequencies[i];
        }
    }

    private int handleOriginNode(Node node) {
        double[] oPs = new double[m_nStateCount];
        System.arraycopy(originPartials, 0, oPs, 0, originPartials.length);
        
        double[] rootFrequencies = getRootFrequencies();
        applyRootFrequencies(oPs, rootFrequencies);
        
        int parentState = Randomizer.randomChoicePDF(oPs);
        probabilities = rootTransitionMatrix;
        
        return parentState;
    }

    private void calculateConditionalProbabilities(double[] partialLikelihood, int parentState, double[] conditionalProbabilities) {
        int parentIndex = parentState * stateCount;
        for (int i = 0; i < stateCount; i++) {
            conditionalProbabilities[i] = partialLikelihood[i] * probabilities[parentIndex + i];
        }
    }

    private void traverseChildren(Node node) {
        Node child1 = node.getChild(0);
        traverseSample(treeInput.get(), child1, reconstructedStates[node.getNr()][0]);

        Node child2 = node.getChild(1);
        traverseSample(treeInput.get(), child2, reconstructedStates[node.getNr()][0]);
    }

    private void handleLeafNode(Node node) {
        System.arraycopy(tipStates[node.getNr()], 0, reconstructedStates[node.getNr()], 0, reconstructedStates[node.getNr()].length);
    }
}
