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

import beam.likelihood.BeamBeagleTreeLikelihood;

/**
 * @author Stephen Staklinski
 */
@Description("Ancestral State Tree Likelihood that extends the custom Beagle likelihood with origin input" +
        " for the branch above the root.")
public class BeamBeagleAncestralTissueLikelihood extends BeamBeagleTreeLikelihood implements TreeTraitProvider {

    public Input<String> tagInput = new Input<String>("tag","label used to report trait in the output tree posterior samples", Validate.REQUIRED);
	public Input<List<LeafTrait>> leafTraitsInput = new Input<List<LeafTrait>>("leaftrait", "list of leaf traits", new ArrayList<LeafTrait>());
    

    @Override
    public void initAndValidate() {

        // Verify that there is only one site in the input data
        if (dataInput.get().getPatternCount() != 1) {
            throw new RuntimeException("Pattern count must be 1, representing a single site/label in the tissue labeling alignment per cell.");
        }

        // Initialize the superclass
        super.initAndValidate();

        // Retrieve and store the tag input for the posterior tree sample outputs
        this.tag = tagInput.get();

        // Initialize data type, and state count
        Alignment data = dataInput.get();
        dataType = data.getDataType();
        stateCount = dataType.getStateCount();

        // Initialize tree model and node count
        TreeInterface treeModel = treeInput.get();
        int nodeCount = treeModel.getNodeCount();

        // Initialize reconstructed states arrays
        reconstructedStates = new int[nodeCount][1];
        storedReconstructedStates = new int[nodeCount][1];

        // Add tree trait for states
        treeTraits.addTrait(STATES_KEY, new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public int[] getTrait(TreeInterface tree, Node node) {
                return getStatesForNode(tree, node);
            }

            public String getTraitString(TreeInterface tree, Node node) {
                return formattedState(getStatesForNode(tree, node), dataType);
            }
        });

        // Validate only one category rate for sites
        if (m_siteModel.getCategoryCount() > 1) {
            throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");
        }

        // Set up the tip states
        tipStates = new int[treeModel.getLeafNodeCount()][1];
        for (Node node : treeModel.getExternalNodes()) {
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
            int[] states = tipStates[node.getNr()];
            int code = data.getPattern(taxonIndex, 0);
            states[0] = data.getDataType().getStatesForCode(code)[0];
        }
    }

    @Override
    public void store() {
        super.store();

        for (int i = 0; i < reconstructedStates.length; i++) {
            System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
        }

        storedAreStatesRedrawn = areStatesRedrawn;
    }

    @Override
    public void restore() {

        super.restore();

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	likelihoodKnown = false;
    	boolean isDirty = super.requiresRecalculation();
    	int hasDirt = Tree.IS_CLEAN;
		isDirty |= super.requiresRecalculation();
		this.hasDirt |= hasDirt;

		return isDirty;
    }
    

    public int[] getStatesForNode(TreeInterface tree, Node node) {
        if (tree != treeInput.get()) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
        	calculateLogP();
        }

        if (!areStatesRedrawn) {
            // redraw ancestral states for each node
            traverseSample(tree, tree.getRoot(), null);
            areStatesRedrawn = true;
        }

        return reconstructedStates[node.getNr()];
    }

    
    @Override
    public double calculateLogP() {
        areStatesRedrawn = false;
        double marginalLogLikelihood = super.calculateLogP();
        likelihoodKnown = true;
        return logP;
    }

    protected TreeTraitProvider.Helper treeTraits = new Helper();

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }


    private static String formattedState(int[] state, DataType dataType) {
        StringBuffer sb = new StringBuffer();
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
     * Traverse (pre-order) the tree sampling the internal node states, assuming all partial likelihoods are already calculated.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeInterface tree, Node node, int[] parentState) {

        int nodeNum = node.getNr();
        Node parent = node.getParent();
        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[1];

        if (!node.isLeaf()) {

            // only use root frequencies before sampling if the root is truly the start (no origin being used)
            if (parent == null && !useOrigin) {

                beagle.getPartials(getPartialBufferHelper().getOffsetIndex(node.getNr()), Beagle.NONE, conditionalProbabilities);

                double[] rootFrequencies = substitutionModel.getFrequencies();
                if (rootFrequenciesInput.get() != null) {
                    rootFrequencies = rootFrequenciesInput.get().getFreqs();
                }

                for (int i = 0; i < stateCount; i++) {
                    conditionalProbabilities[i] *= rootFrequencies[i];
                }

                state[0] = Randomizer.randomChoicePDF(conditionalProbabilities);
                reconstructedStates[nodeNum][0] = state[0];
            } else {
                // This is an internal node (not the root) or it is the root but there is an origin so the root has a transition probability matrix
                double[] partialLikelihood = new double[stateCount];
                beagle.getPartials(getPartialBufferHelper().getOffsetIndex(node.getNr()), Beagle.NONE, partialLikelihood);
                beagle.getTransitionMatrix(getMatrixBufferHelper().getOffsetIndex(nodeNum), probabilities);

                // If this is the root then we need to sample the parent state at the origin
                if (parent == null && useOrigin) {

                    // Get the origin partials
                    double[] oPs = new double[m_nStateCount * categoryCount];
                    System.arraycopy(originPartials, 0, oPs, 0, originPartials.length);

                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }

                    for (int i = 0; i < stateCount; i++) {
                        oPs[i] *= rootFrequencies[i];
                    }

                    parentState = new int[1];
                    parentState[0] = Randomizer.randomChoicePDF(oPs);
                }

                int parentIndex = parentState[0] * stateCount;
                int childIndex = 0;

                for (int i = 0; i < stateCount; i++) {
                    conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                }

                state[0] = Randomizer.randomChoicePDF(conditionalProbabilities);
                reconstructedStates[nodeNum][0] = state[0];
            }

            // Traverse down the two child nodes
            Node child1 = node.getChild(0);
            traverseSample(tree, child1, state);

            Node child2 = node.getChild(1);
            traverseSample(tree, child2, state);
        } else {
            // This is an external leaf
            System.arraycopy(tipStates[nodeNum], 0, reconstructedStates[nodeNum], 0, reconstructedStates[nodeNum].length);
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	// useful when logging on a fixed tree in an AncestralTreeLikelihood that is logged, but not part of the posterior
    	hasDirt = Tree.IS_FILTHY;
    	calculateLogP();
        out.print(getCurrentLogP() + "\t");
    }


    protected DataType dataType;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;

    private String tag;
    private boolean areStatesRedrawn = false;
    private boolean storedAreStatesRedrawn = false;

    boolean likelihoodKnown = false;

    int stateCount;
    int[][] tipStates;

    public static final String STATES_KEY = "states";
}
