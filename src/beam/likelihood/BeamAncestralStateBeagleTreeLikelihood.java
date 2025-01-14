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
public class BeamAncestralStateBeagleTreeLikelihood extends BeamBeagleTreeLikelihood implements TreeTraitProvider {

    public Input<String> tagInput = new Input<String>("tag","label used to report trait in the output tree posterior samples", Validate.REQUIRED);
	public Input<List<LeafTrait>> leafTraitsInput = new Input<List<LeafTrait>>("leaftrait", "list of leaf traits", new ArrayList<LeafTrait>());
    

    @Override
    public void initAndValidate() {

        // Initialize the superclass
        super.initAndValidate();

        // Check for BEAGLE library
        beagle = getBeagle();
        if (beagle == null) {
            throw new IllegalArgumentException("BEAGLE library not found. Please install BEAGLE library.");
        }

        // Retrieve and store the tag input for the posterior tree sample outputs
        this.tag = tagInput.get();

        // Initialize pattern count, data type, and state count
        patternCount = dataInput.get().getPatternCount();
        dataType = dataInput.get().getDataType();
        stateCount = dataType.getStateCount();

        // Initialize tree model and node count
        TreeInterface treeModel = treeInput.get();
        int nodeCount = treeModel.getNodeCount();

        // Initialize reconstructed states arrays
        reconstructedStates = new int[nodeCount][patternCount];
        storedReconstructedStates = new int[nodeCount][patternCount];

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

        // Validate site model input
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        if (m_siteModel.getCategoryCount() > 1) {
            throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");
        }

        // Initialize substitution model
        substitutionModel = (SubstitutionModel.Base) m_siteModel.substModelInput.get();

        // Initialize alignment data and probabilities array
        Alignment data = dataInput.get();
        int dim = data.getMaxStateCount() + 1;
        probabilities = new double[dim * dim];

        // Set up the tip states
        tipStates = new int[treeModel.getLeafNodeCount()][patternCount];
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
            for (int i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                states[i] = (statesForCode.length == 1) ? statesForCode[0] : code; // Causes ambiguous states to be ignored.
            }
        }

        // Initialize trait dimension
        traitDimension = tipStates[0].length;

        // Initialize leaf traits and parameters
        List<LeafTrait> leafTraits = leafTraitsInput.get();
        int leafTraitSize = leafTraits.size();
        leafNr = new int[leafTraitSize];

        List<String> taxaNames = dataInput.get().getTaxaNames();
        for (int i = 0; i < leafTraitSize; i++) {
            LeafTrait leafTrait = leafTraits.get(i);

            // Identify node
            String taxon = leafTrait.taxonName.get();
            int k = taxaNames.indexOf(taxon);
            if (k == -1) {
                throw new IllegalArgumentException("Could not find taxon '" + taxon + "' in tree");
            }
            leafNr[i] = k;
        }

        // Store initial tip states
        storedTipStates = Arrays.stream(tipStates).map(int[]::clone).toArray(int[][]::new);
    }

    @Override
    public void store() {
        super.store();

        for (int i = 0; i < reconstructedStates.length; i++) {
            System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
        }

        storedAreStatesRedrawn = areStatesRedrawn;
        storedJointLogLikelihood = jointLogLikelihood;
        
    }

    @Override
    public void restore() {

        super.restore();

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        jointLogLikelihood = storedJointLogLikelihood;
        
        // deal with ambiguous tips
        if (leafNr != null) {
			for (int i = 0; i < leafNr.length; i++) {
				int k = leafNr[i];
				int[] tmp = tipStates[k];
				tipStates[k] = storedTipStates[k];
				storedTipStates[k] = tmp;
				// Does not handle ambiguities or missing taxa
				likelihoodCore.setNodeStates(k, tipStates[k]);
			}
        }
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
    

    public DataType getDataType() {
        return dataType;
    }

    public int[] getStatesForNode(TreeInterface tree, Node node) {
        if (tree != treeInput.get()) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
        	try {
        		 calculateLogP();
        	} catch (Exception e) {
				throw new RuntimeException(e.getMessage());
			}
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNr()];
    }


    public void redrawAncestralStates() {
        jointLogLikelihood = 0;
        TreeInterface tree = treeInput.get();
        traverseSample(tree, tree.getRoot(), null);
        areStatesRedrawn = true;
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

    public void getStates(int tipNum, int[] states)  {
        // Saved locally to reduce BEAGLE library access
        System.arraycopy(tipStates[tipNum], 0, states, 0, states.length);
    }

	public void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        beagle.getPartials(getPartialBufferHelper().getOffsetIndex(number), cumulativeBufferIndex, partials);
	}

	public void getTransitionMatrix(int matrixNum, double[] probabilities) {

		beagle.getTransitionMatrix(getMatrixBufferHelper().getOffsetIndex(matrixNum), probabilities);
	}
    
    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeInterface tree, Node node, int[] parentState) {

        int nodeNum = node.getNr();

        Node parent = node.getParent();

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!node.isLeaf()) {

            // only use root frequencies before sampling if the root is truly the start (no origin being used)
            if (parent == null && !useOrigin) {

                double[] rootPartials = m_fRootPartials;

                double[] rootFrequencies = substitutionModel.getFrequencies();
                if (rootFrequenciesInput.get() != null) {
                    rootFrequencies = rootFrequenciesInput.get().getFreqs();
                }

                // This is the root node
                for (int j = 0; j < patternCount; j++) {
                	if (beagle != null) {
                		getPartials(node.getNr(), conditionalProbabilities);
                	} else {
                		System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);
                	}

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] *= rootFrequencies[i];
                    }
                    try {
                        state[j] = Randomizer.randomChoicePDF(conditionalProbabilities);
                    } catch (Error e) {
                        System.err.println(e.toString());
                        System.err.println("Please report error to Marc");
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];

                    jointLogLikelihood += Math.log(rootFrequencies[state[j]]);
                }

            } else {

                // This is an internal node, but not the root ... or it is the root but there is an origin so the root has a transition probability matrix
                double[] partialLikelihood = new double[stateCount * patternCount];

                // assumes beagle is available as a requirement for BEAM
                getPartials(node.getNr(), partialLikelihood);
                getTransitionMatrix(nodeNum, probabilities);

                if (parent == null && useOrigin) {
                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }

                    parentState = new int[patternCount];

                    for (int i = 0; i < patternCount; i++) {
                        double[] originPartials = originPartialsGlobal;

                        for (int j = 0; j < stateCount; j++) {
                            originPartials[j] *= rootFrequencies[j];
                        }
                        parentState[i] = Randomizer.randomChoicePDF(originPartials); 

                        // add the likelihood contribution of the origin
                        jointLogLikelihood += Math.log(rootFrequencies[parentState[i]]);
                    }
                }

                for (int j = 0; j < patternCount; j++) {
                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                    }

                    state[j] = Randomizer.randomChoicePDF(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];
                    double contrib = probabilities[parentIndex + state[j]];
                    jointLogLikelihood += Math.log(contrib);
                }
            }

            // Traverse down the two child nodes
            Node child1 = node.getChild(0);
            traverseSample(tree, child1, state);

            Node child2 = node.getChild(1);
            traverseSample(tree, child2, state);
        } else {
            // This is an external leaf
        	getStates(nodeNum, reconstructedStates[nodeNum]);
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

    private double jointLogLikelihood;
    private double storedJointLogLikelihood;

    boolean likelihoodKnown = false;

    int[][] storedTipStates;
	int[] leafNr;
	int traitDimension;
    int patternCount;
    int stateCount;
    int[][] tipStates;

    public static final String STATES_KEY = "states";
}
