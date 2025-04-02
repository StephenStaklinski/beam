package beam.substitutionmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;

/**
 * TideTree substitution model that can be used with the modified BEAGLE tree likelihood
 * under the assumption that editing happens during the entire experiment.
 *
 * @author Stephen Staklinski
 */
@Description("TideTree substitution model that can be used with the modified BEAGLE tree likelihood" +
            "under the assumption that editing happens during the entire experiment.")
public class BeamMutationSubstitutionModel extends SubstitutionModel.Base {

    /** Input for edit rates during the editing window */
    public final Input<List<RealParameter>> editRatesInput = new Input<>("editRates",
            "Rates at which edits are introduced into the genomic barcode during the editing window",
            new ArrayList<>(), Input.Validate.REQUIRED);

    /** Input for the silencing rate throughout the experiment */
    public final Input<RealParameter> silencingRateInput = new Input<>("silencingRate",
            "Rate at which barcodes are silenced throughout the entire experiment",
            Input.Validate.REQUIRED);

    /** Frequencies of states */
    private double[] frequencies;

    /** Edit rate parameter */
    private RealParameter editRate_;

    /** Silencing rate parameter */
    private RealParameter silencingRate_;

    /** Array of edit rates */
    private Double[] editRates;

    /** State representing missing data */
    private int missingDataState;

    /** Stored silencing rate for dirty state checking */
    private double storedSilencingRate;

    @Override
    public void initAndValidate() {
        // One state for each edit type + unedited + lost
        nrOfStates = editRatesInput.get().get(0).getDimension() + 2;

        editRate_ = editRatesInput.get().get(0);
        editRates = new Double[editRate_.getDimension()];
        for (int i = 0; i < editRate_.getDimension(); i++) {
            editRates[i] = editRate_.getValue(i);
        }

        // Add edit rates to rate matrix
        Double editRateSum = 0.0;
        for (double editRate : editRates) {
            if (editRate < 0) {
                throw new RuntimeException("All edit rates must be positive!");
            }
            editRateSum += editRate;
        }
        if (Math.abs(editRateSum - 1.0) > 1e-5) {
            throw new RuntimeException("Sum of edit rates must be 1.0, but it is " + editRateSum + "!");
        }

        silencingRate_ = silencingRateInput.get();
        double silencingRate = silencingRate_.getValue();
        storedSilencingRate = silencingRate;

        if (silencingRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }

        // Missing data state is the last edit
        missingDataState = nrOfStates - 1;

        // Center root frequency on the unedited first state, regardless of input frequencies
        // as this is a property of the barcodes
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // Calculate key parameters
        double silencingRate = silencingRate_.getValue();
        double delta = (startTime - endTime) * rate;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);
        double c = expOfDeltaLoss * (1 - Math.exp(-delta));

        // Initialize matrix to zeros
        Arrays.fill(matrix, 0.0);

        // Set absorbing state (bottom right corner)
        matrix[nrOfStates * nrOfStates - 1] = 1.0;

        // Set diagonal elements and loss probabilities in one pass
        for (int i = 0; i < nrOfStates - 1; i++) {
            int rowOffset = i * nrOfStates;
            matrix[rowOffset + i] = (i == 0) ? Math.exp(-delta * (1 + silencingRate)) : expOfDeltaLoss;
            matrix[rowOffset + (nrOfStates - 1)] = 1 - expOfDeltaLoss;
        }

        // Set edit probabilities (first row)
        for (int j = 1; j < nrOfStates - 1; j++) {
            matrix[j] = editRates[j - 1] * c;
        }
    }

    @Override
    public void store() {
        storedSilencingRate = silencingRate_.getValue();
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        // Since edit rates are fixed hyperparameters, the only thing that can change is the silencing rate
        return silencingRate_.getValue() != storedSilencingRate;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        // Return true for BEAGLE compatibility
        return true;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    /**
     * Returns the state representing missing data.
     *
     * @return The missing data state
     */
    public int getMissingState() {
        return missingDataState;
    }
}

