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
 * @author Stephen Staklinski
 **/
@Description("TideTree substitution model that can be used with the modified BEAGLE tree likelihood" +
            "under the assumption that editing happens during the entire experiment.")
public class BeamMutationSubstitutionModel extends SubstitutionModel.Base {

    final public Input<List<RealParameter>> editRatesInput = new Input<>("editRates",
                "Rates at which edits are introduced into the genomic barcode during the editing window", new ArrayList<>(), Input.Validate.REQUIRED);

        final public Input<RealParameter> silencingRateInput = new Input<>("silencingRate",
                "Rate at which barcodes are silenced throughout the entire experiment", Input.Validate.REQUIRED);


     @Override
    public void initAndValidate() {

        // one state for each edit type + unedited + lost
        nrOfStates = editRatesInput.get().get(0).getDimension() + 2;

        editRate_ = editRatesInput.get().get(0);
        editRates = new Double[editRate_.getDimension()];
        for (int i = 0; i < editRate_.getDimension(); i++) {
            editRates[i] = editRate_.getValue(i);
        }

        // add edit rates to rate matrix
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

        // missing data state is the last edit
        missingDataState = nrOfStates - 1;

        // center root frequency on the unedited first state, irregardless of input frequencies as this is a property of the barcodes
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;
    }


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // we only get here if either the silencing rate or the time has changed, so we need to recalculate the matrix
        // edit rates must be fixed hyperparameters for this setup

        // Calculate key parameters
        double silencingRate = silencingRate_.getValue();
        double delta = (startTime - endTime) * rate;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);
        double c = expOfDeltaLoss * (1 - Math.exp(-delta));

        // Initialize matrix to zeros
        Arrays.fill(matrix, 0.0);

        // Set diagonal elements
        matrix[0] = Math.exp(-delta * (1 + silencingRate));  // top left corner
        for (int i = 1; i < nrOfStates - 1; i++) {
            matrix[i * nrOfStates + i] = expOfDeltaLoss;  // middle diagonal elements
        }
        matrix[nrOfStates * nrOfStates - 1] = 1.0;  // bottom right corner (absorbing state)

        // Set loss probabilities (final column)
        for (int i = 0; i < nrOfStates - 1; i++) {
            matrix[i * nrOfStates + (nrOfStates - 1)] = 1 - expOfDeltaLoss;
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
        // since edit rates are fixed hyperparemeters, the only thing that can change is the silencing rate
        if (silencingRate_.getValue() != storedSilencingRate) {
            return true;
        }
        return false;
    }

    // return true for BEAGLE compatibility
    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {return null;}

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    public int getMissingState() {
        return missingDataState;
    }


    double[] frequencies;
    RealParameter editRate_;
    RealParameter silencingRate_;
    Double[] editRates;
    int missingDataState;

    private double storedSilencingRate;
}

