package beam.substitutionmodel;

import java.util.ArrayList;
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
        editRates = editRate_.getValues();

        // add edit rates to rate matrix
        Double editRateSum = 0.0;
        for (double editRate : editRates) {
            if (editRate < 0) {
            throw new RuntimeException("All edit rates must be positive!");
            }
            editRateSum += editRate;
        }
        if (Math.abs(editRateSum - 1.0) > 1e-5) {
            throw new RuntimeException("Sum of edit rates must be 1.0!");
        }

        silencingRate_ = silencingRateInput.get();
        double silencingRate = silencingRate_.getValue();
        storedSilencingRate = silencingRate;

        if (silencingRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }

        // center root frequency on the unedited first state, irregardless of input frequencies as this is a property of the barcodes
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // we only get here if either the silencing rate or the time has changed, so we need to recalculate the matrix
        // edit rates must be fixed hyperparameters for this setup

        // silencing rate can change so get it each time
        double silencingRate = silencingRate_.getValue();

        // calculate the branch time and scale it by the joint site model rate and clock rate
        double delta = (startTime - endTime) * rate;

        // calculate probability of no loss
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);

        // compute a constant across all edit states
        double c = expOfDeltaLoss * (1 - Math.exp(-delta));

        // fill the transition probability matrix
        for (int i=0; i<nrOfStates; i++){
            for (int j=0; j<nrOfStates; j++){
                if ( i==j ){
                    // setup the diagonal elements
                    if (i == 0) {
                        // top left corner
                        matrix[i] = Math.exp(-delta * (1 + silencingRate)); // assumes the edit rates sum to 1, which is checked in initAndValidate
                    } else if (i == nrOfStates - 1) {
                        // bottom right corner
                        matrix[i * nrOfStates + j] = 1; // absorbing loss state
                    } else {
                        // all other diagonal elements
                        matrix[i * nrOfStates + j] = expOfDeltaLoss;
                    }
                }else if(j == nrOfStates - 1){
                    // final column
                    matrix[i * nrOfStates + j] = 1 - expOfDeltaLoss;    // probability of loss
                }else if (i == 0){
                    // first row (excluding top left corner and top right corner)
                    matrix[j] = editRates[j - 1] * c;   // edit probabilities
                }else{
                    // all other elements
                    matrix[i*nrOfStates + j] = 0;
                }
            }
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


    double[] frequencies;
    RealParameter editRate_;
    RealParameter silencingRate_;
    Double[] editRates;

    private double storedSilencingRate;
}

