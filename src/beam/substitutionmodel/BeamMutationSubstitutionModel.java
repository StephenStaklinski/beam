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
        inputEditRateSum = 0.0;
        for (int i=0; i<editRate_.getDimension(); i++){
            double editRate = editRate_.getValues()[i];
            if (editRate < 0) {
                throw new RuntimeException("All edit rates must be positive!");
            }
            inputEditRateSum += editRate;
        }
        if (inputEditRateSum < 0) {
            throw new RuntimeException("Sum of edit rates must be positive!");
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

        // calculate transition probabilities for loss process
        for (int i=0; i<nrOfStates; i++){
            for (int j=0; j<nrOfStates; j++){
                if ( i==j ){
                    // probability of staying in state
                    matrix[i*nrOfStates + j] = expOfDeltaLoss;
                }else if(j == nrOfStates-1){
                    // probability of loss
                    matrix[i*nrOfStates + j] = 1 - expOfDeltaLoss;
                }else{
                    // all other transitions are 0 probability
                    matrix[i*nrOfStates + j] = 0;
                }
            }
        }
        // set final diagonal element to 1 to ensure loss is an absorbing state, so once loss occurs it cannot be undone
        matrix[nrOfStates * nrOfStates - 1] = 1;
        
        // rate can change so scale the sum each time
        Double editRateSum = inputEditRateSum;
        editRateSum *= rate;

        // fill first row
        matrix[0] = Math.exp(-delta * (silencingRate + editRateSum));
        for (int i = 0; i < nrOfStates - 2; i++) {
            matrix[i + 1] = (editRates[i] * expOfDeltaLoss - editRates[i] * Math.exp(-delta * (silencingRate + editRateSum))) / editRateSum;
        }

        matrix[nrOfStates - 1] = 1 - expOfDeltaLoss;
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
    Double inputEditRateSum;

    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;

    private double storedSilencingRate;
}

