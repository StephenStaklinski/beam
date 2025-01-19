package beam.substitutionmodel;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Stephen Staklinski
 **/

@Description("Logger for BeamGtiTissueSubstitutionModel.")
public class BeamTissueSubstitutionModelLogger extends BEASTObject implements Loggable{

    public Input<GeneralSubstitutionModel> modelInput = new Input<>(
            "model",
            "Beam general substitution model.",
            Input.Validate.REQUIRED);

    public Input<UserDataType> dataTypeInput = new Input<>(
            "dataType",
            "User data type for the location data.  Used to generate " +
                    "more readable logs.",
            Input.Validate.REQUIRED);

    public Input<Boolean> useLocationNamesInput = new Input<>(
            "useLocationNames",
            "When true, use full names of locations in log rather than " +
                    "rate matrix indices.",
            true);

    private int nrOfStates;

    protected GeneralSubstitutionModel model;

    public BeamTissueSubstitutionModelLogger() { }

    @Override
    public void initAndValidate() {
        model = modelInput.get();
        nrOfStates = model.getStateCount();
    }

    /**
     * If available, retrieve string representation of location, otherwise
     * simply return the index as a string.
     *
     * @param i index of location
     * @return string representation
     */
    protected String getLocationString(int i) {
        if (useLocationNamesInput.get())
            return dataTypeInput.get().getCode(i);
        else
            return String.valueOf(i);
    }

    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*"))
                ? "geoSubstModel"
                : getID();
        String relRatePrefix = mainID + ".relGeoRate_";

        UserDataType dataType = dataTypeInput.get();

        for (int i=0; i<nrOfStates; i++) {
            String iStr = getLocationString(i);

            for (int j=0 ; j<nrOfStates; j++) {
                if (j==i) {
                    continue;
                }

                String jStr = getLocationString(j);
                out.print(relRatePrefix + iStr + "_" + jStr + "\t");
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {

        // Logging normalized rates directly from the rate matrix, not the input rate parameters used to setup the rate matrix.
        model.setupRateMatrix();
        double[][] rateMatrix = model.getRateMatrix();

        for (int i=0; i<nrOfStates; i++) {
            for (int j=0; j<nrOfStates; j++) {
                if (j==i) {
                    continue;
                }

                out.print(rateMatrix[i][j] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
