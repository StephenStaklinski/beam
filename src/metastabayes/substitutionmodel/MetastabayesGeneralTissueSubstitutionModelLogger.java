package metastabayes.substitutionmodel;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.UserDataType;

import java.io.PrintStream;

/**
 * @author Stephen Staklinski
 **/

@Description("Metastabayes General Substitution Model Logger. This code is modified from the FixedTreeAnalysis package by Remco Bouckaert")
public class MetastabayesGeneralTissueSubstitutionModelLogger extends BEASTObject implements Loggable{

    public Input<MetastabayesGeneralTissueSubstitutionModel> modelInput = new Input<>(
            "model",
            "Metastabayes general substitution model.",
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

    protected MetastabayesGeneralTissueSubstitutionModel model;

    public MetastabayesGeneralTissueSubstitutionModelLogger() { }

    @Override
    public void initAndValidate() {
        model = modelInput.get();
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

        for (int i=0; i<model.getStateCount(); i++) {
            String iStr = getLocationString(i);

            for (int j=0 ; j<model.getStateCount(); j++) {
                if (j==i)
                    continue;

                String jStr = getLocationString(j);

                out.print(relRatePrefix + iStr + "_" + jStr + "\t");
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        int count = 0;
        int nrOfStates = (int) Math.sqrt(model.structure.get().getDimension());

        for (int i=0; i<nrOfStates; i++) {
            for (int j=0; j<nrOfStates; j++) {
                if (j==i) {
                    count += 1;
                    continue;
                }

                int ms = model.structure.get().getValue(count);
                out.print(model.ratesInput.get().getArrayValue(ms) + "\t");

                count += 1;
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
