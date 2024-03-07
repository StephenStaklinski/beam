package metastabayes.substitutionmodel;

import java.io.PrintStream;

import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger;

public class OneRateAllTissuesLogger extends SVSGeneralSubstitutionModelLogger{

	// no need to log all tissue transition names since there is only one rate parameter for all
    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*"))
                ? "geoSubstModel"
                : getID();
        String relRatePrefix = mainID + ".relGeoRate_";

        String iStr = "all_tissues";
        out.print(relRatePrefix + iStr);

    }
	
	// only log one parameter since all tissue transitions are fixed to one rate
    @Override
    public void log(long nSample, PrintStream out) {
        out.print(model.ratesInput.get().getArrayValue(0) + "\t");
    }
	
}
