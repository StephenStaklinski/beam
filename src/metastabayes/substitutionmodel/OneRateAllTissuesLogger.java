package metastabayes.substitutionmodel;

import java.io.PrintStream;

import beast.base.evolution.datatype.UserDataType;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger;

public class OneRateAllTissuesLogger extends SVSGeneralSubstitutionModelLogger{

	// no need to log all tissue transition names since there is only one rate parameter for all
    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*"))
                ? "geoSubstModel"
                : getID();
        String relRatePrefix = mainID + ".relGeoRate_";

        UserDataType dataType = dataTypeInput.get();
        
        String iStr = "all_tissues";
        out.print(relRatePrefix + iStr);

    }
	
	// only log one parameter since all tissue transitions are fixed to one rate
    @Override
    public void log(long nSample, PrintStream out) {
        int count = 0;
        out.print(model.ratesInput.get().getArrayValue(count)
                *model.indicator.get().getArrayValue(count) + "\t");
    }
	
}
