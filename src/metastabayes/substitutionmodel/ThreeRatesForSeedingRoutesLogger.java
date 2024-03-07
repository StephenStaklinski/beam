package metastabayes.substitutionmodel;

import java.io.PrintStream;

import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger;

public class ThreeRatesForSeedingRoutesLogger extends SVSGeneralSubstitutionModelLogger{

	// no need to log all tissue transition names since there is only one rate parameter for each seeding topology
    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*"))
                ? "geoSubstModel"
                : getID();
        String relRatePrefix = mainID + ".relGeoRate_";

        for (int i=0; i<model.getStateCount(); i++) {
            String iStr = getLocationString(i);

            for (int j=model.isSymmetricInput.get() ? i+1 : 0 ; j<model.getStateCount(); j++) {
                if (j==i)
                    continue;

                String jStr = getLocationString(j);

                out.print(relRatePrefix + iStr + "_" + jStr + "\t");
            }
        }
    }
	
	// only log three parameters since all tissue transitions are fixed based on topology
    @Override
    public void log(long nSample, PrintStream out) {
        for (int i=0; i<model.getStateCount(); i++) {
            for (int j=0; j<model.getStateCount(); j++) {
                if (j==i) {
                	continue;
                } else if (i == 0) {
                    out.print(model.ratesInput.get().getArrayValue(0) + "\t"); 
                } else if (j == 0) {
                    out.print(model.ratesInput.get().getArrayValue(2) + "\t");
                } else {
                    out.print(model.ratesInput.get().getArrayValue(1) + "\t");
                }


            }
        }
    }
	
}
