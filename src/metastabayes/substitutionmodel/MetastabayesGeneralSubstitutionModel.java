package metastabayes.evolution.substitutionmodel;




import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import beastclassic.inference.BayesianStochasticSearchVariableSelection;
import beast.base.core.Log;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.DefaultEigenSystem;



/**
 * @author dkuh004
 *         Date: Sep 18, 2011
 *         Time: 6:03:49 PM
 *
 *  * ported from beast1 - author: Marc Suchard
 */
@Description("SVS General Substitution Model")
public class MetastabayesGeneralSubstitutionModel extends ComplexSubstitutionModel implements BayesianStochasticSearchVariableSelection {

    public Input<BooleanParameter> indicator = new Input<BooleanParameter>("rateIndicator",
            "rates to indicate the presence or absence of transition matrix entries", Validate.REQUIRED);

    public Input<Boolean> isSymmetricInput = new Input<Boolean>("symmetric",
    		"Indicates the rate matrix is symmetric. " +
            "If true (default) n(n-1)/2 rates and indicators need to be specified. " +
            "If false, n(n-1) rates and indicators need to be specified.", Boolean.TRUE);

    private BooleanParameter rateIndicator;
    private boolean isSymmetric = false;

    @Override
    public void initAndValidate(){

        frequencies = frequenciesInput.get();

        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
        if (isSymmetricInput.get() && ratesInput.get().getDimension() != nrOfStates * (nrOfStates-1)/2) {
            throw new IllegalArgumentException("Dimension of input 'rates' is " + ratesInput.get().getDimension() + " but a " +
                    "rate matrix of dimension " + nrOfStates + "x" + (nrOfStates -1) + "/2" + "=" + nrOfStates * (nrOfStates -1) / 2 + " was " +
                    "expected");
        }
        if (!isSymmetricInput.get() && ratesInput.get().getDimension() != nrOfStates * (nrOfStates-1)) {
            int dim = nrOfStates * (nrOfStates -1);
            Log.warning.println("Dimension of input 'rates' is " + ratesInput.get().getDimension() + ". " +
            		"Changing dimension to " + nrOfStates + "x" + (nrOfStates -1)  + "=" + dim);
            if (ratesInput.get() instanceof Parameter.Base) {
            	((Parameter.Base)ratesInput.get()).setDimension(dim);
            }
            Log.warning.println("Setting dimension of indicators to " + dim);
            indicator.get().setDimension(dim);
            isSymmetric = true;
        }

        if (!isSymmetricInput.get() && eigenSystemClass.equals(DefaultEigenSystem.class.getName())) {
        	Log.warning.println("WARNING: eigenSystemClass is DefautlEigneSystem, which may cause trouble with asymtric analysis. "
        			+ "You may want to consider eigensystem='beast.evolution.substitutionmodel.RobustEigenSystem' instead.");
        }
        //try {
			eigenSystem = createEigenSystem();
		//} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
		//		| InvocationTargetException e) {
		//	throw new IllegalArgumentException(e);
		//}
        
//        if (robust.get()){
//            eigenSystem = new RobustEigenSystem(m_nStates);
//        }   else {
//            eigenSystem = new DefaultEigenSystem(m_nStates);
//        }
        
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];

        rateIndicator = indicator.get();
        super.initAndValidate();
    }


    public Parameter<?> getIndicators() {
        return rateIndicator;
    }

    public boolean validState() {
        return !updateMatrix || Utils.connectedAndWellConditioned(probability,this);
    }


    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        updateMatrix = true;
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

//    @Override
//    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
//        double distance = (fStartTime - fEndTime) * fRate;
//
//        int i, j, k;
//        double temp;
//
//        // this must be synchronized to avoid being called simultaneously by
//        // two different likelihood threads - AJD
//        synchronized (this) {
//            if (updateMatrix) {
//                setupRelativeRates();
//                setupRateMatrix();
//                try{
//                    eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
//                }catch(Exception e){
//                    updateMatrix = false;
//                    Arrays.fill(matrix, 0);
//                    return;
//
//                }
//                updateMatrix = false;
//            }
//        }
//
//        // is the following really necessary?
//        // implemented a pool of iexp matrices to support multiple threads
//        // without creating a new matrix each call. - AJD
//        // a quick timing experiment shows no difference - RRB
//        double[] iexp = new double[nrOfStates * nrOfStates];
//        // Eigen vectors
//        double[] Evec = eigenDecomposition.getEigenVectors();
//        // inverse Eigen vectors
//        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
//        // Eigen values
//        double[] Eval = eigenDecomposition.getEigenValues();
//        for (i = 0; i < nrOfStates; i++) {
//            temp = Math.exp(distance * Eval[i]);
//            for (j = 0; j < nrOfStates; j++) {
//                iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
//            }
//        }
//
//        int u = 0;
//        for (i = 0; i < nrOfStates; i++) {
//            for (j = 0; j < nrOfStates; j++) {
//                temp = 0.0;
//                for (k = 0; k < nrOfStates; k++) {
//                    temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
//                }
//
//                matrix[u] = Math.abs(temp);
//                u++;
//            }
//        }
//        int h = 3;
//        h++;
//    } // getTransitionProbabilities


    private double[] probability = null;

    @Override
    public void setupRelativeRates() {

        Function rates = this.ratesInput.get();
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i) * (rateIndicator.getValue(i)?1.:0.);
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {
    	if (!isSymmetricInput.get()) {
    		super.setupRateMatrix();
    		return;
    	}
        double [] fFreqs = frequencies.getFreqs();
        int count = 0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            for (int j = i+1; j <  nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[count];
               	rateMatrix[j][i] = relativeRates[count];
                count++;
            }
//             for (int j = i+1; j < m_nStates; j++) {
//                 m_rateMatrix[i][j] = relativeRates[i * (m_nStates -1) + j-1];
//             }
        }
        // bring in frequencies
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] *= fFreqs[j];
                rateMatrix[j][i] *= fFreqs[i];
            }
        }
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }
        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * fFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
    } // setupRateMatrix


    @Override
    protected boolean requiresRecalculation() {
    	// if the rate is only dirty for a value that the indicators block out,
    	// no recalculation is required, so check this first.
    	Function v = ratesInput.get(); 
    	if (v instanceof Parameter<?>) {
    		Parameter.Base<?> p = (Parameter.Base<?>) v;
    		if (p.somethingIsDirty()) {
    			if (frequencies.isDirtyCalculation()) {
			    	return super.requiresRecalculation();
    			}
        		Parameter<Boolean> indicator2 = indicator.get(); 
    			for (int i = 0; i < p.getDimension(); i++) {
    				if (indicator2.getValue(i) && p.isDirty(i)) {
    			    	return super.requiresRecalculation();
    				}
    			}
    			// no calculation is affected
    			return false;
    		}
    	}

    	return super.requiresRecalculation();
    }
    
    @Override
    public boolean canReturnComplexDiagonalization() {
    	return isSymmetric;
    }
}
