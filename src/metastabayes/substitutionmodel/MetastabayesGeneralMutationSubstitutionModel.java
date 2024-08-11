package metastabayes.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import tidetree.substitutionmodel.EditAndSilencingModel;
import beast.base.evolution.tree.Node;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Stream;


/**
 * @author Stephen Staklinski
 **/
@Description("TideTree model that can be used with the modified BEAGLE tree likelihood" +
            "under the assumption that editing happens during the entire experiment.")
public class MetastabayesGeneralMutationSubstitutionModel extends EditAndSilencingModel {

    RealParameter editHeightP;
    RealParameter editDurationP;

     @Override
    public void initAndValidate() {

        // run normal TideTree setup
        super.initAndValidate();

        // override the edit height and edit duration to be the same
        editHeightP = editHeightInput.get();
        editDurationP = editHeightP;
    }

    // return true for BEAGLE compatibility
    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }

}

