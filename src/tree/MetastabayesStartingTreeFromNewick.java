package metastabayes.tree;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.TreeParser;

import java.util.List;
import java.util.stream.IntStream;

/**
 * @author Stephen Staklinski
 */
@Description("This class uses the normal TreeParser abilities of reading in a starting tree" +
            " while also taking in the required inputs for the TideTree model.")
public class MetastabayesStartingTreeFromNewick extends TreeParser {
    
    final public Input<RealParameter> rootHeightInput = new Input<>("rootHeight", "Time from beginning of the experiment until sequencing", Input.Validate.REQUIRED);
    final public Input<RealParameter> editDurationInput = new Input<>("editDuration", "Time duration from edit start to edit stop", Input.Validate.REQUIRED);
    final public Input<RealParameter> editHeightInput = new Input<>("editHeight", "Time from the onset of edit until sequencing", Input.Validate.REQUIRED);

    // set up useful parameters
    double editHeight;
    double editDuration;
    double rootHeight;

    @Override
    public void initAndValidate() {

        rootHeight = rootHeightInput.get().getValue();
        editDuration = editDurationInput.get().getValue();
        editHeight = editHeightInput.get().getValue();

        if (editHeight > rootHeight){
            throw new RuntimeException("editHeight has to be smaller or equal than rootHeight");
        }
        super.initAndValidate();
    }
}
