package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.jas.structure.Value;
import parser.State;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismNotSupportedException;


/**
 * Abstract class that restricts bisimulation minimisation to specific models, 
 * such as {@link DTMCSimple}, and prevents the use of unsupported models 
 * (e.g., Continuous-Time Markov Chains (CTMC) or other types of DTMCs). 
 * This class provides a framework where additional minimisation algorithms 
 * can be added by extending this class.
 * 
 * <p>The goal is to ensure that only the {@link DTMCSimple} model is eligible for 
 * minimisation, while other models will raise exceptions. Extensions of this class 
 * can implement the minimisation logic for supported models.</p>
 */
public abstract class AbstractBisimulation<Value> extends Bisimulation<Value> {

    public AbstractBisimulation(PrismComponent parent) throws PrismException {
        super(parent);
    }

    @Override
    public Model<Value> minimise(Model<Value> model, List<String> propNames, List<BitSet> propBSs) throws PrismException
	{
    	if(model instanceof DTMCSimple)
			return minimiseDTMC((DTMC<Value>) model, propNames, propBSs);
    	else
			throw new PrismNotSupportedException("Bisimulation minimisation not yet supported for " + model.getModelType() + "s");
	}
   
   

    @Override
    protected CTMC<Value> minimiseCTMC(CTMC<Value> ctmc, List<String> propNames, List<BitSet> propBSs) {
        throw new UnsupportedOperationException("Bisimulation minimisation not yet supported for CTMCs");
    }

    
   
	

   
}
