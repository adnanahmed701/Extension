package Extension.explicit;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import edu.jas.structure.Value;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;
/**
 * Decides which states of a labelled Markov chain are probabilistic bisimilar.  The implementation
 * is based on the bisimilarity algorithm from the paper "Efficient computation of
 * equivalent and reduced representations for stochastic automata" by Peter Buchholz.
 */
public class Buchholz<Value> extends AbstractBisimulation<Value>{

	public Buchholz(PrismComponent parent) throws PrismException {
		super(parent);
	}
	
	public static final double ACCURACY = 1E-5;	
	public static final int PRECISION = 3;
	
	private static class EquivalenceClass {
		private boolean initialized;
		private double value;
		private int next;

		/**
		 * Initializes this equivalence class as uninitialized (no state belonging to this equivalence
		 * class has been found yet).
		 */
		public EquivalenceClass() {
			this.initialized = false;
			this.value = 0;
			this.next = 0;
		}
	}
	/**
	 * Decides probabilistic bisimilarity for the given labelled Markov chain.
	 * This method calculates equivalence classes of states in a discrete-time Markov chain (DTMC) where
	 * bisimilar states are grouped together in the same set
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @return A list of sets, where each set represents an equivalence class of bisimilar states.
	 */
	public List<Set<Integer>> decide(DTMCSimple<Value> dtmc, List<BitSet> propBSs) {
	
		Evaluator<Value> eval = dtmc.getEvaluator();
		numStates = dtmc.getNumStates();
		List<Integer> indices = new ArrayList<Integer>();
		for (int state = 0; state < numStates; state++) {	
			int label = partition[state];
			if (!indices.contains(label)) {
				indices.add(label);
			}
		}

		int numberOfEquivalenceClasses = indices.size(); // number of equivalence classes
		List<Set<Integer>> classes = new ArrayList<Set<Integer>>(); // equivalence classes
		TreeSet<Integer> splitters = new TreeSet<Integer>(); // potential splitters
		int[] clazzOf = new int[numStates]; // for each state ID, the index of its equivalence class
		for (int clazz = 0; clazz < numberOfEquivalenceClasses; clazz++) {
			classes.add(new HashSet<Integer>());
			splitters.add(clazz);
		}
		
		for (int state = 0; state < numStates; state++) {
			int label = partition[state];
			int index = indices.indexOf(label);
			clazzOf[state] = index;
			classes.get(index).add(state);
		}

		double[] values = new double[numStates];
		while (!splitters.isEmpty()) {
			List<EquivalenceClass> split = new ArrayList<EquivalenceClass>();
			for (int clazz = 0; clazz < numberOfEquivalenceClasses; clazz++) {
				split.add(new EquivalenceClass());
			}
			int splitter = splitters.first();
			splitters.remove(splitter);
			
			// computing values
			Arrays.fill(values, 0);
			for (int target : classes.get(splitter)) {
				for (int source = 0; source < numStates; source++) { //source -> target 
					values[source] += eval.toDouble(dtmc.getProbability(source, target));
				}
			}
			

			for (int state = 0; state < numStates; state++) {
				int clazz = clazzOf[state];
				if (!split.get(clazz).initialized) {
					classes.set(clazz, new HashSet<Integer>());
					classes.get(clazz).add(state);
					split.get(clazz).initialized = true;
					split.get(clazz).value = values[state];
				} else {					
					if (Math.abs(split.get(clazz).value - values[state]) >= ACCURACY && split.get(clazz).next == 0) {
						splitters.add(clazz);
					}
					while (Math.abs(split.get(clazz).value - values[state]) >= ACCURACY && split.get(clazz).next != 0) {
						clazz = split.get(clazz).next;
					}
					if (Math.abs(split.get(clazz).value - values[state]) < ACCURACY) {
						clazzOf[state] = clazz;
						classes.get(clazz).add(state);
					} else {
						splitters.add(numberOfEquivalenceClasses);
						clazzOf[state] = numberOfEquivalenceClasses;
						split.get(clazz).next = numberOfEquivalenceClasses;
						split.add(new EquivalenceClass());
						split.get(numberOfEquivalenceClasses).initialized = true;
						split.get(numberOfEquivalenceClasses).value = values[state];
						classes.add(new HashSet<Integer>());
						classes.get(numberOfEquivalenceClasses).add(state);
						numberOfEquivalenceClasses++;
					}
				}
			}
		}
		
		return classes;	
	}
	
	
	
	@Override
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs){
		
		if (!(dtmc instanceof DTMCSimple)) 
			throw new IllegalArgumentException("Expected an instance of DTMCSimple.");
		initialisePartitionInfo(dtmc, propBSs); 
		List<Set<Integer>> classes = decide((DTMCSimple<Value>) dtmc, propBSs);
				
		
		numStates = dtmc.getNumStates();
		numBlocks = classes.size();
		
		int[] stateOf = new int[numStates];
		int id = 0;
		for (Set<Integer> clazz : classes) {
			for (Integer s : clazz) {
				partition[s] = id;
				stateOf[id] = s;
			}
			id++;
		}
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States ");
		DTMCSimple<Value> dtmcNew = new DTMCSimple<Value>(numBlocks);
		for(int b = 0; b < numBlocks; b++) {
			int s = stateOf[b];
			Iterator<Map.Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(s);
			while (iter.hasNext()) {
				Map.Entry<Integer, Value> e = iter.next();
				dtmcNew.addToProbability(b, partition[e.getKey()], e.getValue());
			}
			
		}
		attachStatesAndLabels(dtmc, dtmcNew, propNames, propBSs);

		return dtmcNew;
		
	}
	
	
	
	@Override
	public boolean[] bisimilar(DTMC<Value> dtmc, List<BitSet> propBSs){
		
		if (!(dtmc instanceof DTMCSimple)) 
			throw new IllegalArgumentException("Expected an instance of DTMCSimple.");
		   
		initialisePartitionInfo(dtmc, propBSs); 
		
		List<Set<Integer>> classes = decide((DTMCSimple<Value>) dtmc, propBSs);
		
	
		boolean[] bisimilar = new boolean[numStates * numStates];
		for (Set<Integer> clazz : classes) {
			for (Integer s : clazz) {
				for (Integer t : clazz) {
					bisimilar[s * numStates + t] = true;
				}
			}
		}
		
		return bisimilar;
	}
	
}
