package Extension.explicit;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import edu.jas.structure.Value;
import Extension.explicit.ArraysSort;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;
/**
 * Implementation of the bisimilarity algorithm from the paper "Simple
 * O(m.log(n)) Time Markov Chain Lumping" by Antti Valmari and Giuliana
 * Franceschinis.
 * 
 * Uses the second presented data structure (i.e. the more recent one): States
 * and blocks are represented by numbers. All states (that is, their numbers)
 * are in an array elems so that states that belong to the same block are next
 * to each other. The segment for a block is further divided to a first part
 * that contains the marked states and second part that contains the rest. There
 * is another array that, given the number of a state, returns its location in
 * elems. A third array denotes the block that each state belongs to. Three
 * arrays are indexed by block numbers; they tell where the segment for the
 * block in elems starts and ends, and where is the borderline between the
 * marked and other states.
 * 
 * @author Zainab Fatmi
 */
public class Valmari<Value> extends AbstractBisimulation<Value>{

	public Valmari(PrismComponent parent) throws PrismException {
		super(parent);
	}
	
	public static final double ACCURACY = 1E-5;

	private int[] elems; // a list of states, with states in the same block next to each other
	private int[] location; // a state's location in elems
	private int[] block; // the block that each state belongs to
	private ArrayList<Integer> start; // the index of the start of the block in elems (inclusive)
	private ArrayList<Integer> end; // the index of the end of the block in elems (exclusive)
	private ArrayList<Integer> borderline; // the index of the first unmarked state of the block

	
	

	/**
	 * Computes the probabilistic bisimlarity for the labelled Markov Chain.
	 * 
	 * @return an array, indexed by state ID, containing the block ID each state
	 *         belongs to, i.e. the resulting partition of the states of the
	 *         labelled Markov Chain
	 */
	public int[] evaluate(DTMCSimple<Value> dtmc, List<BitSet> propBSs){ 
		
		Evaluator<Value> eval = dtmc.getEvaluator();
		this.elems = new int[numStates];
		this.location = new int[numStates];
		this.block = partition.clone();
		this.start = new ArrayList<Integer>(numBlocks);
		this.end = new ArrayList<Integer>(numBlocks);
		this.borderline = new ArrayList<Integer>(numBlocks);
	
		
		int[] count = new int[numBlocks];
		for (int i = 0; i < numStates; i++) {
			count[this.block[i]]++;
		}
		this.end.add(count[0]);
		for (int i = 1; i < numBlocks; i++) {
			count[i] += count[i - 1];
			this.end.add(count[i]);
		}
		for (int i = numStates - 1; i >= 0; i--) {
			count[this.block[i]]--;
			this.elems[count[this.block[i]]] = i;
			this.location[i] = count[this.block[i]];
		}
		for (int i = 0; i < numBlocks; i++) {
			this.start.add(count[i]);
			this.borderline.add(count[i]);
			
		}
	
	
		HashSet<Integer> UB = new HashSet<Integer>(); // potential splitters
		for (int i = 0; i < this.start.size(); i++) {
			UB.add(i);
		}
		ArrayList<Integer> BT = new ArrayList<Integer>(10); // touched blocks
		ArrayList<Integer> ST = new ArrayList<Integer>(10); // touched states
		double[] w = new double[numStates]; // total probability of going to the splitter
		while (!UB.isEmpty()) {
			
			int splitter =  (int) UB.iterator().next();//UB.first();
			UB.remove(splitter);
			ST.clear();
			for (int i = (int) this.start.get(splitter); i < (int) this.end.get(splitter); i++) {
				for (int j = 0; j < numStates; j++) {
					double prob = eval.toDouble(dtmc.getProbability(j, this.elems[i]));
					if (prob > 0) { //if (this.W[j][this.elems[i]] > 0) { // predecessors
						if (w[j] == 0) {
							ST.add(j); // first transition to the splitter
							w[j] = prob;
						} else {
							w[j] += prob;
						}
					}
				}
			}
			for (int j = 0; j < ST.size(); j++) {
				int s = (int) ST.get(j);
				int b = this.block[s];
				if (this.start.get(b) == this.borderline.get(b)) {
					BT.add(b); // no marked states yet
				}
				this.mark(s, b);
			}
			for (int j = 0; j < BT.size(); j++) {
				int b = (int) BT.get(j);
				int b1;
				int l = 0; // number of new blocks
				if (this.borderline.get(b) == this.end.get(b)) {
					this.borderline.set(b, this.start.get(b));
					b1 = b; // b was empty so transfer identity
				} else {
					// create b1
					b1 = this.start.size();
					this.split(b);
					l++;
				}
				double pmc = this.pmc(w, b1);
				for (int i = (int) this.start.get(b1); i < (int) this.end.get(b1); i++) {
					if (!isEqual(w[elems[i]], pmc)) {
						mark(elems[i], b1);
					}
				}
				if (this.borderline.get(b1) != this.start.get(b1)) {
					// create b2
					int b2 = this.start.size();
					this.split(b1);
					l++;
					ArraysSort sorter = new ArraysSort(w);
					sorter.sort(this.elems, this.start.get(b2), this.end.get(b2));
					this.location[this.elems[this.start.get(b2)]] = this.start.get(b2);
					for (int i = this.start.get(b2) + 1; i < this.end.get(b2); i++) {
						this.location[this.elems[i]] = i; // update locations
						if (!isEqual(w[this.elems[i]], w[this.elems[i - 1]])) {
							this.borderline.set(this.block[this.elems[i]], i);
							this.split(this.block[this.elems[i]]);
							l++;
						}
					}
				}
				// add the new blocks as potential splitters
				int max = b; // the largest block
				for (int i = this.start.size() - l; i < this.start.size(); i++) {
					UB.add(i);
					if (((int) this.end.get(i) - (int) this.start.get(i)) > ((int) this.end.get(max) - (int) this.start.get(max))) {
						max = i;
					}
				}
				// remove the maximum
				if (max != b && !UB.contains(b)) {
					UB.add(b);
					UB.remove(max);
				}
			}
			
			Collections.copy(borderline, start); 
			
			// clear all
			BT.clear();
			for (int j = 0; j < ST.size(); j++) {
				w[(int)ST.get(j)] = 0;
			}
		}	
		return this.block;
	}
	

	/**
	 * Marks the given state.
	 * 
	 * @param s the state to mark
	 * @param b the block which the state belongs to
	 */
	private void mark(int s, int b) {
		int border = (int) this.borderline.get(b);
		int t = this.elems[border];
		this.elems[this.location[s]] = t;
		this.elems[border] = s;
		this.location[t] = this.location[s];
		this.location[s] = border;
		this.borderline.set(b, border + 1);
			
	}

	/**
	 * Splits a block, by creating a new block with the marked states.
	 * 
	 * @param b the block to split
	 */
	private void split(int b) {
		for (int i = (int) this.start.get(b); i < (int) this.borderline.get(b); i++) {
			this.block[this.elems[i]] = this.start.size(); // update the block of the marked elements
		}
		// create the new block
		this.start.add(this.start.get(b));
		this.borderline.add(this.start.get(b));
		this.end.add(this.borderline.get(b));
		this.start.set(b, this.borderline.get(b));
	}

	/**
	 * Estimates the possible majority of the w[s] for all states, s, in the block.
	 * 
	 * @param w     an array of size n, which stores the total probability of going
	 *              to the current splitter for each state
	 * @param block the block for which to find the possible majority candidate
	 * @return the possible majority candidate
	 */
	private double pmc(double[] w, int block) {
		int count = 0;
		double pmc = 0;
		for (int i = (int) this.start.get(block); i < (int) this.end.get(block); i++) {
			if (count == 0) {
				pmc = w[elems[i]];
				count++;
			} else if (isEqual(pmc, w[elems[i]])) {
				count++;
			} else {
				count--;
			}
		}
		return pmc;
	}

	/**
	 * Determines whether the given double values are within epsilon of each other
	 * (absolute error).
	 * 
	 * @param value1 the first double value.
	 * @param value2 the second double value.
	 * @return true if the values are equal, false otherwise.
	 */
	private boolean isEqual(double value1, double value2) {
		return (Math.abs(value1 - value2) < ACCURACY);
	}

	/**
	 * Converts the current partition from the six arrays data structure to a set of
	 * sets (containing state IDs). Used for testing.
	 * 
	 * @return the current partition of the states of the labelled Markov Chain
	 */
	public Set<Set<Integer>> getPartition() {
		Set<Set<Integer>> partition = new java.util.HashSet<Set<Integer>>();
		for (int i = 0; i < this.start.size(); i++) {
			Set<Integer> block = new java.util.HashSet<Integer>();
			for (int j = (int) this.start.get(i); j < (int) this.end.get(i); j++) {
				block.add(this.elems[j]);
			}
			partition.add(block);
		}
		return partition;
	}

	
	public int[] compressArray(int[] arr) {
		Set<Integer> uniqueSet = new TreeSet<>();
        for (int num : arr) {
            uniqueSet.add(num);
        }
        
        numBlocks = uniqueSet.size();
        
        Map<Integer, Integer> mapping = new HashMap<>();
        int index = 0;
        for (int num : uniqueSet) {
            mapping.put(num, index++);
        }

        int[] result = new int[arr.length];
        for (int i = 0; i < arr.length; i++) {
            result[i] = mapping.get(arr[i]);
        }

        return result;
    }
	
	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	@Override
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs)
	{
		if (!(dtmc instanceof DTMCSimple)) 
			throw new IllegalArgumentException("Expected an instance of DTMCSimple.");
		   
		initialisePartitionInfo(dtmc, propBSs); 
		partition = evaluate((DTMCSimple<Value>) dtmc, propBSs);	
		partition = compressArray(partition);
		numStates = dtmc.getNumStates();
		
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States ");
		
		
		int[] stateOf = new int[numStates];
		for(int s = 0; s < numStates; s++)
			stateOf[partition[s]] = s;
			
		
		DTMCSimple<Value> dtmcNew = new DTMCSimple<Value>(numBlocks);
		for(int b = 0; b < numBlocks; b++) {
			//mainLog.println(numStates + "    " + this.start.get(b));
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
		

	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc The DTMC
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 * @return a two dimensional boolean array that captures for each state pair whether
	 * the states are probabilistic bisimilar
	 */
	@Override
	public boolean[] bisimilar(DTMC<Value> dtmc, List<BitSet> propBSs){
		
		
		if (!(dtmc instanceof DTMCSimple)) 
			throw new IllegalArgumentException("Expected an instance of DTMCSimple.");
		   
		initialisePartitionInfo(dtmc, propBSs); 
		int[] bl = evaluate((DTMCSimple<Value>) dtmc, propBSs);	
		
		int NumberOfStates = dtmc.getNumStates();
		boolean[] bisimilar = new boolean[NumberOfStates * NumberOfStates];
		
		for(int i = 0; i < numStates; i++) {
			for(int j = 0; j < numStates; j++) {
				bisimilar[i*numStates + j] = (bl[i] == bl[j]);
			}
		}
		
		
		return bisimilar;
	}
	
	
	
}
