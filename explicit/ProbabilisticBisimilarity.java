package explicit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

import explicit.DTMCSimple;
import explicit.Distribution;
import explicit.ModelSimple;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;



/**
 * Decides probabilistic bisimilarity for labelled Markov chains.
 * 
 * @author Eric Ruppert
 * @author Franck van Breugel
 * @author Hiva Karami
 */
public class ProbabilisticBisimilarity<Value> extends AbstractBisimulation<Value>{

	public ProbabilisticBisimilarity(PrismComponent parent) throws PrismException {
		super(parent);
	}
	
	
	protected static Partition Partition;
	protected static int numberOfLabels;
	/**
	 * Decides probabilistic bisimilarity for the given labelled Markov chain.
	 * 
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * 
	 */
	public void decide(DTMC<Value> dtmc, List<BitSet> propBSs) {
		Evaluator<Value> eval = dtmc.getEvaluator();
		numStates = dtmc.getNumStates(); 
		numberOfLabels = propBSs.size();

		/*
		 * initial partition: states are in the same block iff they have the same labelling
		 */
		final List<BitSet> initial = new ArrayList<BitSet>();

		final BitSet one = (BitSet) propBSs.get(0).clone();
		final BitSet notOne = (BitSet) one.clone();
		notOne.flip(0, numStates);
		if (!one.isEmpty()) { // ADDED
			initial.add(one);
		}
		if (!notOne.isEmpty()) { // ADDED
			initial.add(notOne);
		}
		for (int l = 1; l < numberOfLabels; l++) {
			BitSet set = propBSs.get(l);
			int size = initial.size();
			for (int i = 0; i < size; i++) {
				BitSet intersection = (BitSet) initial.get(i).clone();
				BitSet difference = (BitSet) intersection.clone();
				difference.andNot(set);
				intersection.and(set);
				if (!intersection.isEmpty()) {
					initial.set(i, intersection);
					if (!difference.isEmpty()) {
						initial.add(difference);
					}
				}
			}
		}

		int numberOfBlocks = initial.size();

		// partition
		Partition = new Partition(numStates, initial);

		// first potential splitter
		int first = 0; 

		// last potential splitter
		int last = numberOfBlocks - 1;

		/*
		 * States that transition to potential splitters and their signatures
		 */
		final SignatureList toCheck = new SignatureList();

		/*
		 * for all 0 <= s < numberOfStates : hasBeenChecked[s] = state s has been checked in the current round 
		 */
		final boolean[] hasBeenChecked = new boolean[numStates];

		/*
		 * for all 0 <= s < numberOfStates : index[state] = index of state in toCheck
		 */
		final int[] index = new int[numStates];

		/*
		 * for all 0 <= target < numberOfStates : predecessors[target] contains (source, probability)
		 * if source transitions to target with probability. 
		 */
		final List<Edge>[] predecessors = new List[numStates]; // use LinkedHashMap?
		for (int target = 0; target < numStates; target++) {
			predecessors[target] = new ArrayList<Edge>();
		}		
		for (int source = 0; source < numStates; source++) {

			Iterator<Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(source);
			while (iter.hasNext()) {
				Map.Entry<Integer, Value> e = iter.next();
				predecessors[e.getKey()].add(new Edge(source, eval.toDouble(e.getValue())));
			}
		}

		while (first <= last) {

			Arrays.fill(hasBeenChecked, false);
			toCheck.clearisFirst();
			toCheck.clear();
			for (int block = first; block <= last; block++) { // loop through new blocks created in previous round
				for (int target : Partition.getStates(block)) { // loop through states of block permutation[block]
					for (Edge edge : predecessors[target]) { // loop through incoming transitions of state
						int source = edge.getSource();
						double weight = edge.getWeight();
						if (!hasBeenChecked[source]) { // first time we are visiting source state in this round; create a new signature for it
							hasBeenChecked[source] = true;
							index[source] = toCheck.size();
							toCheck.add(new Signature(source, Partition.getBlock(source), block, weight), false);
						} else { // already visited source in this round; so simply add edge.weight to weight of block b in e.source's signature
							toCheck.get(index[source]).add(block, weight);
						}
					}
				}
			}

			toCheck.quicksort(0, toCheck.size()-1, 0);

			// split the blocks
			first = last + 1;

			int numberToCheck = toCheck.size();
			int numberChecked = 0;
			while (numberChecked < numberToCheck) {
				int maxSize = 0;
				int maxBlock = 0;
				int sumOfSizes = 0;
				int oldBlock = toCheck.get(numberChecked).getOldBlock();
				while (numberChecked < numberToCheck && toCheck.get(numberChecked).getOldBlock() == oldBlock) { // out of bounds
					ArrayList<Integer> newBlock = new ArrayList<Integer>();
					do { // loop through nodes with same signature, adding them to new block
						int state = toCheck.get(numberChecked).getState(); 
						newBlock.add(state);
						numberChecked++;
					} while (numberChecked < numberToCheck && !toCheck.isFirst(numberChecked));

					if (newBlock.size() != Partition.getStates(oldBlock).size()) { 
						// create new block
						Partition.createNewBlock(newBlock);
						//partition.refine(oldBlock, newBlock);
						last++;
						int size = newBlock.size();
						sumOfSizes += size;
						if (size > maxSize) { // keep track of largest sub-block of split block
							maxBlock = last;
							maxSize = size;
						}
					}
				}
				if (maxSize > Partition.getStates(oldBlock).size()-sumOfSizes) {

					Partition.refine(oldBlock);
					Partition.swap(maxBlock, oldBlock);
				}
			}
		}


		return;
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

		//initialisePartitionInfo(dtmc, propBSs); 
		decide(dtmc, propBSs);
		for(int i = 0; i < Partition.getNumofBlocks(); i++)
		{
			Partition.refine(i);
		}
		final boolean[] bisimilar = new boolean[numStates*numStates];
		for (List<Integer> block : Partition) {
			for (Integer s : block) {
				for (Integer t : block) {
					bisimilar[s*numStates+t] = true;
				}
			}
		}

		return bisimilar;
	}

	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	@Override
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs){

		double totalTime = 0;
		long startTimeTotal = System.nanoTime();
		decide(dtmc, propBSs);

		numBlocks = 0;
		partition = new int[numStates];	
		int numStates = dtmc.getNumStates();
		int[] stateOf = new int[numStates];
		int[] index = new int[numStates];
		for(int i = 0; i < numStates; i++) {
			int bl = Partition.getBlock(i);
			if(index[bl] == 0) {
				numBlocks++;
				index[bl] = numBlocks;
			}
			partition[i] = index[bl]-1;
			stateOf[index[bl]-1] = i;
				
		}

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
		
		long endTimeTotal = System.nanoTime();
		totalTime += (endTimeTotal - startTimeTotal) / 1_000_000_000.0;
		System.out.println("!!!!!!!!!!!!Total time taken for the newalgo : " + totalTime + " seconds");
		return dtmcNew;
	}




}
