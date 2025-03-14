//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package Extension.explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import parser.State;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismNotSupportedException;

/**
 * Class to perform bisimulation minimisation for explicit-state models.
 */
public class Bisimulation<Value> extends PrismComponent
{
	// Local storage of partition info
	protected int numStates;
	protected int[] partition;
	protected int numBlocks;
	protected MDPSimple<Value> mdp;

	/**
	 * Construct a new Bisimulation object.
	 */
	public Bisimulation(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Perform bisimulation minimisation on a model.
	 * @param model The model
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	public Model<Value> minimise(Model<Value> model, List<String> propNames, List<BitSet> propBSs) throws PrismException
	{
		switch (model.getModelType()) {
		case DTMC:
			return minimiseDTMC((DTMC<Value>) model, propNames, propBSs);
		case CTMC:
			return minimiseCTMC((CTMC<Value>) model, propNames, propBSs);
		default:
			throw new PrismNotSupportedException("Bisimulation minimisation not yet supported for " + model.getModelType() + "s");
		}
	}

	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs)
	{
		double totalTime = 0;
		long startTimeTotal = System.nanoTime();
		// Create initial partition based on propositions
		initialisePartitionInfo(dtmc, propBSs);
		//printPartition(dtmc);

		// Iterative splitting
		boolean changed = true;
		while (changed)
			changed = splitDTMC(dtmc);
		//printPartition(dtmc);

		// Build reduced model
		DTMCSimple<Value> dtmcNew = new DTMCSimple<>(numBlocks);
		for (int i = 0; i < numBlocks; i++) {
			for (Map.Entry<Integer, Value> e : mdp.getChoice(i, 0)) {
				dtmcNew.setProbability(i, e.getKey(), e.getValue());
			}
		}
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States " + "and " + dtmcNew.getNumTransitions());
		attachStatesAndLabels(dtmc, dtmcNew, propNames, propBSs);
		
		long endTimeTotal = System.nanoTime();
		totalTime += (endTimeTotal - startTimeTotal) / 1_000_000_000.0;
		System.out.println("Total time taken for the bisim : " + totalTime + " seconds");
		return dtmcNew;
	}

	
	public boolean[] bisimilar(DTMC<Value> dtmc, List<BitSet> propBSs) {
		
		initialisePartitionInfo(dtmc, propBSs);

		boolean changed = true;
		while (changed)
			changed = splitDTMC(dtmc);
		
		boolean[] result = new boolean[numStates * numStates];
		for (int s = 0; s < numStates; s++) {
			for (int t = 0; t < numStates; t++) {
				result[s*numStates + t] = (partition[s] == partition[t]);
			}
		}

		return result;
	}

	
	/**
	 * Perform bisimulation minimisation on a CTMC.
	 * @param ctmc The CTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	protected CTMC<Value> minimiseCTMC(CTMC<Value> ctmc, List<String> propNames, List<BitSet> propBSs)
	{
		// Create initial partition based on propositions
		initialisePartitionInfo(ctmc, propBSs);
		//printPartition(ctmc);

		// Iterative splitting
		boolean changed = true;
		while (changed)
			changed = splitDTMC(ctmc);
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States");
		//printPartition(ctmc);

		// Build reduced model
		CTMCSimple<Value> ctmcNew = new CTMCSimple<>(numBlocks);
		for (int i = 0; i < numBlocks; i++) {
			for (Map.Entry<Integer, Value> e : mdp.getChoice(i, 0)) {
				ctmcNew.setProbability(i, e.getKey(), e.getValue());
			}
		}
		attachStatesAndLabels(ctmc, ctmcNew, propNames, propBSs);

		return ctmcNew;
	}

	/**
	 * Construct the initial partition based on a set of proposition bitsets.
	 * Store info in {@code numStates}, {@code numBlocks} and {@code partition}.
	 */
	protected void initialisePartitionInfo(Model<Value> model, List<BitSet> propBSs)
	{
		BitSet bs1, bs0;
		numStates = model.getNumStates();
		partition = new int[numStates];

		// Compute all non-empty combinations of propositions
		List<BitSet> all = new ArrayList<BitSet>();
		bs1 = (BitSet) propBSs.get(0).clone();
		bs0 = (BitSet) bs1.clone();
		bs0.flip(0, numStates);
		all.add(bs1);
		all.add(bs0);
		int n = propBSs.size();
		for (int i = 1; i < n; i++) {
			BitSet bs = propBSs.get(i);
			int m = all.size();
			for (int j = 0; j < m; j++) {
				bs1 = all.get(j);
				bs0 = (BitSet) bs1.clone();
 				bs0.andNot(bs);
				bs1.and(bs);
				if (bs1.isEmpty()) {
					all.set(j, bs0);
				} else {
					if (!bs0.isEmpty())
						all.add(bs0);
				}
			}
		}

		// Construct initial partition
		all.removeIf(BitSet::isEmpty);
		numBlocks = all.size();
		for (int j = 0; j < numBlocks; j++) {
			BitSet bs = all.get(j);
			for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
				partition[i] = j;
			}
		}
		
		
//		System.out.println("partition:");
//		for(int i = 0; i < numStates; i++)
//			System.out.print(partition[i] + " ");
//		System.out.println(" ");
	}

	/**
	 * Perform a split of the current partition, if possible, updating {@code numBlocks} and {@code partition}.
	 * @return whether or not the partition was split 
	 */
	private boolean splitDTMC(DTMC<Value> dtmc)
	{
		int s, a, i, numBlocksNew, numChoicesOld;
		Distribution<Value> distrNew;
		int partitionNew[];

		partitionNew = new int[numStates];
		numBlocksNew = 0;
		// Compute the signature for each state (i.e. the distribution for outgoing
		// transitions, lifted to the current partition)
		// For convenience, we just store them as an MDP, with action label equal to the index of the block
		mdp = new MDPSimple<>(numBlocks);
		for (s = 0; s < numStates; s++) {
			// Build lifted distribution
			Iterator<Map.Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(s);
			distrNew = new Distribution<>(dtmc.getEvaluator());
			while (iter.hasNext()) {
				Map.Entry<Integer, Value> e = iter.next();
				distrNew.add(partition[e.getKey()], e.getValue());
			}
			//mainLog.println(s + " " + distrNew.toString());
			// Store in MDP, update new partition
			a = partition[s];
			numChoicesOld = mdp.getNumChoices(a);
			i = mdp.addChoice(a, distrNew);
			if (i == numChoicesOld)
				mdp.setAction(a, i, numBlocksNew++);
			partitionNew[s] = (Integer) mdp.getAction(a, i);
		}
		// Debug info
		//mainLog.println("New partition: " + java.util.Arrays.toString(partitionNew));
		//mainLog.println("Signatures MDP: " + mdp.infoString());
		//mainLog.println("Signatures MDP: " + mdp);
		//mainLog.println("numBlocksNew:" + numBlocksNew);
		//mainLog.println("numBlocks:" + numBlocks);
		//try { mdp.exportToDotFile("mdp.dot"); } catch (PrismException e) {}
		// Update info
		boolean changed = numBlocks != numBlocksNew;
		if(changed)
		{
			numBlocks = numBlocksNew;
			partition = partitionNew;
		}
		return changed;
	}

	/**
	 * Display the current partition, showing the states in each block.
	 */
	@SuppressWarnings("unused")
	private void printPartition(Model<Value> model)
	{
		for (int i = 0; i < numBlocks; i++) {
			mainLog.print(i + ":");
			for (int j = 0; j < numStates; j++)
				if (partition[j] == i)
					if (model.getStatesList() != null)
						mainLog.print(" " + model.getStatesList().get(j));
					else
						mainLog.print(" " + j);
			mainLog.println();
		}
	}

	/**
	 * Attach a list of states to the minimised model by adding a representative state
	 * from the original model.
	 * Also attach information about the propositions (used for bisimulation minimisation)
	 * to the minimised model, in the form of labels (stored as BitSets).
	 * @param model The original model
	 * @param modelNew The minimised model
	 * @param propNames The names of the propositions
	 * @param propBSs Satisfying states (of the minimised model) for the propositions
	 */
	protected void attachStatesAndLabels(Model<Value> model, ModelExplicit<Value> modelNew, List<String> propNames, List<BitSet> propBSs)
	{
		// Attach states
		if (model.getStatesList() != null) {
			List<State> statesList = model.getStatesList();
			List<State> statesListNew = new ArrayList<State>(numBlocks);
			for (int i = 0; i < numBlocks; i++) {
				statesListNew.add(null);
			}
			for (int i = 0; i < numStates; i++) {
				if (statesListNew.get(partition[i]) == null)
					statesListNew.set(partition[i], statesList.get(i));
			}
			modelNew.setStatesList(statesListNew);
		}

		// Build/attach new labels
		int numProps = propBSs.size();
		for (int i = 0; i < numProps; i++) {
			String propName = propNames.get(i);
			BitSet propBS = propBSs.get(i);
			BitSet propBSnew = new BitSet();
			for (int j = propBS.nextSetBit(0); j >= 0; j = propBS.nextSetBit(j + 1))
				propBSnew.set(partition[j]);
			modelNew.addLabel(propName, propBSnew);
		}
	}

	
}
