package explicit;

/*
 * Copyright (C)  2020  Zainab Fatmi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import explicit.DTMCSimple;
import explicit.ModelSimple;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;


/**
 * Decides which states of a labelled Markov chain have probabilistic bisimilarity 
 * distance zero, that is, they are probabilistic bisimilar.  The implementation is 
 * based on the algorithm from the paper "Optimal State-Space Lumping in Markov Chains" 
 * by Salem Derisavi, Holger Hermanns, and William Sanders.
 * 
 * @author Zainab Fatmi
 */
public class DerisaviRedBlack<Value> extends AbstractBisimulation<Value>{
	
	public DerisaviRedBlack(PrismComponent parent) throws PrismException {
		super(parent);
	}
	
	public static final double ACCURACY = 1E-8;
	public static final int PRECISION = 3;
	
	/**
	 * A class to represent the nodes of a splay tree.  Each node of the tree stores 
	 * a block and its probability of transitioning to the current splitter.
	 */
	enum Color {
	    RED, BLACK
	}
	private static class Node {
		private Block block;
		private double probability;
		private Color color;
		private Node parent;
		private Node left; // left child
		private Node right; // right child
		

		/**
		 * Initializes this node with an empty block, the given probability and the given parent node.
		 * 
		 * @param probability the probability of the states in the block of this node transitioning to the 
		 * current splitter
		 * @param parent the parent node
		 */
		public Node(double probability, Node parent) {
			this.block = new Block();
			this.probability = probability;
			this.parent = parent;
			this.left = null;
			this.right = null;
			this.color = Color.RED;
		}
		
		public boolean isRed() {
	        return color == Color.RED;
	    }
	}
	
	/**
	 * A splay tree.  Each node of the tree stores a block and its probability of 
	 * transitioning to the current splitter.
	 */
	private static class RedBlackTree {
		private Node root;
		private static final Color RED = Color.RED;
	    private static final Color BLACK = Color.BLACK;
		/**
		 * Initializes this splay tree as empty.
		 */
		public RedBlackTree() {
			this.root = null;
		}

		/**
		 * Inserts the state in its appropriate position in this splay tree. If the probability
		 * exists in the splay tree, adds the state to the block associated to the probability,
		 * otherwise creates a new node in this splay tree.
		 * 
		 * @param probability the probability of the state transitioning to the current splitter
		 * @param state a state
		 */
		public void insert(double probability, State state) {
		Node currentNode = root;
		    Node parent = null;
		    while (currentNode != null) {
		        parent = currentNode;
		        // Compare probabilities considering the ACCURACY
		        if (Math.abs(probability - currentNode.probability) < ACCURACY) {
		            // Probability is similar, add state to its block
		            currentNode.block.elements.add(state);
		            state.block = currentNode.block;
		            return;
		        } else if (probability < currentNode.probability) {
		            currentNode = currentNode.left;
		        } else {
		            currentNode = currentNode.right;
		        }
		    }

		    // If no similar probability node found, create a new node and insert it
		    Node newNode = new Node(probability, parent);
		    newNode.block.elements.add(state);
		    state.block = newNode.block;

		    // Insert the new node into the tree
		    if (parent == null) {
		        root = newNode;
		    } else if (probability < parent.probability) {
		        parent.left = newNode;
		    } else {
		        parent.right = newNode;
		    }

		    fixInsert(newNode);
		}

		/**
		 * Moves the given node to the root of the splay tree.
		 * 
		 * @param node a node
		 */
		private void fixInsert(Node node) {
	        Node parent, grandParent;

	        while (node != root && node.parent.isRed()) {
	            parent = node.parent;
	            grandParent = parent.parent;

	            if (parent == grandParent.left) {
	                Node uncle = grandParent.right;
	                if (uncle != null && uncle.isRed()) {
	                    parent.color = BLACK;
	                    uncle.color = BLACK;
	                    grandParent.color = RED;
	                    node = grandParent;
	                } else {
	                    if (node == parent.right) {
	                        rotateLeft(parent);
	                        node = parent;
	                        parent = node.parent;
	                    }
	                    rotateRight(grandParent);
	                    parent.color = BLACK;
	                    grandParent.color = RED;
	                    node = parent;
	                }
	            } else {
	                Node uncle = grandParent.left;
	                if (uncle != null && uncle.isRed()) {
	                    parent.color = BLACK;
	                    uncle.color = BLACK;
	                    grandParent.color = RED;
	                    node = grandParent;
	                } else {
	                    if (node == parent.left) {
	                        rotateRight(parent);
	                        node = parent;
	                        parent = node.parent;
	                    }
	                    rotateLeft(grandParent);
	                    parent.color = BLACK;
	                    grandParent.color = RED;
	                    node = parent;
	                }
	            }
	        }
	        root.color = BLACK;
	    }

		/**
		 * Rotates left at the given node.
		 * 
		 * @param node a node
		 */
		private void rotateLeft(Node node) {
	        Node child = node.right;
	        node.right = child.left;
	        if (child.left != null) {
	            child.left.parent = node;
	        }
	        child.parent = node.parent;
	        if (node.parent == null) {
	            root = child;
	        } else if (node == node.parent.left) {
	            node.parent.left = child;
	        } else {
	            node.parent.right = child;
	        }
	        child.left = node;
	        node.parent = child;
	    }

	    private void rotateRight(Node node) {
	        Node child = node.left;
	        node.left = child.right;
	        if (child.right != null) {
	            child.right.parent = node;
	        }
	        child.parent = node.parent;
	        if (node.parent == null) {
	            root = child;
	        } else if (node == node.parent.right) {
	            node.parent.right = child;
	        } else {
	            node.parent.left = child;
	        }
	        child.right = node;
	        node.parent = child;
	    }
	}

	/**
	 * A class to represent the blocks of the partition.
	 */
	private static class Block {
		private static int numberOfBlocks = 0;
		
		private int id; // for easier hashCode and equals methods
		private LinkedList<State> elements;
		private RedBlackTree tree;

		/**
		 * Initializes this block as empty and adds it to the partition.
		 */
		public Block() {
			this.id = numberOfBlocks++;
			this.elements = new LinkedList<State>();
			this.tree = new RedBlackTree();
			classes.add(this);
		}

		@Override
		public int hashCode() {
			return this.id;
		}

		@Override
		public boolean equals(Object object) {
			if (this != null && this.getClass() == object.getClass()) {
				Block other = (Block) object;
				return this.id == other.id;
			} else {
				return false;
			}
		}
	}

	/**
	 * A class to represent the states of the labelled Markov chain.
	 */
	private static class State {
		private int id;
		private Block block; // needed by the splay tree
		private double sum;
		private LinkedHashMap<State, Double> predecessors; // no need for successors

		/**
		 * Initializes this state with the given index.
		 * 
		 * @param id the non-negative ID of the state
		 */
		public State(int id) {
			this.id = id;
			this.sum = 0;
			this.predecessors = new LinkedHashMap<State, Double>();
		}

		@Override
		public int hashCode() {
			return this.id;
		}

		@Override
		public boolean equals(Object object) {
			if (this != null && this.getClass() == object.getClass()) {
				State other = (State) object;
				return this.id == other.id;
			} else {
				return false;
			}
		}
	}

	/**
	 * Partition of the states into blocks.
	 */
	private static LinkedList<Block> classes;
	
	/**
	 * Decides probabilistic bisimilarity distance zero for the given labelled Markov chain.
	 * 
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 */
	public void decide(DTMCSimple<Value> dtmc, List<BitSet> propBSs) {
	
		Evaluator<Value> eval = dtmc.getEvaluator();

		
		
	
		// start with an empty partition
		classes = new LinkedList<Block>();
		
		// create an empty block for each label and add it to the partition
		for (int i = 0; i < numBlocks; i++) {
			new Block();
		}
		
		// add the states to the blocks corresponding to the label of the state 
		State[] idToState = new State[numStates]; // map id to State
		for (int id = 0; id < numStates; id++) {
			State state = new State(id);
			idToState[id] = state;
			Block block = classes.get(partition[id]); //Block block = partition.get(chain.getLabel(id));
			block.elements.add(state);
			state.block = block;
		}
		for (int source = 0; source < numStates; source++) {
			for (int target = 0; target < numStates; target++) {
				double prob = eval.toDouble(dtmc.getProbability(source, target));
				if (prob != 0.0) {	
					idToState[target].predecessors.put(idToState[source], prob);
				}
			}
		}
		
		LinkedList<Block> potentialSplitters = new LinkedList<Block>(classes); // potential splitters
		Set<State> predecessors = new HashSet<State>(); // states that have a transition to the current splitter
		LinkedList<Block> partitioned = new LinkedList<Block>(); // blocks which will be partitioned

		while (!potentialSplitters.isEmpty()) {
			Block splitter = potentialSplitters.pop();

			predecessors.clear();
			for (State state : splitter.elements) {
				for (State predecessor : state.predecessors.keySet()) {
					predecessor.sum = 0;
				}
			}
			for (State state : splitter.elements) {
				for (Map.Entry<State, Double> entry : state.predecessors.entrySet()) {
					State predecessor = entry.getKey();
					predecessor.sum += entry.getValue();
					predecessors.add(predecessor);
				}
			}
			
			partitioned.clear();
			for (State state : predecessors) {
				Block block = state.block;
				block.elements.remove(state);
				block.tree.insert(state.sum, state);
				if (!partitioned.contains(block)) {
					partitioned.add(block);
				}
			}
			
			for (Block block : partitioned) {
				// traverse the subblock tree, adding subblocks and keeping track of the maximum
				Block max = block;
				LinkedList<Node> queue = new LinkedList<Node>();
				queue.add(block.tree.root);
				while (!queue.isEmpty()) {
					Node node = queue.removeFirst();
					if (node.left != null) {
						queue.add(node.left);
					}
					if (node.right != null) {
						queue.add(node.right);
					}
					if (node.block.elements.size() > max.elements.size()) {
						max = node.block;
					}
					potentialSplitters.add(node.block);
				}
				
				if (!potentialSplitters.contains(block) && !(max == block)) {
					potentialSplitters.add(block);
					potentialSplitters.remove(max);
				}

				if (block.elements.isEmpty()) {
					classes.remove(block);
					potentialSplitters.remove(block);
				} else {
					block.tree.root = null; // reset the splay tree
				}
			}
		}
		
	
	}
	
	
	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	@Override
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs){
		
		if (!(dtmc instanceof DTMCSimple)) 
			throw new IllegalArgumentException("Expected an instance of DTMCSimple.");
		initialisePartitionInfo(dtmc, propBSs); 
		decide((DTMCSimple<Value>) dtmc, propBSs);
		
		// Remove Blocks with empty elements list
        Iterator<Block> iterator = classes.iterator();
        while (iterator.hasNext()) {
            Block block = iterator.next();
            if (block.elements.size() == 0) {
                iterator.remove();
            }
        }
		
		numStates = dtmc.getNumStates();
		numBlocks = classes.size();
		
		int[] stateOf = new int[numStates];
		int id = 0;
		for(Block block : classes) {
			for (State s : block.elements) {
				partition[s.id] = id;
				stateOf[id] = s.id;
			}
			id++;
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
		
		decide((DTMCSimple<Value>) dtmc, propBSs);
		
		boolean[] bisimilar = new boolean[numStates * numStates];
		for (Block block : classes) {
			//System.out.println(block.elements.size());
			for (State s : block.elements) {
				for (State t : block.elements) {
					assert t.id < numStates;
					bisimilar[s.id * numStates + t.id] = true;
				}
			}
		}
		return bisimilar;
	}
	
}
