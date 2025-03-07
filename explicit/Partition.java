package explicit;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * A partition.
 * 
 * @author Eric Ruppert
 * @author Franck van Breugel
 * @author Hiva Karami
 */
public class Partition extends ArrayList<ArrayList<Integer>> {

	/*
	 * mapping of block IDs to indices
	 */
	private ArrayList<Integer> blockToIndex;

	/*
	 * mapping of indices to block IDs
	 */
	private ArrayList<Integer> indexToBlock;

	/*
	 * @inv. this.blockToIndex.size() == this.indexToBlock.size() == this.size()
	 * @inv. for all 0 <= i < this.size() : this.blockToIndex.get(this.indexToBlock.get(i)) == i
	 * @inv. for all 0 <= b < this.size() : this.indexToBlock.get(this.blockToIndex.get(b)) == b
	 */

	/*
	 * for all 0 <= s < this.index.size() : state s belongs to this.get(index[s])
	 */
	private int[] index;

	/*
	 * @inv. for all 0 <= b < this.size() : for all 0 <= s < this.index.size() : s in this.get(this.blockToIndex.get(b)) iff this.indexToBlock.get(this.index[s]) == b
	 */

	/**
	 * Initializes this partition of the given number of states
	 * from the given partition.
	 * 
	 * @param numberOfStates the number of states
	 * @pre. numberOfStates &gt; 0
	 * @param partition a partition of the states
	 * @pre. partition forms a partition of 0..(numberOfStates - 1)
	 */
	public Partition(int numberOfStates, List<BitSet> partition) {
		this.blockToIndex = new ArrayList<Integer>();
		this.indexToBlock = new ArrayList<Integer>();
		this.blockToIndex.ensureCapacity(partition.size());
		this.indexToBlock.ensureCapacity(partition.size());

		this.index = new int[numberOfStates];

		for (int block = 0; block < partition.size(); block++) {
			ArrayList<Integer> newBlock = new ArrayList<Integer>();
			this.add(newBlock);

			this.blockToIndex.add(block);
			this.indexToBlock.add(block);

			BitSet blockSet = partition.get(block);
			newBlock.ensureCapacity(blockSet.size());
			for (int state = blockSet.nextSetBit(0); state >= 0; state = blockSet.nextSetBit(state + 1)) {
				newBlock.add(state);
				this.index[state] = block;				
			}
		}
	}

	/**
	 * Returns the collection of states that belong to the block with the given ID.
	 * 
	 * @param block a block ID
	 * @pre. 0 &le; block &lt; this.size()
	 * @return the collection of states that belong to the block with the given ID
	 */
	public ArrayList<Integer> getStates(int block) {
		return this.get(this.blockToIndex.get(block));
	}

	/**
	 * Returns the ID of the block to which the given state belongs.
	 * 
	 * @param state a state
	 * @return the ID of the block to which the given state belongs
	 */
	public int getBlock(int state) {
		return this.indexToBlock.get(this.index[state]);
	}




	/**
	 * add all states of newBlock as a new block to this partition.
	 * @param newBlock states of the block
	 * 
	 */
	public void createNewBlock(ArrayList<Integer> newBlock) {
		int last = this.size();
		this.add(newBlock);
		this.blockToIndex.add(last);
		this.indexToBlock.add(last);
		for (Integer state : newBlock) {
			this.index[state] = last; 
		}
	}


	/**
	 * Removes all states of the given block whose ID doesn't match the index of the block.
	 * @param block ID of a block
	 */
	public void refine(int block) {

		int id = blockToIndex.get(block);
		ArrayList<Integer> refinedBlock = new ArrayList<Integer>();
		for (int state : this.getStates(block)) { 
			if (index[state] == id) { 
				refinedBlock.add(state); 
			}
		}
		this.setStates(block, refinedBlock);

	}

	/**
	 * 
	 * @param block	The identifier for the block whose states are to be updated.
	 * @param newStates	The new list of states to associate with the specified block.
	 * This method retrieves the index corresponding to the specified block from
	 * the `blockToIndex` list. It then replaces the current list of states at
	 * the retrieved index with the provided `newStates` list. Both the index lookup
	 * and the state replacement are efficient operations, assuming the underlying
	 * data structure is properly maintained.
	 */
	public void setStates(int block, ArrayList<Integer> newStates) {
		int index = this.blockToIndex.get(block);
		this.set(index, newStates); 
	}



	/**
	 * Removes all states of the given new block from the block with the given ID
	 * and add all those states as a new block to this partition.
	 * 
	 * @param block ID of a block
	 * @pre. 0 &le; block &lt; this.size()
	 * @param newBlock a block
	 * @pre. newBlock is a strict subset of this.getStates(block)
	 */
	public void refine(int block, ArrayList<Integer> newBlock) {
		int last = this.size();
		this.add(newBlock);
		this.blockToIndex.add(last);
		this.indexToBlock.add(last);
		this.getStates(block).removeAll(newBlock);
		for (Integer state : newBlock) {
			this.index[state] = last; 
		}

	}

	/**
	 * Swaps the blocks with the given IDs.
	 * 
	 * @param first a block ID
	 * @pre. 0 &le; first &lt; this.size()
	 * @param second a block ID
	 * @pre. 0 &le; second &lt; this.size()
	 */
	public void swap(int first, int second) {
		Collections.swap(this.blockToIndex, first, second);
		Collections.swap(this.indexToBlock, this.blockToIndex.get(first), this.blockToIndex.get(second));
	}

	/**
	 * Returns a string representation of this partition.
	 * 
	 * @return a string representation of this partition
	 */
	public String toString() {
		StringBuffer representation = new StringBuffer();

		representation.append("Block number\tIndex\t\tStates\n");
		for (int block = 0; block < this.size(); block++) {
			representation.append(block + "\t\t" + this.blockToIndex.get(block) + "\t\t" + this.get(this.blockToIndex.get(block)) + "\n");
		}

		representation.append("State\tBlock number\n");
		for (int state = 0; state < this.index.length; state++) {
			representation.append(state + "\t" + this.indexToBlock.get(this.index[state]) + "\n");
		}

		return representation.toString();
	}

	/**
	 * Returns the number of blocks
	 * @return Number of blocks
	 */
	public int getNumofBlocks() {
		return blockToIndex.size();
	}
}
