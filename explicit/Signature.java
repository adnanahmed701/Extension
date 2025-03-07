package explicit;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * The signature of a state.  The signature contains
 * <ul>
 * <li>the ID of the state,</li>
 * <li>the ID of a block (the block to which the state belonged previously),</li>
 * <li>the IDs of blocks (those blocks that have just been split to which the state can transition), and<li>
 * <li>probabilities (the probabilities of the state transitioning to those blocks).</li>
 * </ul>
 * 
 * @author Eric Ruppert
 * @author Franck van Breugel
 */
public class Signature implements Comparable<Signature> {

	private final int state;
	private final int oldBlock;

	/**
	 * @inv. this.blocks.size() &gt; 0
	 * @inv. this.blocks.size() == this.probabilities.size()
	 * @inv. &forall; 0 &le; i &lt; this.blocks.size() - 1 : this.blocks.get(i) &lt; this.blocks.get(i + 1)
	 * @inv. &forall; 0 &le; i &lt; this.probabilities.size() : this.probabilities.get(i) &gt; 0
	 */
	private final List<Integer> blocks;
	private final List<Double> probabilities;

	static final double EPSILON = 1E-14;
	/* If difference between two probabilities is less than EPSILON, they are treated as equal.
	 */

	/**
	 * Initializes this signature with
	 * <ul>
	 * <li>the ID of the state,</li>
	 * <li>the ID of the block to which the state belonged previously,</li>
	 * <li>the ID of a block, and</li>
	 * <li>the probability of the state transitioning to that block.</li>
	 * </ul>
	 * 
	 * @param state ID of a state
	 * @param oldBlock ID of a block
	 * @param block ID of a block
	 * @param probability probability
	 * @pre. probability &gt; 0
	 */
	public Signature(int state, int oldBlock, int block, double probability) {
		this.state = state;
		this.oldBlock = oldBlock;

		this.blocks = new ArrayList<Integer>(1);
		this.blocks.add(block);

		this.probabilities = new ArrayList<Double>(1);
		this.probabilities.add(probability);
	}

	/**
	 * Initializes this signature randomly.
	 */
	public Signature() {
		final int STATES = 4; // maximum number of states
		final int BLOCKS = 5; // maximum number of blocks 
		final int PROBABILITIES = 3; // maximum number of different probabilities

		final Random RANDOM = new Random();

		this.state = RANDOM.nextInt(STATES);
		this.oldBlock = RANDOM.nextInt(BLOCKS);

		int size = 1 + RANDOM.nextInt(BLOCKS);
		this.blocks = new ArrayList<Integer>(size);
		this.probabilities = new ArrayList<Double>(size);
		for (int s = 0; s < size; s++) {
			int block;
			do {
				block = RANDOM.nextInt(BLOCKS);
			} while (this.blocks.contains(block));
			this.blocks.add(block);
			double probability = (1 + RANDOM.nextInt(PROBABILITIES)) / (double) PROBABILITIES;
			this.probabilities.add(probability);
		}

		Collections.sort(this.blocks);
	}

	/**
	 * Returns the state of this signature.
	 * 
	 * @return the state of this signature
	 */
	public int getState() {
		return this.state;
	}

	/**
	 * Returns the old block, that is, the block to which the state belonged previously, of this signature.
	 * 
	 * @return the old block of this signature
	 */
	public int getOldBlock() {
		return this.oldBlock;
	}

	/**
	 * Returns the probability of transitioning from the state to the block with the given index.
	 *
	 * @param index an index
	 * @pre. 0 &le; index &lt; this.probabilities.size();
	 * @return the probability of transitioning from the state to the block with the given index
	 */
	public double getProbability(int index) {
		return this.probabilities.get(index);
	}

	/**
	 * Returns the block with the given index.
	 *
	 * @param index an index
	 * @pre. 0 &le; index &lt; this.blocks.size();
	 * @return the block with the given index
	 */
	public int getBlock(int index) {
		return this.blocks.get(index);
	}

	/**
	 * Adds the given block and its probability to this signature.
	 * 
	 * @param block a block
	 * @pre. block &ge; this.blocks.get(this.blocks.size() - 1)
	 * @param probability the probability of the block
	 * @pre. probability &gt; 0
	 */
	public void add(int block, double probability) {
		int size = this.blocks.size() - 1;
		if (this.blocks.get(size) == block) {
			this.probabilities.set(size, this.probabilities.get(size) + probability);
		} else {
			this.blocks.add(block);
			this.probabilities.add(probability);
		}   
	}

	/**
	 * Returns the size of this signature.
	 *
	 * @return the size of this signature.
	 */
	public int size() {
		return this.blocks.size();
	}

	// FvB: now that I have removed blocks[0], this needs to be adjusted
	// ER: routine may not be needed anymore
	public boolean posNull(int pos) {
		// This method says whether the pos'th position exists
		return ((pos & 1) == 1 ? 2 * this.probabilities.size() <= pos : 2 * this.blocks.size() <= pos);
	}

	// Following sequence of comparison functions compare signatures according to following sequence of fields
	// blocks[0], ..., blocks[size-1], probabilitie[0], ... probabilities[size-1]
	// Positions 0...2*size-1 refer to positions within this sequence.
	// We have functions to compare at a particular position (either to a constant or to another signature)
	// and to compare from some position onwards.

	public int compareToAt(int val, int d) {
		// precondition 0<= d < size
		return this.blocks.get(d) - val;
	}

	public int compareToAt(double val, int d) {
		// precondition size <= d < 2*size
		Double diff = probabilities.get(d) - val;
		if (Math.abs(diff) < EPSILON)
			return 0;
		else if (diff < 0)
			return -1;
		else
			return 1;
	}	

	public int compareToAt(Signature other, int d) {
		// compare this to other at position d
		// precondition:  d should be between 0 and 2*size-1

		if (d <  size())
			return compareToAt(other.blocks.get(d), d);
		else 
			return compareToAt(other.probabilities.get(d - size()), d);
	}

	//Remark:  would be better structure to rewrite following routine using preceding routine, but that might be inefficient--lots of unnecessary tests

	public int compareToFrom(Signature other, int d) {
		// compare this to other from position d onwards in lexicographic ordering of following fields
		// precondition:  d should be between 0 and 2*size-1

		for (int i = d; i < this.blocks.size(); i++) {
			if (this.blocks.get(i) != other.blocks.get(i)) {
				return this.blocks.get(i) - other.blocks.get(i);
			}
		}

		for (int i = Math.max(d - this.probabilities.size(), 0); i < this.probabilities.size(); i++) {
			Double diff = this.probabilities.get(i) - other.probabilities.get(i);
			if (Math.abs(diff) > EPSILON) {
				if (diff < 0) {
					return -1;
				} else { // this.probabilities.get(i) > other.probabilities.get(i))
					return 1;
				} 
			} 
		}

		return 0;
	}

	public int compareTo(Signature other) {
		// When comparing two signatures, we use lexicographic order on following fields
		// oldBlock, size, blocks[0], ..., blocks[size-1], probabilities[0], ... probabilities[size-1];

		if (this.oldBlock != other.oldBlock) {
			return this.oldBlock - other.oldBlock;
		}

		if (this.blocks.size() != other.blocks.size()) {
			return this.blocks.size() - other.blocks.size();
		}

		return compareToFrom(other, 0);
	}

	/**
	 * Returns a string representation of this signature.
	 * This string has the following form:
	 * <ul>
	 * <li>the ID of the state,</li>
	 * <li>the ID of the block to which the state belonged previously,</li>
	 * <li>(b, p) pairs, where p is the probability that the state transitions to the block with ID b.</li>
	 * </ul> 
	 *
	 * @return a string representation of this signature.
	 */
	public String toString() {
		StringBuffer representation = new StringBuffer();
		representation.append(this.state + "\t" + this.oldBlock + "\t");
		for (int i = 0; i < this.blocks.size(); i++) {
			representation.append("(" + this.blocks.get(i) + ", " + String.format("%.2f", this.probabilities.get(i)) + ")\t");
		}
		return representation.toString();
	}

	/**
	 * Returns a string that can be used to sort signatures.  It consists of
	 * oldBlock, size, blocks[0] ... blocks[size-1] probabilities[0] ... probabilities[size-1]
	 * May be useful for debugging sorting (if list of signatures is sorted, then these strings
	 * are printed, the resulting strings should be in sorted order if you treat each tab field 
	 * as a single unit)
	 */
	public String toSortString() {
		StringBuffer representation = new StringBuffer();
		representation.append(oldBlock + "\t" + blocks.size());
		for (int i = 0; i < blocks.size(); i++)
			representation.append("\t" + blocks.get(i));
		for (int i = 0; i < probabilities.size(); i++)
			representation.append("\t" + probabilities.get(i));
		return representation.toString();
	}

}
