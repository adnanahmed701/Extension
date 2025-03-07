package explicit;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;


/**
 * A list of signatures.
 * 
 * @author Eric Ruppert
 * @author Franck van Breugel
 * @author Hiva Karami
 */
public class SignatureList extends ArrayList<Signature> {

	/*
	 * @inv. this.size() == isFirst.size()
	 */

	/**
	 * Whether the signature at the index of this list is the first of a block.
	 */
	private ArrayList<Boolean> isFirst;

	/**
	 * Initializes this list as empty.
	 */
	public SignatureList() {
		super();
		this.isFirst = new ArrayList<Boolean>();
	}

	/**
	 * Initializes this list as nonempty and random.
	 * 
	 * @param size maximal size of the list
	 */
	public SignatureList(int size) {
		final Random RANDOM = new Random();

		size = 1 + RANDOM.nextInt(size);
		for (int s = 0; s < size; s++) {
			this.add(new Signature());
		}

		this.isFirst = new ArrayList<Boolean>(size);
		for (int i = 0; i < size; i++) {
			this.isFirst.add(false);
		}
	}

	/**
	 * Initializes this list as a copy of the given list.
	 * 
	 * @param copy a list.
	 */
	public SignatureList(SignatureList copy) {
		super(copy);
		this.isFirst = new ArrayList<Boolean>(copy.isFirst);
	}

	/**
	 * Adds the given signature at the end of this list.
	 * 
	 * @param signature a signature
	 * @param isFirst whether the signature is first
	 */
	public void add(Signature signature, boolean isFirst) { // ER:  probably don't need second param which should always be false
		this.add(signature);
		this.isFirst.add(isFirst);
	}

	/**
	 * Swaps the elements at the given indices of this list.
	 * 
	 * @param i an index
	 * @pre. 0 &le; i &lt; this.size()
	 * @param j an index
	 * @pre. 0 &le; j &lt; this.size()
	 */
	private void swap(int i, int j) {
		Signature temp = this.get(i);
		this.set(i, this.get(j));
		this.set(j, temp);
	}

	/**
	 * Constant that determines when to switch from radix quick sort to standard sort.
	 * This constant should be greater or equal to zero.  (Zero means never switch.)
	 */
	private static final int THRESHOLD = 0;

	/**
	 * Sorts the sublist that starts at index low and ends at index high of this list.
	 * This sublist is denoted as [low..high].
	 * 
	 * 
	 * 
	 * @param low the index at which the sublist starts (inclusive)
	 * @pre. 0 &le; low &lt; this.size()
	 * @param high the index at which the sublist ends (inclusive)
	 * @pre. 0 &le; high &lt; this.size()
	 * @param d
	 * @pre. all elements of the sublist [low..high] have the same length, which is at least d
	 * @pre. all elements of the sublist [low..high] are identical up to position d - 1
	 * @pre. if low &gt; 0 then this.get(low) &gt; this.get(low-1) and this.isFirst(low)
	 */
	 void threeWayQuickSort(int low, int high, int d) {
		// sorts subarray A[lo..hi], assuming all elements of A[lo..hi] are identical up to position d-1 and all elements
		// of A[lo..hi] are of length at least d+1
		// Another precondition:  if lo>0 then A[lo] > A[lo-1] (and A[lo].first has already been set).
		// Also, identifies signatures that are the first of their group (i.e. the ones that differ from previous elt of A)
		// (The positions of a signature are block[0] weight[1] block[1] weight[2] block[2]...)
		if (high <= low + THRESHOLD) { // Use a standard sort routine for small sublists
			if (low < high) { 
				Collections.sort(this.subList(low, high));
				for (int i = low + 1; i <= high; i++) 
					this.isFirst.set(i, this.get(i - 1).compareToFrom(this.get(i), d) < 0);
			}
		} else {
			int size = get(low).size(); // size of block lists of all the signatures to be sorted
			//System.out.println("low " + low + " hight: " + high + " d: " + d + " size: " + size);
			int less = low;
			int greater = high;
			int i = low + 1;
			if (d < size) { // use blocks[d] to partition
				int pivot = this.get(low).getBlock(d);  // later, might want to choose random pivot instead 
				while (i <= greater) {
					// invariant lo <= less <= i-1 <= greater <= hi
					// invariant:  in position pos, A[lo..less-1] are < pivot and A[less..i] are = pivot and A[greater+1..hi] are > pivot
					int comparison = this.get(i).compareToAt(pivot, d);
					if (comparison < 0) {
						this.swap(less++, i++);
					} else if (comparison > 0) {
						this.swap(i, greater--);
					} else {
						i++;
					}
				}
			} else { // use probabilities[d-size] to partition
				int pos = d - size;
				double pivot = this.get(low).getProbability(pos); // later, might want to choose random pivot instead
				while (i <= greater) {
					// invariant lo <= less < i-1 <= greater <= hi
					// invariant:  in position pos, A[lo..less-1] are < pivot and A[less..i] are = pivot and A[greater+1..hi] are > pivot
					int comparison = this.get(i).compareToAt(pivot, pos);
					if (comparison < 0) {
						this.swap(less++, i++);
					} else if (comparison > 0) {
						this.swap(i, greater--);
					} else {
						i++;
					}
				}
			}
			// In position pos,
			// entries lo..less-1 are less than v,
			// entries less..greater are equal to v,
			// entries greater+1..hi are greater than v

			this.threeWayQuickSort(low, less - 1, d);    // sort first third recursively
			if (d+1 < 2*size)
				this.threeWayQuickSort(less, greater, d+1); // sort middle third recursively
			this.threeWayQuickSort(greater + 1, high, d); // sort last third recursively

			this.isFirst.set(less, true);
			if (greater < high) {
				this.isFirst.set(greater + 1, true);
			}

			// Postcondition:  A[lo..hi] now in sorted order and for lo<=i<=hi,  isFirst(i) is true
			// iff either A[i] is different from A[i-1] or i=0.
		}
	}

	 /**
	  * sort based on oldblock and size
	  * @param low the index at which the sublist starts (inclusive)
	  * @param high the index at which the sublist ends (inclusive)
	  * @param d is 0 or 1. 
	  */
	 void quicksort(int low, int high, int d) {
		 if(high < low)
			 return;
		 this.isFirst.set(low, true);
		 if(high == low)
			 return;
		
		int less = low;
		int greater = high;
		int i = low + 1;
		if(d == 0) {	// oldblock
			int pivot = this.get(low).getOldBlock();
			while (i <= greater) {
				int comparison = this.get(i).getOldBlock() - pivot;
				if (comparison < 0) {
					this.swap(less++, i++);
				} else if (comparison > 0) {
					this.swap(i, greater--);
				} else {
					i++;
				}
			}
			this.quicksort(low, less - 1, d);   
			this.quicksort(less, greater, d+1);
			this.quicksort(greater + 1, high, d); 
			this.isFirst.set(less, true);
			if (greater < high) {
				this.isFirst.set(greater + 1, true);
			}
		}
		else {			// size
			int pivot = this.get(low).size();
			//System.out.println(low + " " + high + " " + d + " " + pivot);
			while (i <= greater) {
				int comparison = this.get(i).size() - pivot;
				if (comparison < 0) {
					this.swap(less++, i++);
				} else if (comparison > 0) {
					this.swap(i, greater--);
				} else {
					i++;
				}
			}
			this.quicksort(low, less - 1, d);    
			this.threeWayQuickSort(less, greater, 0); 
			this.quicksort(greater + 1, high, d); 
			this.isFirst.set(less, true);
			if (greater < high) {
				this.isFirst.set(greater + 1, true);
			}
		}
		
		
		 
	 }
	
	 
	 /**
	 * Clears the isFirst ArrayList.
	 */
	public void clearisFirst() {
		this.isFirst.clear();
	}
	 
	/**
	 * Sorts this list.
	 */
	public void sort() {
		Collections.sort(this);
		this.isFirst.set(0, true);
		for (int i = 1; i < this.size() ; i++) {
			this.isFirst.set(i, this.get(i - 1).compareTo(this.get(i)) != 0);
		}
	}

	/**
	 * Tests whether the signature at the given index of this list is the first of a block.
	 * 
	 * @param index an index
	 * @pre. 0 &le; index &lt; this.size()
	 * @return true if the signature at the given index of this list is the first of a block,
	 * false otherwise
	 */
	public boolean isFirst(int index) {
		return this.isFirst.get(index);
	}

	/**
	 * Returns a string representation of this list of signatures.
	 * Each signature occurs on a separate line.
	 * 
	 * @return a string representation of this list of signatures
	 */
	public String toString() {
		StringBuffer representation = new StringBuffer("\n");
		for (int i = 0; i < this.size() ; i++) {
			if (this.isFirst.get(i)) {
				representation.append("* ");
			} else {
				representation.append("  ");
			}
			representation.append(this.get(i).toString() + "\n");
		}
		return representation.toString();
	}
}
