package explicit;

import java.util.Arrays;

/**
 * Re-implementation of Java's Arrays.sort for primitive integers with a
 * comparator.
 * 
 * @author Zainab Fatmi
 */
public class ArraysSort {
	private final double[] valueArray;

	/**
	 * Initializes the comparator with the indexed array.
	 * 
	 * @param array the array which is indexed by the array to be sorted.
	 */
	public ArraysSort(double[] array) {
		this.valueArray = array;
	}

	public int compare(int index1, int index2) {
		if (Math.abs(this.valueArray[index1] - this.valueArray[index2]) > 1e-10) {
			return this.valueArray[index1] > this.valueArray[index2] ? +1 : -1;
		} else {
			return 0;
		}
	}

	// From Arrays

	/**
	 * Sorts the specified range of the array into ascending order. The range to be
	 * sorted extends from the index {@code fromIndex}, inclusive, to the index
	 * {@code toIndex}, exclusive. If {@code fromIndex == toIndex}, the range to be
	 * sorted is empty.
	 *
	 * @implNote The sorting algorithm is a Dual-Pivot Quicksort by Vladimir
	 *           Yaroslavskiy, Jon Bentley, and Joshua Bloch. This algorithm offers
	 *           O(n log(n)) performance on all data sets, and is typically faster
	 *           than traditional (one-pivot) Quicksort implementations.
	 *
	 * @param a         the array to be sorted
	 * @param fromIndex the index of the first element, inclusive, to be sorted
	 * @param toIndex   the index of the last element, exclusive, to be sorted
	 *
	 * @throws IllegalArgumentException       if {@code fromIndex > toIndex}
	 * @throws ArrayIndexOutOfBoundsException if {@code fromIndex < 0} or
	 *                                        {@code toIndex > a.length}
	 */
	public void sort(int[] a, int fromIndex, int toIndex) {
		rangeCheck(a.length, fromIndex, toIndex);
		sort(a, 0, fromIndex, toIndex);
	}

	/**
	 * Checks that {@code fromIndex} and {@code toIndex} are in the range and throws
	 * an exception if they aren't.
	 */
	static void rangeCheck(int arrayLength, int fromIndex, int toIndex) {
		if (fromIndex > toIndex) {
			throw new IllegalArgumentException("fromIndex(" + fromIndex + ") > toIndex(" + toIndex + ")");
		}
		if (fromIndex < 0) {
			throw new ArrayIndexOutOfBoundsException(fromIndex);
		}
		if (toIndex > arrayLength) {
			throw new ArrayIndexOutOfBoundsException(toIndex);
		}
	}

	// from DualPivotQuicksort

	/**
	 * Max array size to use mixed insertion sort.
	 */
	private static final int MAX_MIXED_INSERTION_SORT_SIZE = 65;

	/**
	 * Max array size to use insertion sort.
	 */
	private static final int MAX_INSERTION_SORT_SIZE = 44;

	/**
	 * Min array size to try merging of runs.
	 */
	private static final int MIN_TRY_MERGE_SIZE = 4 << 10;

	/**
	 * Min size of the first run to continue with scanning.
	 */
	private static final int MIN_FIRST_RUN_SIZE = 16;

	/**
	 * Min factor for the first runs to continue scanning.
	 */
	private static final int MIN_FIRST_RUNS_FACTOR = 7;

	/**
	 * Max capacity of the index array for tracking runs.
	 */
	private static final int MAX_RUN_CAPACITY = 5 << 10;

	/**
	 * Threshold of mixed insertion sort is incremented by this value.
	 */
	private static final int DELTA = 3 << 1;

	/**
	 * Max recursive partitioning depth before using heap sort.
	 */
	private static final int MAX_RECURSION_DEPTH = 64 * DELTA;

	/**
	 * Sorts the specified array using the Dual-Pivot Quicksort and/or other sorts
	 * in special-cases, possibly with parallel partitions.
	 *
	 * @param a    the array to be sorted
	 * @param bits the combination of recursion depth and bit flag, where the right
	 *             bit "0" indicates that array is the leftmost part
	 * @param low  the index of the first element, inclusive, to be sorted
	 * @param high the index of the last element, exclusive, to be sorted
	 */
	@SuppressWarnings("unused")
	void sort(int[] a, int bits, int low, int high) {
		while (true) {
			int end = high - 1, size = high - low;

			/*
			 * Run mixed insertion sort on small non-leftmost parts.
			 */
			if (size < MAX_MIXED_INSERTION_SORT_SIZE + bits && (bits & 1) > 0) {
				mixedInsertionSort(a, low, high - 3 * ((size >> 5) << 3), high);
				return;
			}

			/*
			 * Invoke insertion sort on small leftmost part.
			 */
			if (size < MAX_INSERTION_SORT_SIZE) {
				insertionSort(a, low, high);
				return;
			}

			/*
			 * Check if the whole array or large non-leftmost parts are nearly sorted and
			 * then merge runs.
			 */
			if ((bits == 0 || size > MIN_TRY_MERGE_SIZE && (bits & 1) > 0) && tryMergeRuns(a, low, size)) {
				return;
			}

			/*
			 * Switch to heap sort if execution time is becoming quadratic.
			 */
			if ((bits += DELTA) > MAX_RECURSION_DEPTH) {
				heapSort(a, low, high);
				return;
			}

			/*
			 * Use an inexpensive approximation of the golden ratio to select five sample
			 * elements and determine pivots.
			 */
			int step = (size >> 3) * 3 + 3;

			/*
			 * Five elements around (and including) the central element will be used for
			 * pivot selection as described below. The unequal choice of spacing these
			 * elements was empirically determined to work well on a wide variety of inputs.
			 */
			int e1 = low + step;
			int e5 = end - step;
			int e3 = (e1 + e5) >>> 1;
			int e2 = (e1 + e3) >>> 1;
			int e4 = (e3 + e5) >>> 1;
			int a3 = a[e3];

			/*
			 * Sort these elements in place by the combination of 4-element sorting network
			 * and insertion sort.
			 *
			 * 5 ------o-----------o------------ | | 4 ------|-----o-----o-----o------ | | |
			 * 2 ------o-----|-----o-----o------ | | 1 ------------o-----o------------
			 */
			if (compare(a[e5], a[e2]) < 0) {
				int t = a[e5];
				a[e5] = a[e2];
				a[e2] = t;
			}
			if (compare(a[e4], a[e1]) < 0) {
				int t = a[e4];
				a[e4] = a[e1];
				a[e1] = t;
			}
			if (compare(a[e5], a[e4]) < 0) {
				int t = a[e5];
				a[e5] = a[e4];
				a[e4] = t;
			}
			if (compare(a[e2], a[e1]) < 0) {
				int t = a[e2];
				a[e2] = a[e1];
				a[e1] = t;
			}
			if (compare(a[e4], a[e2]) < 0) {
				int t = a[e4];
				a[e4] = a[e2];
				a[e2] = t;
			}

			if (compare(a3, a[e2]) < 0) {
				if (compare(a3, a[e1]) < 0) {
					a[e3] = a[e2];
					a[e2] = a[e1];
					a[e1] = a3;
				} else {
					a[e3] = a[e2];
					a[e2] = a3;
				}
			} else if (compare(a3, a[e4]) > 0) {
				if (compare(a3, a[e5]) > 0) {
					a[e3] = a[e4];
					a[e4] = a[e5];
					a[e5] = a3;
				} else {
					a[e3] = a[e4];
					a[e4] = a3;
				}
			}

			// Pointers
			int lower = low; // The index of the last element of the left part
			int upper = end; // The index of the first element of the right part

			/*
			 * Partitioning with 2 pivots in case of different elements.
			 */
			if (compare(a[e1], a[e2]) < 0 && compare(a[e2], a[e3]) < 0 && compare(a[e3], a[e4]) < 0
					&& compare(a[e4], a[e5]) < 0) {

				/*
				 * Use the first and fifth of the five sorted elements as the pivots. These
				 * values are inexpensive approximation of tertiles. Note, that pivot1 < pivot2.
				 */
				int pivot1 = a[e1];
				int pivot2 = a[e5];

				/*
				 * The first and the last elements to be sorted are moved to the locations
				 * formerly occupied by the pivots. When partitioning is completed, the pivots
				 * are swapped back into their final positions, and excluded from the next
				 * subsequent sorting.
				 */
				a[e1] = a[lower];
				a[e5] = a[upper];

				/*
				 * Skip elements, which are less or greater than the pivots.
				 */
				while (compare(a[++lower], pivot1) < 0)
					;
				while (compare(a[--upper], pivot2) > 0)
					;

				/*
				 * Backward 3-interval partitioning
				 *
				 * left part central part right part
				 * +------------------------------------------------------------+ | < pivot1 | ?
				 * | pivot1 <= && <= pivot2 | > pivot2 |
				 * +------------------------------------------------------------+ ^ ^ ^ | | |
				 * lower k upper
				 *
				 * Invariants:
				 *
				 * all in (low, lower] < pivot1 pivot1 <= all in (k, upper) <= pivot2 all in
				 * [upper, end) > pivot2
				 *
				 * Pointer k is the last index of ?-part
				 */
				for (int unused = --lower, k = ++upper; --k > lower;) {
					int ak = a[k];

					if (compare(ak, pivot1) < 0) { // Move a[k] to the left side
						while (lower < k) {
							if (compare(a[++lower], pivot1) >= 0) {
								if (compare(a[lower], pivot2) > 0) {
									a[k] = a[--upper];
									a[upper] = a[lower];
								} else {
									a[k] = a[lower];
								}
								a[lower] = ak;
								break;
							}
						}
					} else if (compare(ak, pivot2) > 0) { // Move a[k] to the right side
						a[k] = a[--upper];
						a[upper] = ak;
					}
				}

				/*
				 * Swap the pivots into their final positions.
				 */
				a[low] = a[lower];
				a[lower] = pivot1;
				a[end] = a[upper];
				a[upper] = pivot2;

				/*
				 * Sort non-left parts recursively (possibly in parallel), excluding known
				 * pivots.
				 */
				sort(a, bits | 1, lower + 1, upper);
				sort(a, bits | 1, upper + 1, high);

			} else { // Use single pivot in case of many equal elements

				/*
				 * Use the third of the five sorted elements as the pivot. This value is
				 * inexpensive approximation of the median.
				 */
				int pivot = a[e3];

				/*
				 * The first element to be sorted is moved to the location formerly occupied by
				 * the pivot. After completion of partitioning the pivot is swapped back into
				 * its final position, and excluded from the next subsequent sorting.
				 */
				a[e3] = a[lower];

				/*
				 * Traditional 3-way (Dutch National Flag) partitioning
				 *
				 * left part central part right part
				 * +------------------------------------------------------+ | < pivot | ? | ==
				 * pivot | > pivot | +------------------------------------------------------+ ^
				 * ^ ^ | | | lower k upper
				 *
				 * Invariants:
				 *
				 * all in (low, lower] < pivot all in (k, upper) == pivot all in [upper, end] >
				 * pivot
				 *
				 * Pointer k is the last index of ?-part
				 */
				for (int k = ++upper; --k > lower;) {
					int ak = a[k];

					if (compare(ak, pivot) != 0) {
						a[k] = pivot;

						if (compare(ak, pivot) < 0) { // Move a[k] to the left side
							while (compare(a[++lower], pivot) < 0)
								;

							if (compare(a[lower], pivot) > 0) {
								a[k] = a[--upper]; // wasn't there
								a[upper] = a[lower]; // used to be --upper
							} else {
								a[k] = a[lower]; // wasn't there
							}
							a[lower] = ak;
						} else { // ak > pivot - Move a[k] to the right side
							a[k] = a[--upper]; // wasn't there
							a[upper] = ak; // usd to be --upper
						}
					}
				}

				/*
				 * Swap the pivot into its final position.
				 */
				a[low] = a[lower];
				a[lower] = pivot;

				/*
				 * Sort the right part (possibly in parallel), excluding known pivot. All
				 * elements from the central part are equal and therefore already sorted.
				 */
				sort(a, bits | 1, upper, high);
			}
			high = lower; // Iterate along the left part
		}
	}

	/**
	 * Sorts the specified range of the array using mixed insertion sort.
	 *
	 * Mixed insertion sort is combination of simple insertion sort, pin insertion
	 * sort and pair insertion sort.
	 *
	 * In the context of Dual-Pivot Quicksort, the pivot element from the left part
	 * plays the role of sentinel, because it is less than any elements from the
	 * given part. Therefore, expensive check of the left range can be skipped on
	 * each iteration unless it is the leftmost call.
	 *
	 * @param a    the array to be sorted
	 * @param low  the index of the first element, inclusive, to be sorted
	 * @param end  the index of the last element for simple insertion sort
	 * @param high the index of the last element, exclusive, to be sorted
	 */
	private void mixedInsertionSort(int[] a, int low, int end, int high) {
		if (end == high) {

			/*
			 * Invoke simple insertion sort on tiny array.
			 */
			for (int i; ++low < end;) {
				int ai = a[i = low];

				while (compare(ai, a[--i]) < 0) {
					a[i + 1] = a[i];
				}
				a[i + 1] = ai;
			}
		} else {

			/*
			 * Start with pin insertion sort on small part.
			 *
			 * Pin insertion sort is extended simple insertion sort. The main idea of this
			 * sort is to put elements larger than an element called pin to the end of array
			 * (the proper area for such elements). It avoids expensive movements of these
			 * elements through the whole array.
			 */
			int pin = a[end];

			for (int i, p = high; ++low < end;) {
				int ai = a[i = low];

				if (compare(ai, a[i - 1]) < 0) { // Small element

					/*
					 * Insert small element into sorted part.
					 */
					a[i] = a[--i];

					while (compare(ai, a[--i]) < 0) {
						a[i + 1] = a[i];
					}
					a[i + 1] = ai;

				} else if (p > i && compare(ai, pin) > 0) { // Large element

					/*
					 * Find element smaller than pin.
					 */
					while (compare(a[--p], pin) > 0)
						;

					/*
					 * Swap it with large element.
					 */
					if (p > i) {
						ai = a[p];
						a[p] = a[i];
					}

					/*
					 * Insert small element into sorted part.
					 */
					while (compare(ai, a[--i]) < 0) {
						a[i + 1] = a[i];
					}
					a[i + 1] = ai;
				}
			}

			/*
			 * Continue with pair insertion sort on remain part.
			 */
			for (int i; low < high; ++low) {
				int a1 = a[i = low], a2 = a[++low];

				/*
				 * Insert two elements per iteration: at first, insert the larger element and
				 * then insert the smaller element, but from the position where the larger
				 * element was inserted.
				 */
				if (compare(a1, a2) > 0) {

					while (compare(a1, a[--i]) < 0) {
						a[i + 2] = a[i];
					}
					a[++i + 1] = a1;

					while (compare(a2, a[--i]) < 0) {
						a[i + 1] = a[i];
					}
					a[i + 1] = a2;

				} else if (compare(a1, a[i - 1]) < 0) {

					while (compare(a2, a[--i]) < 0) {
						a[i + 2] = a[i];
					}
					a[++i + 1] = a2;

					while (compare(a1, a[--i]) < 0) {
						a[i + 1] = a[i];
					}
					a[i + 1] = a1;
				}
			}
		}
	}

	/**
	 * Sorts the specified range of the array using insertion sort.
	 *
	 * @param a    the array to be sorted
	 * @param low  the index of the first element, inclusive, to be sorted
	 * @param high the index of the last element, exclusive, to be sorted
	 */
	private void insertionSort(int[] a, int low, int high) {
		for (int i, k = low; ++k < high;) {
			int ai = a[i = k];

			if (compare(ai, a[i - 1]) < 0) {
				while (--i >= low && compare(ai, a[i]) < 0) {
					a[i + 1] = a[i];
				}
				a[i + 1] = ai;
			}
		}
	}

	/**
	 * Sorts the specified range of the array using heap sort.
	 *
	 * @param a    the array to be sorted
	 * @param low  the index of the first element, inclusive, to be sorted
	 * @param high the index of the last element, exclusive, to be sorted
	 */
	private void heapSort(int[] a, int low, int high) {
		for (int k = (low + high) >>> 1; k > low;) {
			pushDown(a, --k, a[k], low, high);
		}
		while (--high > low) {
			int max = a[low];
			pushDown(a, low, a[high], low, high);
			a[high] = max;
		}
	}

	/**
	 * Pushes specified element down during heap sort.
	 *
	 * @param a     the given array
	 * @param p     the start index
	 * @param value the given element
	 * @param low   the index of the first element, inclusive, to be sorted
	 * @param high  the index of the last element, exclusive, to be sorted
	 */
	private void pushDown(int[] a, int p, int value, int low, int high) {
		for (int k;; a[p] = a[p = k]) {
			k = (p << 1) - low + 2; // Index of the right child

			if (k > high) {
				break;
			}
			if (k == high || compare(a[k], a[k - 1]) < 0) {
				--k;
			}
			if (compare(a[k], value) <= 0) {
				break;
			}
		}
		a[p] = value;
	}

	/**
	 * Tries to sort the specified range of the array.
	 *
	 * @param a    the array to be sorted
	 * @param low  the index of the first element to be sorted
	 * @param size the array size
	 * @return true if finally sorted, false otherwise
	 */
	private boolean tryMergeRuns(int[] a, int low, int size) {

		/*
		 * The run array is constructed only if initial runs are long enough to
		 * continue, run[i] then holds start index of the i-th sequence of elements in
		 * non-descending order.
		 */
		int[] run = null;
		int high = low + size;
		int count = 1, last = low;

		/*
		 * Identify all possible runs.
		 */
		for (int k = low + 1; k < high;) {

			/*
			 * Find the end index of the current run.
			 */
			if (compare(a[k - 1], a[k]) < 0) {

				// Identify ascending sequence
				while (++k < high && compare(a[k - 1], a[k]) <= 0)
					;

			} else if (compare(a[k - 1], a[k]) > 0) {

				// Identify descending sequence
				while (++k < high && compare(a[k - 1], a[k]) >= 0)
					;

				// Reverse into ascending order
				for (int i = last - 1, j = k; ++i < --j && compare(a[i], a[j]) > 0;) {
					int ai = a[i];
					a[i] = a[j];
					a[j] = ai;
				}
			} else { // Identify constant sequence
				for (int ak = a[k]; ++k < high && compare(ak, a[k]) == 0;)
					;

				if (k < high) {
					continue;
				}
			}

			/*
			 * Check special cases.
			 */
			if (run == null) {
				if (k == high) {

					/*
					 * The array is monotonous sequence, and therefore already sorted.
					 */
					return true;
				}

				if (k - low < MIN_FIRST_RUN_SIZE) {

					/*
					 * The first run is too small to proceed with scanning.
					 */
					return false;
				}

				run = new int[((size >> 10) | 0x7F) & 0x3FF];
				run[0] = low;

			} else if (compare(a[last - 1], a[last]) > 0) {

				if (count > (k - low) >> MIN_FIRST_RUNS_FACTOR) {

					/*
					 * The first runs are not long enough to continue scanning.
					 */
					return false;
				}

				if (++count == MAX_RUN_CAPACITY) {

					/*
					 * Array is not highly structured.
					 */
					return false;
				}

				if (count == run.length) {

					/*
					 * Increase capacity of index array.
					 */
					run = Arrays.copyOf(run, count << 1);
				}
			}
			run[count] = (last = k);
		}

		/*
		 * Merge runs of highly structured array.
		 */
		if (count > 1) {
			int[] b;
			int offset = low;

			b = new int[size];
			mergeRuns(a, b, offset, 1, run, 0, count);
		}
		return true;
	}

	/**
	 * Merges the specified runs.
	 *
	 * @param a      the source array
	 * @param b      the temporary buffer used in merging
	 * @param offset the start index in the source, inclusive
	 * @param aim    specifies merging: to source ( > 0), buffer ( < 0) or any ( ==
	 *               0)
	 * @param run    the start indexes of the runs, inclusive
	 * @param lo     the start index of the first run, inclusive
	 * @param hi     the start index of the last run, inclusive
	 * @return the destination where runs are merged
	 */
	private int[] mergeRuns(int[] a, int[] b, int offset, int aim, int[] run, int lo, int hi) {
		if (hi - lo == 1) {
			if (aim >= 0) {
				return a;
			}
			for (int i = run[hi], j = i - offset, low = run[lo]; i > low; b[--j] = a[--i])
				;
			return b;
		}

		/*
		 * Split into approximately equal parts.
		 */
		int mi = lo, rmi = (run[lo] + run[hi]) >>> 1;
		while (run[++mi + 1] <= rmi)
			;

		/*
		 * Merge the left and right parts.
		 */
		int[] a1, a2;

		a1 = mergeRuns(a, b, offset, -aim, run, lo, mi);
		a2 = mergeRuns(a, b, offset, 0, run, mi, hi);

		int[] dst = a1 == a ? b : a;

		int k = a1 == a ? run[lo] - offset : run[lo];
		int lo1 = a1 == b ? run[lo] - offset : run[lo];
		int hi1 = a1 == b ? run[mi] - offset : run[mi];
		int lo2 = a2 == b ? run[mi] - offset : run[mi];
		int hi2 = a2 == b ? run[hi] - offset : run[hi];

		mergeParts(dst, k, a1, lo1, hi1, a2, lo2, hi2);
		return dst;
	}

	/**
	 * Merges the sorted parts.
	 *
	 * @param merger parallel context
	 * @param dst    the destination where parts are merged
	 * @param k      the start index of the destination, inclusive
	 * @param a1     the first part
	 * @param lo1    the start index of the first part, inclusive
	 * @param hi1    the end index of the first part, exclusive
	 * @param a2     the second part
	 * @param lo2    the start index of the second part, inclusive
	 * @param hi2    the end index of the second part, exclusive
	 */
	private void mergeParts(int[] dst, int k, int[] a1, int lo1, int hi1, int[] a2, int lo2, int hi2) {
		/*
		 * Merge small parts sequentially.
		 */
		while (lo1 < hi1 && lo2 < hi2) {
			dst[k++] = compare(a1[lo1], a2[lo2]) < 0 ? a1[lo1++] : a2[lo2++];
		}
		if (dst != a1 || k < lo1) {
			while (lo1 < hi1) {
				dst[k++] = a1[lo1++];
			}
		}
		if (dst != a2 || k < lo2) {
			while (lo2 < hi2) {
				dst[k++] = a2[lo2++];
			}
		}
	}
}
