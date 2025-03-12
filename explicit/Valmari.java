package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.jas.structure.Value;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;

public class Valmari<Value> extends AbstractBisimulation<Value> {

    public Valmari(PrismComponent parent) throws PrismException {
        super(parent);
    }
    
    public static final double ACCURACY = 1E-8;
    public static final int PRECISION = 3;

    private int[] elems;
    private int[] location;
    private int[] block;
    private ArrayList<Integer> start;
    private ArrayList<Integer> end;
    private ArrayList<Integer> borderline;

    public int[] evaluate(DTMCSimple<Value> dtmc, List<BitSet> propBSs) {
        Evaluator<Value> eval = dtmc.getEvaluator();
        this.elems = new int[numStates];
        this.location = new int[numStates];
        this.block = partition.clone();
        this.start = new ArrayList<>(numBlocks);
        this.end = new ArrayList<>(numBlocks);
        this.borderline = new ArrayList<>(numBlocks);
    
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
    
        HashSet<Integer> UB = new HashSet<>();
        for (int i = 0; i < this.start.size(); i++) {
            UB.add(i);
        }
        ArrayList<Integer> BT = new ArrayList<>(10);
        ArrayList<Integer> ST = new ArrayList<>(10);
        double[] w = new double[numStates];
    
        while (!UB.isEmpty()) {
            int splitter = UB.iterator().next();
            UB.remove(splitter);
            ST.clear();
            
            for (int i = this.start.get(splitter); i < this.end.get(splitter); i++) {
                for (int j = 0; j < numStates; j++) {
                    double prob = eval.toDouble(dtmc.getProbability(j, this.elems[i]));
                    if (prob > 0) {
                        if (w[j] == 0) {
                            ST.add(j);
                            w[j] = prob;
                        } else {
                            w[j] += prob;
                        }
                    }
                }
            }
            
            for (int j : ST) {
                int b = this.block[j];
                if (this.start.get(b) == this.borderline.get(b)) {
                    BT.add(b);
                }
                this.mark(j, b);
            }
            
            for (int b : BT) {
                int b1;
                int l = 0;
                if (this.borderline.get(b) == this.end.get(b)) {
                    this.borderline.set(b, this.start.get(b));
                    b1 = b;
                } else {
                    b1 = this.start.size();
                    this.split(b);
                    l++;
                }
                
                double pmc = this.pmc(w, b1);
                for (int i = this.start.get(b1); i < this.end.get(b1); i++) {
                    if (!isEqual(w[elems[i]], pmc)) {
                        mark(elems[i], b1);
                    }
                }
                
                if (this.borderline.get(b1) != this.start.get(b1)) {
                    int b2 = this.start.size();
                    this.split(b1);
                    l++;
                    new ArraysSort(w).sort(this.elems, this.start.get(b2), this.end.get(b2));
                }
            }
            
            BT.clear();
            for (int j : ST) {
                w[j] = 0;
            }
        }
        return this.block;
    }
    
    private class ArraysSort {
        private final double[] valueArray;

        public ArraysSort(double[] array) {
            this.valueArray = array;
        }

        public void sort(int[] a, int fromIndex, int toIndex) {
            if (fromIndex > toIndex) {
                throw new IllegalArgumentException("Invalid range: fromIndex > toIndex");
            }
            quickSort(a, fromIndex, toIndex - 1);
        }

        private void quickSort(int[] a, int low, int high) {
            if (low < high) {
                int pi = partition(a, low, high);
                quickSort(a, low, pi - 1);
                quickSort(a, pi + 1, high);
            }
        }

        private int partition(int[] a, int low, int high) {
            int pivot = a[high];
            int i = (low - 1);
            for (int j = low; j < high; j++) {
                if (compare(a[j], pivot) <= 0) {
                    i++;
                    int temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                }
            }
            int temp = a[i + 1];
            a[i + 1] = a[high];
            a[high] = temp;
            return i + 1;
        }

        private int compare(int index1, int index2) {
            return Double.compare(this.valueArray[index1], this.valueArray[index2]);
        }
    }
}
