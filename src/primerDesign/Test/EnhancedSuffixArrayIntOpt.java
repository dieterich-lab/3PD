package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Stack;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.dsc.FastStack;
import primerDesign.util.SimpleTimer;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;


/**
 * This class implements the enhanced suffix array as proposed by Kurtz et.al. 2004
 * 
 * Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 * 
 * @author froehler
 *
 */
public class EnhancedSuffixArrayIntOpt implements DNASequenceIndex, Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private char[] sequence;
	private boolean includesScanRegion;
	private int maxWordSize;
	private int sequenceLength;
	private static String TERMINATION_SYMBOL = "$";
	private int[] suftab;
	private byte[] lcptab;
	private byte[] childtab;
	private IntArrayList lcpmap;
	private IntArrayList childmap;
	private OpenIntIntHashMap bucketTab;
	private int d; // length of prefixes indexed in bucket table
	private static boolean printStatus = false;
	
	public EnhancedSuffixArrayIntOpt(String sequence){
		if(sequence.length() < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		this.sequence = sequence.toCharArray();
		this.sequenceLength = sequence.length();
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Creates an enhanced suffix array for the string 'sequence'.
	 * @param sequence the string to construct the enhanced suffi array for
	 */
	public void createIndex(int maxWordSize, boolean includesScanRegion){

		this.maxWordSize = sequenceLength;
		this.includesScanRegion = includesScanRegion;
		this.suftab = new int[sequenceLength+1];
		this.lcptab = new byte[sequenceLength+1];
		this.lcpmap = new IntArrayList();
		this.childtab = new byte[sequenceLength+1];
		this.childmap = new IntArrayList();
		bucketTab = new OpenIntIntHashMap();
		d = 8;
		
		for(int i=0; i< sequenceLength; i++){
			this.suftab[i] = i;
		}
		this.suftab[sequenceLength] = Integer.MAX_VALUE;
		// sort table
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		SimpleTimer timer = new SimpleTimer();
		if(printStatus){
			System.gc();
			System.out.print("Before sort: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		quicksortSuffices(0, sequenceLength-1);
		if(printStatus){
			System.out.println(" - sorted in " + timer.getTimeString());
			System.gc();
			System.out.println("Before LCP Table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		computeLCPTable();
		if(printStatus){
			System.out.println(" - lcp table computed in " + timer.getTimeString());
			System.gc();
			System.out.print("Before child table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		computeChildTable();
		if(printStatus){
			System.out.println(" - child table computed in " + timer.getTimeString());
			System.gc();
			System.out.print("Before bucket table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		computeBucketTable();
		if(printStatus){
			System.out.println(" - bucket table computed in " + timer.getTimeString());
			System.gc();
			System.out.println("After all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Computes the longest common prefix of a sorted suffix array for element 'i' w.r.t element 'i-1'.
	 *
	 */
	private void computeLCPTable(){
		for(int i=0; i<this.lcptab.length; i++) this.lcptab[i] = Byte.MIN_VALUE;
		int lcp;
		int stop = this.suftab.length;
		for (int i = 0; i < stop; i++) {
			if(i == 0 || i == stop-1) setLcpTab(i, 0);
			else{
				lcp = computeLcp(i, i-1);
				setLcpTab(i, lcp);
//				this.lcptab_compara[i] = lcp;
			}
		}
		if(printStatus){
			System.out.println("Size of LCP Map: " + this.lcpmap.size());
			System.out.println("Size of LCP table: " + this.lcptab.length);
		}
	}
	
	private void computeChildTable(){
		Arrays.fill(this.childtab, Byte.MIN_VALUE);
//		computeChildNextIndex();
//		computeChildUpDown();
		computeCombinedChildTable();
		quicksortChildmap(0, this.childmap.size()-1);
		System.out.println("Size of childmap: " + this.childmap.size());
		System.out.println("Size of childtab: " + this.childtab.length);
	}
	
	protected void computeChildUpDown(){
		int lastIndex = -1;
		Stack<Integer> stack = new Stack<Integer>();
		stack.push(0);
		int top;
		for (int i = 1; i <= sequenceLength; i++) {
			while(getLcpTab(i) < getLcpTab(stack.peek())){
				lastIndex = stack.pop();
				top = stack.peek();
				if(getLcpTab(i) <= getLcpTab(top) && getLcpTab(top) != getLcpTab(lastIndex) && !containsNextIndex(top)){
					setChildTab(top, lastIndex);
					quicksortChildmap(0, this.childmap.size()-1);
				}
			}
			if(lastIndex != -1 && !containsNextIndex(i-1) && !containsDownIndex(i-1)){
				setChildTab(i-1, lastIndex);
				quicksortChildmap(0, this.childmap.size()-1);
				lastIndex = -1;
			}	
			stack.push(i);
		}
	}
	
	protected void computeChildNextIndex(){
		int lastIndex;
		Stack<Integer> stack = new Stack<Integer>();
		stack.push(0);
		for (int i = 1; i <= sequenceLength; i++) {
			while(getLcpTab(i) < getLcpTab(stack.peek())){
				stack.pop();
			}
			if(getLcpTab(i) == getLcpTab(stack.peek())){
				lastIndex = stack.pop();
				setChildTab(lastIndex, i);
			}
			stack.push(i);
		}
	}
	
	protected void computeCombinedChildTable(){
		int lastIndex = -1;
		//Stack<Integer> stack = new Stack<Integer>();
		FastStack stack = new FastStack();
		stack.push(0);
		int top;
		for (int i = 1; i <= sequenceLength; i++) {
			while(getLcpTab(i) < getLcpTab(stack.peek())){
				lastIndex = stack.pop();
				top = stack.peek();
				if(getLcpTab(i) < getLcpTab(top) && getLcpTab(top) != getLcpTab(lastIndex)){
					setChildTab(top, lastIndex);
//					this.childtab_compara_down[top] = lastIndex;
				}
			}

			if(lastIndex != -1 && getLcpTab(i-1) > getLcpTab(i)){
				setChildTab(i-1, lastIndex);
//				this.childtab_compara_up[i] = lastIndex;
				lastIndex = -1;
			}
			if(getLcpTab(i) == getLcpTab(stack.peek()) && stack.peek() < i){
				lastIndex = stack.peek();
				setChildTab(lastIndex, i);
//				this.childtab_compara_next[lastIndex] = i;
			}
			stack.push(i);
		}
	}
	
	protected byte[] getChildTab(){
		return this.childtab;
	}
	
	protected void initChildTab(){
		this.childtab = new byte[sequenceLength+1];
		this.childmap = new IntArrayList();
		Arrays.fill(this.childtab, Byte.MIN_VALUE);
	}
	
	protected void setChildTab(byte[] tab){
		this.childtab = tab;
	}
	
	/**
	 * Computes the bucket table of prefix->first index associations.
	 * 
	 * The bucket table contains mappings of the form:
	 * pattern -> start index of interval in suffix array, where pattern is one prefix
	 * the hashCode of the pattern is stored in the bucket table instead of its string value!
	 *
	 */
	private void computeBucketTable(){
		//String currentString;
		int suftab;
		for(int i=sequenceLength-1; i>=0; i--){
			suftab = getSufTab(i);
			if(sequenceLength-suftab < d) continue;
			else if(getLcpTab(i+1) >= d){
//				currentString = getSubstring(suftab, suftab+d);
//				bucketTab.put(currentString.hashCode(), i);
				bucketTab.put(getHashCode(suftab, suftab+d), i);
			}
		}
	}
	
	/**
	 * Computes the hash code of sequence[i]...sequence[j-1] as String.hashCode().
	 * 
	 * @param i the start index
	 * @param j the stop index +1 s.th. stop-start=length
	 * @return the hash code of sequence[i]...sequence[j-1] as String.hashCode()
	 */
	private int getHashCode(int i, int j){
		assert i>=0 && j>i && i<=this.sequence.length && j<=this.sequence.length;
		int result = 0;
		int length = j-i;
		for(int k=0; k<length; k++){
			//result += sequence[i+k] * 31^(length-1-k);
			result = 31*result + this.sequence[i+k];
		}
		return result;
	}
	
	/**
	 * Checks whether a pattern can be found in the bucket table.
	 * 
	 * @param pattern the pattern to search
	 * @return true iff the pattern can be found in the bucket table
	 */
	private boolean isInBucketTab(String pattern){
		assert pattern.length() > 0 && pattern != null;
		return bucketTab.containsKey(pattern.hashCode());
	}
	
	/**
	 * Returns the start index of an interval where pattern is one prefix.
	 * 
	 * @param pattern the pattern
	 * @return the start index of an interval where pattern is one prefix
	 */
	private int getFromBucketTab(String pattern){
		assert pattern.length() > 0 && pattern != null;
		if(isInBucketTab(pattern)) return bucketTab.get(pattern.hashCode());
		else return -1;
	}
	
	/**
	 * Returns an interval where pattern is one prefix.
	 * 
	 * @param pattern the pattern
	 * @return an int[start,stop] interval where pattern is one prefix
	 */
	private int[] getIntervalFromBucketTable(String pattern){
		assert pattern.length() >= d && pattern != null;
		pattern = pattern.substring(0, d);
		if(!isInBucketTab(pattern)) return new int[]{0,-1}; // pattern of length 'd' was not put into bucket table -> has NO associated interval in the interval tree!
		int start = this.bucketTab.get(pattern.hashCode());
		int stop = start -1; 
		for(int i=start; i<sequenceLength; i++){
			if(getLcpTab(i+1) < d){
				stop = i;
				break;
			}else stop = i;
		}
		return new int[]{start, stop};
	}
	
	/**
	 * Convenience function to return the entry 'SUF_TAB' for index i.
	 * 
	 * @param i the index
	 * @return the entry 'SUF_TAB' for index i
	 */
	protected int getSufTab(int i){
		assert(i>=0 && i<this.suftab.length);
		return this.suftab[i];
	}
	
	private void setSufTab(int i, int j){
		assert(i>=0 && j>=0 && i<this.suftab.length && j<this.suftab.length);
		this.suftab[i] = j;
	}
	
	/**
	 * Convenience function to return the entry 'LCP_TABLE' for index 'i'.
	 * 
	 * @param i the index
	 * @return the entry 'LCP_TABLE' for index 'i'
	 */
	protected int getLcpTab(int i){
		assert(i>=0 && i<this.lcptab.length);
		return (this.lcptab[i] < Byte.MAX_VALUE) ? this.lcptab[i] - Byte.MIN_VALUE : getMap(this.lcpmap, i, 0, this.lcpmap.size()-1);
	}
	
	protected void setLcpTab(int i, int j){
		assert(i>=0 && j>=0 && i<this.lcptab.length && j<=this.sequenceLength);
		if(j + Byte.MIN_VALUE < Byte.MAX_VALUE) this.lcptab[i] = (byte) (j + Byte.MIN_VALUE);
		else{
			this.lcptab[i] = Byte.MAX_VALUE;
			putLcpMap(i, j);  // since computeLCPTable calls setLcpTable only one per index i, no 'contains' lookup needs to be done in lcpmap!!
		}
	}
	
	/**
	 * Stores the real value of lcp(i,i-1)=j.
	 * 
	 * Format: [index, value]_n -> position%2=0 -> index, position%2=1 -> value
	 * 
	 * @param i the index of the lcp value
	 * @param j the lcp value for index i
	 */
	protected void putLcpMap(int i, int j){
		assert i>=0 && j>=0 && i<sequenceLength && j<=this.sequenceLength;
		this.lcpmap.add(i);
		this.lcpmap.add(j);
	}
	
	/**
	 * Binary searches the real lcp value for index i from lcptab.
	 * 
	 * This function is only called iff the corresponding child tab index for 'i' is Byte.MAX_VALUE
	 * and therefore an entry in the childmap table has been created previously.
	 * 
	 * @param i the index from the lcp tab where lcptab(i) = Byte.MAX_VALUE
	 * @param from the start position for the binary search
	 * @param to the end position of the binary search
	 * @return the real lcp value for index i from lcptab
	 */
	protected int getLcpMap(int i, int from, int to){
		assert(i>=0 && from>=0 && to>=0 && from<=to && i<sequenceLength && from<this.lcpmap.size() && to<this.lcpmap.size());
		int result = 0;
		int mean = Math.round((0.0f+to-from)/2) + from;
		int element;
		if(mean%2 == 1){
			element = this.lcpmap.getQuick(--mean);
		}else{
			element = this.lcpmap.getQuick(mean);
		}
		if(element == i) return this.lcpmap.getQuick(mean+1);
		else if(element < i && mean+2<=to) result =  getLcpMap(i, mean+2, to);
		else if(element > i && mean-2>=from) result = getLcpMap(i, from, mean-2);
		else throw new IllegalStateException("Unhandled case!");
		return result; // never returned!
	}
	
	protected int getMap(IntArrayList map, int i, int from, int to){
		assert(i>=0 && from>=0 && to>=0 && from<=to && i<sequenceLength && from<map.size() && to<map.size());
		int result = 0;
		int mean = Math.round((0.0f+to-from)/2) + from;
		int element;
		if(mean%2 == 1){
			element = map.getQuick(--mean);
		}else{
			element = map.getQuick(mean);
		}
		if(element == i) return map.getQuick(mean+1);
		else if(element < i && mean+2<=to) result =  getMap(map, i, mean+2, to);
		else if(element > i && mean-2>=from) result = getMap(map, i, from, mean-2);
		else throw new IllegalStateException("Unhandled case!");
		return result; // never returned!
	}
	
	/**
	 * Convenience function to return the entry 'CHILD_TABLE' for index 'i'.
	 * 
	 * @param i the index
	 * @return the entry 'CHILD_TABLE' for index 'i', -1 iff it needs to be recomputed
	 */
	protected int getChildTab(int i){
		assert(i>=0 && i<this.childtab.length);
		return (this.childtab[i] < Byte.MAX_VALUE) ? this.childtab[i] - Byte.MIN_VALUE + i : getMap(this.childmap, i, 0, this.childmap.size()-1); 
	}
	
	private void setChildTab(int i, int j){
		assert i>=0 && i<this.suftab.length && j>=0 && j<this.suftab.length : "i and j must be in interval [0,sequenceLength]";
		if(j-i + Byte.MIN_VALUE < Byte.MAX_VALUE && j-i >= 0) this.childtab[i] = (byte) (j-i+Byte.MIN_VALUE);
		else{
			this.childtab[i] = Byte.MAX_VALUE;
			putChildMap(i,j);
		}
	}
	
	private void putChildMap(int i, int j){
		assert i>=0 && j>=0 && i<sequenceLength && j<=this.sequenceLength;
		this.childmap.add(i);
		this.childmap.add(j);
	}
	
	// for debugging purpose
	protected void setChildMap(IntArrayList map){
		this.childmap = map;
	}
	protected IntArrayList getChildMap(){
		return this.childmap;
	}
	
	protected void quicksortChildmap(int from, int to){
		assert(from>=0 && to>=0 && from<this.childmap.size() && to<this.childmap.size() && from < to) : "Childmap index 'from' has to be smaller than 'to'!";
	    int i = (from%2 == 1) ? from-1 : from;
	    int j = (to%2 == 1) ? to-1 : to;
	    int pivot = (from+to)/2;
		if(pivot%2 == 1){
			pivot--;
		}
	    
	    //  partition
		int tempKey;
		int tempVal;
	    do{    
	        while (i < pivot && this.childmap.getQuick(i) < this.childmap.getQuick(pivot)) i+=2; 
	        while (j > pivot && this.childmap.getQuick(j) > this.childmap.getQuick(pivot)) j-=2;
	        if (i<j)
	        {
	        	tempKey = this.childmap.getQuick(i);
	        	tempVal = this.childmap.getQuick(i+1);
	        	this.childmap.setQuick(i, this.childmap.getQuick(j));
	        	this.childmap.setQuick(i+1, this.childmap.getQuick(j+1));
	        	this.childmap.setQuick(j, tempKey);
	        	this.childmap.setQuick(j+1, tempVal);
	        	if(i == pivot){
	        		pivot = j;
	        		i+=2;
	        	}
	        	else if(j == pivot){
	        		pivot = i;
	        		j-=2;
	        	}
	        	else{
	        		i+=2; 
		            j-=2;
	        	}
	        }
	    }while(i<j);

	    //  recursion
	    if (from<pivot-2) quicksortChildmap(from, pivot-1);
	    if (pivot+2<to) quicksortChildmap(pivot+2, to);
	}
	
	public int computeLcp(int i, int j){
		int idxI = this.suftab[i];
		int idxJ = this.suftab[j];
		int end = sequenceLength - Math.max(idxI, idxJ);
		int lcp = 0;
		for(int k = 0; k < end; k++) {
			if(this.sequence[idxI + k] == this.sequence[idxJ + k]) lcp++;
			else break;
		}
		return lcp;
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is an up index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is an up index
	 */
	protected boolean containsUpIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		if(i == sequenceLength) return true;
		else return getLcpTab(i) > getLcpTab(i+1);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is a down index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is a down index
	 */
	protected boolean containsDownIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		return getLcpTab(getChildTab(i)) > getLcpTab(i);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is a next index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is the next index
	 */
	protected boolean containsNextIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		return getLcpTab(getChildTab(i)) == getLcpTab(i) && getChildTab(i) > i; // if no child value can be found during ChildTabNext-recomputation, value '-1' is returned to signal this
	}
	
	/**
	 * Returns a list of child intervals of the interval [i,j]
	 * 
	 * @param i the start position of the 'parent' interval
	 * @param j the stop position of the 'parent' interval
	 * @return a list of child intervals of the 'parent' interval [i,j], format: [DOWN,UP]_x
	 */
	protected int[] getChildIntervals(int i, int j){
		assert i>=0 && j>=0 && i<this.suftab.length && j<this.suftab.length;
		// init intervals list
		int[] result = new int[0]; // list of child intervals: i%2=0 -> lower bound, i%2=1 -> upper bound
		
		if(i>=0 && j<sequenceLength){
			int i1;
			
			if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) i1 = getChildTabUP(j+1);
			else i1 = getChildTabDown(i);
			// add i, i1-1 to intervals list
			result = addInterval(result, i, i1-1);
			
			while(containsNextIndex(i1)){
				int i2 = getChildTabNext(i1);
				// add i1, i2-1 to intervals list
				result = addInterval(result, i1, i2-1);
				i1 = i2;
			}
			// add i1, j to intervals list
			result = addInterval(result, i1, j);
		}else{
			int k=i;
			while(containsNextIndex(k)){
				result = addInterval(result, k, getChildTabNext(k)-1);
				k = getChildTabNext(k);
			}
			if(j>=k) result = addInterval(result, k, j);
		}
		
		return result;
	}
	
	private int[] addInterval(int[] interval, int a, int b){
		int[] temp = new int[interval.length+2];
		System.arraycopy(interval, 0, temp, 0, interval.length);
		temp[temp.length-2] = a;
		temp[temp.length-1] = b;
		
		return temp;
	}
	
//	/**
//	 * For debugging purpose...
//	 * @param intervals
//	 * @return
//	 */
//	protected ObjectArrayList getChildIntervalsValues(ObjectArrayList intervals){
//		assert intervals != null;
//		ObjectArrayList values = new ObjectArrayList();
//		for(int i=0; i<intervals.size(); i++){
//			values.add(new Integer[]{((Interval)intervals.get(i)).getLowerBoundary(),((Interval)intervals.get(i)).getUpperBoundary()});
//		}
//		return values;
//	}
//	
//	protected int[] intervalToArray(Interval interval){
//		assert interval != null;
//		return new int[]{interval.getLowerBoundary(), interval.getUpperBoundary()};
//	}
	
	/**
	 * Returns the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval.
	 * 
	 * @param i the start of the interval
	 * @param j the end of the interval
	 * @param p the character of the child interval to be found at position getLCP(i,j) at all suffices in the child interval
	 * @return the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval
	 */
	protected int[] getChildInterval(int i, int j, char p){
		assert i>=0 && j>=0 && i<this.childtab.length && j<this.childtab.length && p >= Character.MIN_VALUE && p<= Character.MAX_VALUE;
		int[] result = new int[]{0,-1};
		int lcpOffset = getLCP(i, j);
		int[] intervals = getChildIntervals(i, j);
		for(int k=0; k<intervals.length; k+=2){
			int lowerBound = intervals[k]; 
			int index = getSufTab(lowerBound) + lcpOffset;
			if(index < sequenceLength && this.sequence[index] == p){
				result = new int[]{intervals[k], intervals[k+1]};
				break;
			}
		}
		return result;
	}
	
	/**
	 * Returns all match positions for pattern 'pattern' in the current suffix array.
	 * 
	 * @param patter the pattern to search in the suffix array
	 * @return all match positions in the string the suffix array was constructed for
	 */
	public int[] findMatchPositions(String patter){
		if(patter.length() < 1) throw new IllegalArgumentException("Invalid pattern to scan with!");
		int c = 0;
		boolean queryFound = true;
		int elements;
		int[] interval;
		// if pattern length >= d -> use bucket table to locate initial interval
		if(patter.length() >= d){
			interval = getIntervalFromBucketTable(patter);
			c = d-1;
		}
		// else de-novo search
		else 
			interval = getChildInterval(0, sequenceLength, patter.charAt(c));
		
		if(interval[0] <= interval[1]){ // an empty interval has length -1: e.g. index a=0, index b=-1
			int start = interval[0];
			int end = interval[1];
			while(interval[0] <= interval[1]  && c<patter.length() && queryFound){
				if(start != end){
					int l = getLCP(start, end);
					int min = Math.min(l,patter.length());
					queryFound = isEqualToSequRegion(patter, getSufTab(start)+c, c, min-c);
					c = min;
					if(c<patter.length()) interval = getChildInterval(start, end, patter.charAt(c));
					else break;
					if(interval[0] > interval[1]) break;
					else{
						start = interval[0];
						end = interval[1];
					}
				}
				else{
					queryFound = isEqualToSequRegion(patter, getSufTab(start)+c, c, patter.length()-c);
					break;
				}
			}
		}
		if(interval[0] > interval[1]) elements = 0;
		else elements = interval[1]-interval[0]+1;
		int[] positions = new int[elements];
		for(int i=0; i<elements; i++){
			positions[i] = getSufTab(interval[0] + i);
		}
		return positions;
	}
	
	/**
	 * Compares a patter region to a sequence region efficiently.
	 * 
	 * @param pattern the pattern to compare to the index sequence region
	 * @param sstart the start position in the index sequence
	 * @param pstart the start position in the pattern
	 * @param length the number of characters to compare srating at the start positions
	 * @return true iff sequence.substring(sstart,sstart+length).equals(pattern.substring(pstart,pstart+length))
	 */
	private boolean isEqualToSequRegion(String pattern, int sstart, int pstart, int length){
		assert(sstart >= 0 && sstart <= sequenceLength && pstart >= 0 && pstart <= pattern.length() && length >=0 && length <= Math.min(pattern.length()-pstart, sequenceLength-sstart));
		if(sstart == sequenceLength || pstart == pattern.length() || length == 0) return false;
		for(int i=0; i<length; i++){
			if(this.sequence[sstart+i] != pattern.charAt(pstart+i)) return false;
		}
		return true;
	}
	
	/**
	 * Returns the length of the longest common prefix of the interval [i,j].
	 * 
	 * @param i the start of the interval
	 * @param j the stop of the interval
	 * @return the length of the longest common prefix of the interval [i,j]
	 */
	private int getLCP(int i, int j){
		assert i>= 0 && j>=0 && i<=sequenceLength && j<=sequenceLength;
		if(j == sequenceLength) return 0; // there is NO LCP between the whole string and the sentinel!
		else if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) return getLcpTab(getChildTabUP(j+1));
		else return getLcpTab(getChildTabDown(i));
	}
	
	protected int getChildTabUP(int i){
		assert containsUpIndex(i-1) : "Child table does not contain an up value for index " + i + "!";
		return getChildTab(i-1);
	}
	
	protected int getChildTabDown(int i){
		assert containsDownIndex(i) : "Child table does not contain a down value for index " + i + "!";
		return getChildTab(i);
	}
	
	protected int getChildTabNext(int i){
		assert containsNextIndex(i) : "Child table does not contain a next value for index " + i + "!";
		return getChildTab(i);
	}
	
//	private String getStringLabel(int element){
//		assert element>= 0 && element<this.sequenceLength;
//		return this.sequence.substring(element);
//	}
	
	/**
	 * Quicksort for an array of suffices.
	 * 
	 * @param from first index to be included in the sorting
	 * @param to last index to be included in the sorting
	 */
	private void quicksortSuffices(int from, int to){
		assert(from>=0 && to>=0 && from<this.sequenceLength && to<this.sequenceLength && from < to) : "Suffix index 'from' has to be smaller than 'to'!";
	    int i = from;
	    int j = to;
	    //int pivot = Math.round((from+to+0.0f)/2);
	    int pivot = (from+to)/2;
	    //int pivot = from;
	    
	    //  partition
	    do{    
	        while (i < pivot && compareSuffix(i, pivot) <= 0) i++; 
	        while (j > pivot && compareSuffix(j, pivot) >= 0) j--;
	        if (i<j)
	        {
	        	swap(i, j);
	        	if(i == pivot){
	        		pivot = j;
	        		i++;
	        	}
	        	else if(j == pivot){
	        		pivot = i;
	        		j--;
	        	}
	        	else{
	        		i++; 
		            j--;
	        	}
	        }
	    }while(i<j);

	    //  recursion
	    if (from<pivot-1) quicksortSuffices(from, pivot-1);
	    if (pivot+1<to) quicksortSuffices(pivot+1, to);
    }
	
	/**
	 * Swaps two elements (suffices) in the suffix array.
	 * 
	 * @param a the first element
	 * @param b the second element
	 */
	private void swap(int a, int b){
		assert a>=0 && b>=0 && a<this.sequenceLength && b<this.sequenceLength;
		int temp = getSufTab(a);
		setSufTab(a, getSufTab(b));
		setSufTab(b, temp);
	}
	
	/**
	 * Compares two suffices - similar to compareTo for strings.
	 * 
	 * @param a the index of the first suffix in SUF_TAB
	 * @param b the index of the second suffix in SUF_TAB
	 * @return -1|0|1 like compareTo for Strings does
	 */
private int compareSuffix(int a, int b){
		assert a>=0 && b>=0 && a<this.suftab.length && b<this.suftab.length;
		int idxA = getSufTab(a);
		int idxB = getSufTab(b);
		int end = Math.min(sequenceLength-idxA, sequenceLength-idxB);
		for (int i = 0; i < end; i++) {
			if(this.sequence[idxA+i] == this.sequence[idxB+i]) continue;
			else if(this.sequence[idxA+i] > this.sequence[idxB+i]) return 1;
			else if(this.sequence[idxA+i] < this.sequence[idxB+i]) return -1;
		}
		if(sequenceLength-idxA > sequenceLength-idxB) return 1;
		else if(sequenceLength-idxA < sequenceLength-idxB) return -1;
		else return 0;
	}
	
	/**
	 * Returns a textual representation of the current suffix array.
	 *
	 */
	public void printArray(){
		for(int i=0; i<this.suftab.length; i++){
			if(getSufTab(i) == Integer.MAX_VALUE) System.out.print(TERMINATION_SYMBOL);
			else System.out.print(getSubstring(getSufTab(i)) + TERMINATION_SYMBOL);
			System.out.println("\t" + getLcpTab(i) + "\t" + ((containsNextIndex(i)) ? getChildTabNext(i) : ((containsUpIndex(i)) ? getChildTabUP(i) : getChildTabDown(i))));
		}
		System.out.println();
	}
	
	private String getSubstring(int i){
		return new String(this.sequence, i, sequenceLength - i);
	}
	
	private String getSubstring(int i, int j){
		return new String(this.sequence, i, Math.min(sequenceLength - i, j));
	}
	
//	/**
//	 * Encapsulates an interval in the enhanced suffix array.
//	 * 
//	 * @author froehler
//	 *
//	 */
//	
//	private class Interval{
//		private int lowerBoundary = 0;
//		private int upperBoundary = Integer.MAX_VALUE;
//		private boolean isEmptyInterval;
//		
//		public Interval(){
//			this.isEmptyInterval = true;
//		}
//		
//		public Interval(int i, int j){
//			setLowerBoundary(i);
//			setUpperBoundary(j);
//			this.isEmptyInterval = false;
//		}
//		
//		public int getLowerBoundary(){
//			assert !isEmptyInterval : "No lower bound has been set for this interval!";
//			return this.lowerBoundary;
//		}
//		
//		public void setLowerBoundary(int i){
//			assert i >= 0 && i <= sequenceLength && i <= this.upperBoundary : "Interval boundaries must be in interval [0," + (sequenceLength-1) + "]";
//			this.lowerBoundary = i;
//		}
//		
//		public int getUpperBoundary(){
//			assert !isEmptyInterval : "No upper bound has been defined for this interval!";
//			return this.upperBoundary;
//		}
//		
//		public void setUpperBoundary(int i){
//			assert  i >= 0 && i <= sequenceLength && i >= this.lowerBoundary : "Interval boundaries must be in interval [0," + (sequenceLength-1) + "]";
//			this.upperBoundary = i;
//		}
//		
//		public boolean isEmptyInterval(){
//			return this.isEmptyInterval;
//		}
//	}

	public int getMaxWordSize() {
		return this.maxWordSize;
	}

	public String getSequence() {
		return new String(this.sequence);
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	public Integer[] searchMatchPositionsInIndex(String searchString) {
		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
		int[] positions = findMatchPositions(searchString);
		Integer[] result = new Integer[positions.length];
		System.arraycopy(positions, 0, result, 0, positions.length);
		return result;
	}

	public int searchNbMatchesInIndex(String searchString) {
		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
		return findMatchPositions(searchString.toUpperCase()).length;
	}
	
	public void serialize(String filename){
		try{
			ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static EnhancedSuffixArrayIntOpt deserialize(String filename){
		EnhancedSuffixArrayIntOpt result = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			result = (EnhancedSuffixArrayIntOpt) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IOException{		
		EnhancedSuffixArrayIntOpt.printStatus = false;
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			EnhancedSuffixArrayIntOpt testIndex;
			Pattern fastaStart = Pattern.compile(">.*");
				System.out.println("Reading sequence " + file);
				{
					String sequence;
					{
						BufferedReader reader = new BufferedReader(new FileReader(file)); //path + sequences[i]));
						//BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(path + sequences[i]), "US-ASCII"));
						StringBuffer seq = new StringBuffer();
						String line;
						while((line = reader.readLine()) != null){
							if(!fastaStart.matcher(line).matches()){
								seq.append(line.trim());
							}
							else seq.append("N");
						}
						sequence = seq.toString().toUpperCase();
						seq = null;  // clear string buffer -> save 2n space!
					}
					System.gc();
					System.out.println("Reading sequence " + file + " (length: " + sequence.length() + ") took " + timer.getTimeString());
					
					System.out.println("Creating EnhancedSuffixArray for sequence " + file);
					
					System.out.println("Constructing EnhancedSuffixArray");
					//sequence = "ACAAACATAT";
					testIndex = new EnhancedSuffixArrayIntOpt(sequence);
					sequence = null; // clear sequence -> save 2n space!
				}
				//System.gc();
				testIndex.createIndex(Integer.MAX_VALUE, true);
				//testIndex = EnhancedSuffixArrayInt.deserialize("EnhancedSuffixArrayInt.ser");
				System.out.println("EnhancedSuffixArray construction complete - " + timer.getTimeString());
	//			String query = "ATGATGATG";
	//			System.out.println("Querying EnhancedSuffixArray for: " + query);
	//			System.out.println(testIndex.searchNbMatchesInIndex(query) + " matches");
	//			System.out.println("Queried EnhancedSuffixArray - " + timer.getTimeString());
	//			System.out.println();
				//testIndex.printArray();
				//testIndex.serialize("/Users/froehler/indices/" + sequences[i] + ".esaidx");
				System.out.print("Serializing index");
				testIndex.serialize(file + ".ESAOptInt.esaidx");
				System.out.println(" - done in " + timer.getTimeString());
				System.out.println("Index construction done in " + timer.getTotalTimestring());
		}
		else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayIntOpt index = EnhancedSuffixArrayIntOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for string " + query);
			int[] result = {};
			for(int i=0; i<numQueries; i++){
				result = index.findMatchPositions(query);
			}
			System.out.println(" " + result.length + " matches - done in " + timer.getTimeString());
		}
		else{
			System.out.println("Usage: EnhancedSuffixArrayIntOpt query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: EnhancedSuffixArrayIntOpt create <serializedIndexFile>");
		}
	}

	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
		// TODO Auto-generated method stub
		
	}
}
