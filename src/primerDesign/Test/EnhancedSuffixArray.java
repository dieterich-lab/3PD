package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Stack;
import java.util.regex.Pattern;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;


/**
 * This class implements the enhanced suffix array as proposed by Kurtz et.al. 2004
 * 
 * Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 * 
 * @author froehler
 *
 */
public class EnhancedSuffixArray implements DNASequenceIndex{
	private String sequence;
	private boolean includesScanRegion;
	private int maxWordSize;
	private int sequenceLength;
	private static String TERMINATION_SYMBOL = "$";
	private int[] suftab;
	private int[] lcptab;
	private int[] childtab;
	
	public EnhancedSuffixArray(){
//		Runtime runtime = Runtime.getRuntime();
//		NumberFormat format = NumberFormat.getInstance();
//		System.gc();
//		System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
	}
	
	/**
	 * Creates an enhanced suffix array for the string 'sequence'.
	 * @param sequence the string to construct the enhanced suffi array for
	 */
	public void createIndex(String sequence){
		this.sequence = sequence;
		this.sequenceLength = sequence.length();
		this.suftab = new int[sequenceLength+1];
		this.lcptab = new int[sequenceLength+1];
		this.childtab = new int[sequenceLength+1];
		
		for(int i=0; i< sequenceLength; i++){
			this.suftab[i] = i;
		}
		this.suftab[sequenceLength] = Integer.MAX_VALUE;
		// sort table
//		Runtime runtime = Runtime.getRuntime();
//		NumberFormat format = NumberFormat.getInstance();
//		System.gc();
//		System.out.println("Before sort: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		quicksortSuffices(0, sequenceLength-1);
//		System.gc();
//		System.out.println("Before LCP Table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		computeLCPTable();
//		System.gc();
//		System.out.println("Before child table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		computeChildTable();
//		System.gc();
//		System.out.println("After all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
	}
	
	/**
	 * Computes the longest common prefix of a sorted suffix array for element 'i' w.r.t element 'i-1'.
	 *
	 */
	private void computeLCPTable(){
		for(int i=0; i<this.lcptab.length; i++) this.lcptab[i] = 0;
		for (int i = 0; i < this.suftab.length; i++) {
			if(i == 0 || i == this.suftab.length-1) setLcpTab(i, 0);
			else{
				for (int j = 0; j < sequenceLength - Math.max(getSufTab(i), getSufTab(i-1)); j++) {
					if(this.sequence.charAt(getSufTab(i) + j) == this.sequence.charAt(getSufTab(i-1) + j)) setLcpTab(i, getLcpTab(i)+1);
					else break;
				}
			}
		}
	}
	
	private void computeChildTable(){
		computeChildNextIndex();
		computeChildUpDown();
	}
	
	private void computeChildUpDown(){
		int lastIndex = -1;
		Stack<Integer> stack = new Stack<Integer>();
		stack.push(0);
		for (int i = 1; i <= sequenceLength; i++) {
			while(getLcpTab(i) < getLcpTab(stack.peek())){
				lastIndex = stack.pop();
				if(getLcpTab(i) <= getLcpTab(stack.peek()) && getLcpTab(stack.peek()) != getLcpTab(lastIndex) && !containsNextIndex(stack.peek())) setChildTab(stack.peek(), lastIndex);
			}
			if(lastIndex != -1 && !containsNextIndex(i-1) && !containsDownIndex(i-1)){
				setChildTab(i-1, lastIndex);
				lastIndex = -1;
			}
			stack.push(i);
		}
	}
	
	private void computeChildNextIndex(){
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
	
	/**
	 * Convenience function to return the entry 'SUF_TAB' for index i.
	 * 
	 * @param i the index
	 * @return the entry 'SUF_TAB' for index i
	 */
	protected int getSufTab(int i){
		assert(i>=0);
		return this.suftab[i];
	}
	
	private void setSufTab(int i, int j){
		assert(i>=0 && j>=0);
		this.suftab[i] = j;
	}
	
	/**
	 * Convenience function to return the entry 'LCP_TABLE' for index 'i'.
	 * 
	 * @param i the index
	 * @return the entry 'LCP_TABLE' for index 'i'
	 */
	protected int getLcpTab(int i){
		assert(i>=0);
		return this.lcptab[i];
	}
	
	protected void setLcpTab(int i, int j){
		assert(i>=0 && j>=0);
			this.lcptab[i] = j;
	}
	
	/**
	 * Convenience function to return the entry 'CHILD_TABLE' for index 'i'.
	 * 
	 * @param i the index
	 * @return the entry 'CHILD_TABLE' for index 'i'
	 */
	protected int getChildTab(int i){
		assert(i>=0);
		return this.childtab[i];
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is an up index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is an up index
	 */
	private boolean containsUpIndex(int i){
		assert(i>=0);
		if(i == sequenceLength) return true;
		else return getLcpTab(i) > getLcpTab(i+1);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is a down index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is a down index
	 */
	private boolean containsDownIndex(int i){
		assert(i>=0);
		return getLcpTab(getChildTab(i)) > getLcpTab(i);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is the next index.
	 * 
	 * @param i the child to test
	 * @return true iff the child entry for suffix 'i' is the next index
	 */
	private boolean containsNextIndex(int i){
		assert(i>=0);
		return getLcpTab(getChildTab(i)) == getLcpTab(i) && getChildTab(i) > i;
	}
	
	/**
	 * Returns a list of child intervals of the interval [i,j]
	 * 
	 * @param i the start position of the 'parent' interval
	 * @param j the stop position of the 'parent' interval
	 * @return a list of child intervals of the 'parent' interval [i,j], format: [DOWN,UP]_x
	 */
	protected ObjectArrayList getChildIntervals(int i, int j){
		// init intervals list
		ObjectArrayList result = new ObjectArrayList(); 
		
		if(i>=0 && j<sequenceLength){
			int i1;
			
			if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) i1 = getChildTabUP(j+1);
			else i1 = getChildTabDown(i);
			// add i, i1-1 to intervals list
			result.add(new Interval(i, i1-1));
			
			while(containsNextIndex(i1)){
				int i2 = getChildTabNext(i1);
				// add i1, i2-1 to intervals list
				result.add(new Interval(i1,i2-1));
				i1 = i2;
			}
			// add i1, j to intervals list
			result.add(new Interval(i1, j));
		}else{
			int k=i;
			while(containsNextIndex(k)){
				result.add(new Interval(k, getChildTabNext(k)-1));
				k = getChildTabNext(k);
			}
			if(j>=k) result.add(new Interval(k,j));
		}
		
		return result;
	}
	
	/**
	 * For debugging purpose...
	 * @param intervals
	 * @return
	 */
	protected ObjectArrayList getChildIntervalsValues(ObjectArrayList intervals){
		ObjectArrayList values = new ObjectArrayList();
		for(int i=0; i<intervals.size(); i++){
			values.add(new Integer[]{((Interval)intervals.get(i)).getLowerBoundary(),((Interval)intervals.get(i)).getUpperBoundary()});
		}
		return values;
	}
	
	protected int[] intervalToArray(Interval interval){
		return new int[]{interval.getLowerBoundary(), interval.getUpperBoundary()};
	}
	
	/**
	 * Returns the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval.
	 * 
	 * @param i the start of the interval
	 * @param j the end of the interval
	 * @param p the character of the child interval to be found at position getLCP(i,j) at all suffices in the child interval
	 * @return the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval
	 */
	protected Interval getChildInterval(int i, int j, char p){
		Interval result = new Interval();
		int lcpOffset = getLCP(i, j);
		ObjectArrayList intervals = getChildIntervals(i, j);
		for(int k=0; k<intervals.size(); k++){
			int lowerBound = ((Interval)intervals.get(k)).getLowerBoundary(); 
			int index = getSufTab(lowerBound) + lcpOffset;
			if(index < sequenceLength && this.sequence.charAt(index) == p){
				result = (Interval) intervals.get(k);
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
		Interval interval = getChildInterval(0, sequenceLength, patter.charAt(c));
		if(!interval.isEmptyInterval){
			int start = interval.getLowerBoundary();
			int end = interval.getUpperBoundary();
			while(!interval.isEmptyInterval  && c<patter.length() && queryFound){
				if(start != end){
					int l = getLCP(start, end);
					int min = Math.min(l,patter.length());
					queryFound = (this.sequence.substring(getSufTab(start)+c, getSufTab(start)+min).equals(patter.substring(c, min)));
					c = min;
					if(c<patter.length()) interval = getChildInterval(start, end, patter.charAt(c));
					else break;
					if(interval.isEmptyInterval) break;
					else{
						start = interval.getLowerBoundary();
						end = interval.getUpperBoundary();
					}
				}
				else{
					queryFound = (this.sequence.substring(getSufTab(start)+c, getSufTab(start)+patter.length()).equals(patter.substring(c, patter.length())));
					break;
				}
			}
		}
		if(interval.isEmptyInterval) elements = 0;
		else elements = interval.getUpperBoundary()-interval.getLowerBoundary()+1;
		int[] positions = new int[elements];
		for(int i=0; i<elements; i++){
			positions[i] = getSufTab(interval.getLowerBoundary() + i);
		}
		return positions;
	}
	
	/**
	 * Returns the length of the longest common prefix of the interval [i,j].
	 * 
	 * @param i the start of the interval
	 * @param j the stop of the interval
	 * @return the length of the longest common prefix of the interval [i,j]
	 */
	private int getLCP(int i, int j){
		if(j == sequenceLength) return 0; // there is NO LCP between the whole string and the sentinel!
		else if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) return getLcpTab(getChildTabUP(j+1));
		else return getLcpTab(getChildTabDown(i));
	}
	
	private int getChildTabUP(int i){
		if(containsUpIndex(i-1)) return getChildTab(i-1);
		else throw new IllegalArgumentException("Child table does not contain an up value for suffix " + i + "!");
	}
	
	private int getChildTabDown(int i){
		if(containsDownIndex(i)) return getChildTab(i);
		else throw new IllegalArgumentException("Child table does not contain a down value for suffix " + i + "!");
	}
	
	private int getChildTabNext(int i){
		if(containsNextIndex(i)) return getChildTab(i);
		else throw new IllegalArgumentException("Child table does not contain a next value for suffix " + i + "!");
	}
	
	private void setChildTab(int i, int j){
		if(i<0 || i>=this.suftab.length || j<0 || j>=this.suftab.length) throw new IllegalArgumentException("i and j must be in interval [0,sequenceLength]");
		this.childtab[i] = j;
	}
	
	private String getStringLabel(int element){
		return this.sequence.substring(element);
	}
	
	/**
	 * Quicksort for an array of suffices.
	 * 
	 * @param from first index to be included in the sorting
	 * @param to last index to be included in the sorting
	 */
	private void quicksortSuffices(int from, int to){
		if(from >= to) throw new IllegalArgumentException("Suffix index 'from' has to be smaller than 'to'!");
	    int i = from;
	    int j = to;
	    //int pivot = Math.round((from+to+0.0f)/2);
	    int pivot = from;
	    
	    //  partition
	    do{    
	        while (i < pivot && compareSuffix(i, pivot) < 0) i++; 
	        while (j > pivot && compareSuffix(j, pivot) > 0) j--;
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
		return this.sequence.substring(getSufTab(a)).compareTo(this.sequence.substring(getSufTab(b)));
//		for (int i = 0; i < Math.min(sequenceLength-a, sequenceLength-b); i++) {
//			if(this.sequence.charAt(a+i) > this.sequence.charAt(b+i)) return 1;
//			else if(this.sequence.charAt(a+i) < this.sequence.charAt(b+i)) return -1;
//		}
//		if(sequenceLength-a > sequenceLength-b) return 1;
//		else if(sequenceLength-a < sequenceLength-b) return -1;
//		else return 0;
	}
	
	/**
	 * Returns a textual representation of the current suffix array.
	 *
	 */
	public void printArray(){
		for(int i=0; i<this.suftab.length; i++){
			if(getSufTab(i) == Integer.MAX_VALUE) System.out.print(TERMINATION_SYMBOL);
			else System.out.print(this.sequence.substring(getSufTab(i)) + TERMINATION_SYMBOL);
			System.out.println("\t" + getLcpTab(i) + "\t" + getChildTab(i));
		}
		System.out.println();
	}
	
	/**
	 * Encapsulates an interval in the enhanced suffix array.
	 * 
	 * @author froehler
	 *
	 */
	
	private class Interval{
		private int lowerBoundary = 0;
		private int upperBoundary = Integer.MAX_VALUE;
		private boolean isEmptyInterval;
		
		public Interval(){
			this.isEmptyInterval = true;
		}
		
		public Interval(int i, int j){
			setLowerBoundary(i);
			setUpperBoundary(j);
			this.isEmptyInterval = false;
		}
		
		public int getLowerBoundary(){
			if(!isEmptyInterval) return this.lowerBoundary;
			else throw new IllegalStateException("No lower bound has been set for this interval!");
		}
		
		public void setLowerBoundary(int i){
			if(i >= 0 && i <= sequenceLength && i <= this.upperBoundary) this.lowerBoundary = i;
			else throw new IllegalArgumentException("Interval boundaries must be in interval [0," + (sequenceLength-1) + "]");
		}
		
		public int getUpperBoundary(){
			if(!isEmptyInterval) return this.upperBoundary;
			else throw new IllegalStateException("No upper bound has been defined for this interval!");
		}
		
		public void setUpperBoundary(int i){
			if(i >= 0 && i <= sequenceLength && i >= this.lowerBoundary) this.upperBoundary = i;
			else throw new IllegalArgumentException("Interval boundaries must be in interval [0," + (sequenceLength-1) + "]");
		}
		
		public boolean isEmptyInterval(){
			return this.isEmptyInterval;
		}
	}
	
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
		createIndex(sequence.toUpperCase());
		this.maxWordSize = sequenceLength;
		this.includesScanRegion = includesScanRegion;
	}

	public int getMaxWordSize() {
		return this.maxWordSize;
	}

	public String getSequence() {
		return this.sequence;
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	public Integer[] searchMatchPositionsInIndex(String searchString) {
		int[] positions = findMatchPositions(searchString.toUpperCase());
		Integer[] result = new Integer[positions.length];
		System.arraycopy(positions, 0, result, 0, positions.length);
		return result;
	}

	public int searchNbMatchesInIndex(String searchString) {
		return findMatchPositions(searchString.toUpperCase()).length;
	}
	
	public static void main(String[] args) throws IOException{		
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		//String[] sequences = {"100kb.dna", "200kb.dna", "400kb.dna", "800kb.dna", "../src/yeast-chrXII.dna", "Cel/Celegans.fa", "Ppa/Ppacificus.fa"};
		String[] sequences = {"../src/yeast-chrXII.dna"};
		String path = "/export/Sebastian/PrimerDesign/Testsequenzen/";
		DNASequenceIndex testIndex;
		BufferedReader reader;
		String sequence;
		String line;
		Pattern fastaStart = Pattern.compile(">.*");
		for(int i=0; i<sequences.length; i++){
			System.out.println("Reading sequence " + sequences[i]);
			reader = new BufferedReader(new FileReader(path + sequences[i]));
			{
				StringBuffer seq = new StringBuffer();
				while((line = reader.readLine()) != null){
					if(!fastaStart.matcher(line).matches()){
						seq.append(line.trim());
					}
				}
				sequence = seq.toString();
			}
			System.out.println("Reading sequence " + sequences[i] + " (length: " + sequence.length() + ") took " + timer.getTimeString());
			
			System.out.println("Creating EnhancedSuffixArray for sequence " + sequences[i]);
			
			System.out.println("Constructing EnhancedSuffixArray");
			testIndex = new EnhancedSuffixArray();
			testIndex.createIndex(sequence, Integer.MAX_VALUE, true);
			System.out.println("EnhancedSuffixArray construction complete - " + timer.getTimeString());
				String query = "ATG";
				System.out.println("Querying EnhancedSuffixArray for: " + query);
				System.out.println(testIndex.searchNbMatchesInIndex(query) + " matches");
				System.out.println("Queried EnhancedSuffixArray - " + timer.getTimeString());
				System.out.println();
		}
	}
}
