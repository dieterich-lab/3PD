/**
 * 
 */
package primerDesign.dsc.indexStructures.esa;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;

import primerDesign.dsc.FastStack;
import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SlimFastaParser;
import cern.colt.list.IntArrayList;
import cern.colt.list.ObjectArrayList;
import cern.colt.map.OpenIntIntHashMap;

/**
 * Implements the basic methods for an enhanced suffix array (ESA) as proposed by Kurtz et.al. 2004.
 * 
 * Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 * 
 * Implementations of this ESA have to provide an appropriate datastructure for storing the tables: suftab, lcptab, childtab,
 * the corresponding accessor methods for those tables and a proper initialization.
 * 
 * If sequence is to be searched case-insensitive, this has to be explicitely specified!
 *  
 * @author Sebastian Fršhler
 *
 */
public abstract class EnhancedSuffixArray implements DNASequenceIndex, Serializable{
	
	private static final long serialVersionUID = -4436486307664490034L;

	protected char[] sequence;
	protected int sequenceLength;
	protected String name;
	protected static final String TERMINATION_SYMBOL = "$";
	protected OpenIntIntHashMap bucketTab;
	protected int d;
	public static boolean printStatus = false;
	
	protected abstract boolean containsNextIndex(int i);
	protected abstract boolean containsDownIndex(int i);
	protected abstract boolean containsUpIndex(int i);
	
	protected abstract int getChildTabNext(int i);
	protected abstract int getChildTabDown(int i);
	protected abstract int getChildTabUP(int i);
	
	protected abstract void setChildTab(int i, int j);
	protected abstract int getChildTab(int i);	
	protected abstract int getChildTabLength();
	
	protected abstract void setLcpTab(int i, int j);
	protected abstract int getLcpTab(int i);	
	protected abstract int getLcpTabLength();
	
	protected abstract void setSufTab(int i, int j);
	protected abstract int getSufTab(int i);	
	protected abstract int getSufTabLength();
	
	public abstract void createIndex();

	/** 
	 * Initializes an enhanced suffix array, sets the max length of prefixes stored in the bucket table.
	 * 
	 * If an implementation of this class is initialized with a multi-sequence fasta file,
	 * only the first sequence is read to construct the index on!
	 * 
	 * @param file the file with the sequence(s) to create the esa from
	 * @throws IOException 
	 */
	protected EnhancedSuffixArray(File file) throws IOException {
		super();
		SlimFastaParser parser = new SlimFastaParser(file);
		SimpleContig contig = parser.parseNextContigIgnoreCase();
		this.sequence = contig.getSequence();
		this.sequenceLength = this.sequence.length;
		this.name = contig.getID();
		if(sequenceLength < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		this.bucketTab = new OpenIntIntHashMap();
		this.d = 8;
	}
	
	protected EnhancedSuffixArray(String sequence, String name){
		this.sequence = sequence.toCharArray();
		this.sequenceLength = this.sequence.length;
		this.name = name;
		if(sequenceLength < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		this.bucketTab = new OpenIntIntHashMap();
		this.d = 8;
	}
	
	protected EnhancedSuffixArray(char[] sequence, String name){
		this.sequence = sequence;
		this.sequenceLength = this.sequence.length;
		this.name = name;
		if(sequenceLength < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		this.bucketTab = new OpenIntIntHashMap();
		this.d = 8;
	}
	
	/**
	 * Initializes the lcp table - data structure dependent.
	 *
	 */
	protected abstract void initLcpTab();

	/**
	 * Computes the longest common prefix of a sorted suffix array for element 'i' w.r.t element 'i-1'.
	 *
	 */
	protected final void computeLCPTable() {
		initLcpTab();
		int lcp;
		int stop = this.getSufTabLength();
		for (int i = 0; i < stop; i++) {
			if(i == 0 || i == stop-1) setLcpTab(i, 0);
			else{
				lcp = computeLcp(i, i-1);
				setLcpTab(i, lcp);
			}
		}
	}
	
	/**
	 * Initializes the table 'childtab' - init method depends on datastructure used for table 'childtab'.
	 *
	 */
	protected abstract void initChildTable();

	/**
	 * Computes the child table.
	 * 
	 * Inits the child table, then calls: computeChildTabUpDown() and computeChildTabNext
	 *
	 */
	protected final void computeChildTable() {
		initChildTable();
		//computeCombinedChildTable();
		computeChildTabUpDown();
		computeChildTabNext();
	}

	/**
	 * Computes the 'next' value of the child table.
	 *
	 */
	private final void computeChildTabNext() {
		FastStack stack = new FastStack();
		stack.push(0);
		for (int i = 1; i <= sequenceLength; i++) {
			while(getLcpTab(i) < getLcpTab(stack.peek())){
				stack.pop();
			}
			assert(getLcpTab(i) >= getLcpTab(stack.peek()));
			if(getLcpTab(i) == getLcpTab(stack.peek())){
				int lastIndex = stack.pop();
				setChildTabNext(lastIndex, i);
			}
			stack.push(i);
		}
	}

	/**
	 * Computes the up and down values of the child table.
	 *
	 */
	private final void computeChildTabUpDown() {
			FastStack stack = new FastStack();
			stack.push(0);
			int lastIndex = -1;
			int top;
			for (int i = 1; i <= sequenceLength; i++) {
				while(getLcpTab(i) < getLcpTab(stack.peek())){
					lastIndex = stack.pop();
					top = stack.peek();
					if(getLcpTab(i) <= getLcpTab(top) && getLcpTab(top) != getLcpTab(lastIndex) && getChildTab(i) == Integer.MIN_VALUE){
						setChildTabDown(top, lastIndex);
						assert(getLcpTab(i) <= getLcpTab(stack.peek()) && getLcpTab(stack.peek()) < getLcpTab(lastIndex) && stack.peek() < lastIndex && lastIndex < i);
					}
				}
				top = stack.peek();
				assert(getLcpTab(i) >= getLcpTab(top));
				if(lastIndex != -1){ //  && getChildTab(i) == Integer.MIN_VALUE){
					//if(i<sequenceLength){
						setChildTabUp(i, lastIndex);
					//}
					//else if(i == sequenceLength){
					//	setChildTabUp(i, i-1);
					//}
					assert(getLcpTab(stack.peek()) <= getLcpTab(i) && getLcpTab(i) < getLcpTab(lastIndex) && stack.peek() < lastIndex && lastIndex < i);
	//				this.childtab_compara_up[i] = lastIndex;
					lastIndex = -1;
				}
				stack.push(i);
				//System.out.println("Stack: " + stack.toString());
			}
		}

	/**
	 * Computes the bucket table of prefix->first-index position associations to increase query performance.
	 * 
	 * The bucket table contains mappings of the form:
	 * pattern -> start index of interval in suffix array, where pattern is one prefix
	 * the hashCode of the pattern is stored in the bucket table instead of its string value!
	 *
	 */
	protected final void computeBucketTable() {
			//String currentString;
			int suftab;
			for(int i=sequenceLength-1; i>=0; i--){
				suftab = getSufTab(i);
				bucketTab.put(getHashCode(suftab, suftab + Math.min(d, sequenceLength-suftab)), i);
			}
		}

	/**
	 * Computes the hash code of sequence[i]...sequence[j-1] as String.hashCode().
	 * 
	 * @param i the start index
	 * @param j the stop index +1 s.th. stop-start=length
	 * 
	 * @return the hash code of sequence[i]...sequence[j-1] as String.hashCode()
	 */
	private final int getHashCode(int i, int j) {
		assert i>=0 && j>=i && i<=this.sequence.length && j<=this.sequence.length : "i: " + i + " j: " + j;
		int result = 0;
		int length = j-i;
		for(int k=0; k<length; k++){
			//result += sequence[i+k] * 31^(length-1-k);
			result = 31*result + this.sequence[i+k];
		}
		return result;
	}
	
	private final int getHashCode(char[] pattern) {
		assert pattern.length > 0;
		int result = 0;
		int length = pattern.length;
		for(int k=0; k<length; k++){
			//result += sequence[i+k] * 31^(length-1-k);
			result = 31*result + pattern[k];
		}
		return result;
	}

	/**
	 * Checks whether a pattern (its d-prefix) can be found in the bucket table.
	 * 
	 * @param pattern the pattern to search the d-prefix for
	 * 
	 * @return true iff the d-preix of the pattern can be found in the bucket table
	 */
	private final boolean isInBucketTab(String pattern) {
		assert pattern.length() > 0 && pattern != null;
		return bucketTab.containsKey(pattern.substring(0, Math.min(d,pattern.length())).hashCode());
	}

	/**
	 * Returns the start index of an interval where pattern is one prefix.
	 * 
	 * @param pattern the pattern
	 * 
	 * @return the start index of an interval where pattern is one prefix
	 */
	private final int getFromBucketTab(String pattern) {
		assert pattern.length() > 0 && pattern != null;
		if(isInBucketTab(pattern)) return bucketTab.get(pattern.substring(0, Math.min(d,pattern.length())).hashCode());
		else throw new IllegalStateException("Presence in bucket table is supposed to be checked before calling getFromBucketTab!");
	}

	/**
	 * Returns an interval where pattern is one prefix.
	 * 
	 * @param pattern the pattern
	 * @return an int[start,stop] interval where pattern is one prefix
	 */
	private final int[] getIntervalFromBucketTable(String pattern) {
		assert pattern.length() >= d && pattern != null;
		//pattern = pattern.substring(0, d);
		if(!isInBucketTab(pattern)) return new int[]{0,-1}; // pattern of length 'd' was not put into bucket table -> has NO associated interval in the interval tree!
		int start = getFromBucketTab(pattern);
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
	 * Convenience function including check whether this entry really contains a next index.
	 * 
	 * @param i the index
	 * @param j the value to set
	 */
	private void setChildTabUp(int i, int j) {
		if(!containsNextIndex(i)) setChildTab(i-1, j);
	}

	/**
	 * Convenience function including check whether this entry really contains a down index.
	 * @param i the index
	 * @param j the value to set
	 */
	private void setChildTabDown(int i, int j) {
		if(!containsNextIndex(i) && !containsUpIndex(i)) setChildTab(i, j);
	}

	/**
	 * Convenience function including check whether this entry really contains a next index.
	 * 
	 * @param i the index
	 * @param j the value to set
	 */
	private void setChildTabNext(int i, int j) {
		setChildTab(i, j);
	}

	/**
	 * Computes the longest common prefix of prefixes 'i' and 'j'.
	 * 
	 * @param i the index of the first prefix
	 * @param j the index of the second prefix
	 * @return the longest common prefix of strings stored in suffix table 'i' and 'j'
	 */
	private final int computeLcp(int i, int j) {
		int idxI = getSufTab(i);
		int idxJ = getSufTab(j);
		int end = sequenceLength - Math.max(idxI, idxJ);
		int lcp = 0;
		for(int k = 0; k < end; k++) {
			if(this.sequence[idxI + k] == this.sequence[idxJ + k]) lcp++;
			else break;
		}
		return lcp;
	}

	/**
	 * Returns a list of child intervals of the interval [i,j].
	 * 
	 * This list contains zero or more pairs of entries of the form: <i=start,i+1=end>
	 * 
	 * @param i the start position of the 'parent' interval
	 * @param j the stop position of the 'parent' interval
	 * 
	 * @return a list of child intervals of the 'parent' interval [i,j], format: [DOWN,UP]_x
	 */
	protected final int[] getChildIntervals(int i, int j) {
		assert i>=0 && j>=0 && i<this.getSufTabLength() && j<this.getSufTabLength();
		// init intervals list
		//int[] result = new int[0]; // list of child intervals: i%2=0 -> lower bound, i%2=1 -> upper bound
		IntArrayList result = new IntArrayList();
		
		if(i>=0 && j<sequenceLength){
			int i1;
			
			if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) i1 = getChildTabUP(j+1);
			else i1 = getChildTabDown(i);
			// add i, i1-1 to intervals list
			//result = addInterval(result, i, i1-1);
			result.add(i);
			result.add(i1-1);
			
			while(containsNextIndex(i1)){
				int i2 = getChildTabNext(i1);
				//if(i2 == j) break;
				// add i1, i2-1 to intervals list
				//result = addInterval(result, i1, i2-1);
				result.add(i1);
				result.add(i2-1);
				i1 = i2;
			}
			// add i1, j to intervals list
			//result = addInterval(result, i1, j);
			result.add(i1);
			result.add(j);
		}else{
			int k=i;
			while(containsNextIndex(k)){
				//result = addInterval(result, k, getChildTabNext(k)-1);
				result.add(k);
				result.add(getChildTabNext(k)-1);
				k = getChildTabNext(k);
			}
			if(j>=k){
				//result = addInterval(result, k, j);
				result.add(k);
				result.add(j);
			}
		}
		
		int[] temp = new int[result.size()];
		for(int z=0; z<result.size(); z++) temp[z] = result.get(z);
		
		return temp; //result;
	}

	/**
	 * Returns the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval.
	 * 
	 * @param i the start of the interval
	 * @param j the end of the interval
	 * @param p the character of the child interval to be found at position getLCP(i,j) at all suffices in the child interval
	 * 
	 * @return the child interval of the interval [i,j] with character p at position getLCP(i,j) at all suffices in that interval
	 */
	protected final int[] getChildInterval(int i, int j, char p) {
		assert i>=0 && j>=0 && i<this.getChildTabLength() && j<this.getChildTabLength() && p >= Character.MIN_VALUE && p<= Character.MAX_VALUE;
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
	 * 
	 * @return all match positions in the string the suffix array was constructed for
	 */
	public final int[] findMatchPositions(String patter) {
		if(patter.length() < 1) throw new IllegalArgumentException("Invalid pattern to scan with!");
		int c = 0;
		boolean queryFound = true;
		int elements;
		int[] interval = new int[0];
		// if pattern length >= d -> use bucket table to locate initial interval
		//if(patter.length() >= d){
		if(isInBucketTab(patter)){
			interval = getIntervalFromBucketTable(patter);
			c = d-1;
		}
		//}
		// else de-novo search
		else{ 
			interval = getChildInterval(0, sequenceLength, patter.charAt(c));
		}
		if(interval[0] <= interval[1]){ // an empty interval has length -1: e.g. index a=0, index b=-1
			int start = interval[0];
			int end = interval[1];
			int l;
			int min;
			while(interval[0] <= interval[1]  && c<patter.length() && queryFound){
				if(start != end){
					l = getLCP(start, end);
					min = Math.min(l,patter.length());
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
		if(interval[0] > interval[1]  || !queryFound) elements = 0;
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
	 * 
	 * @return true iff sequence.substring(sstart,sstart+length).equals(pattern.substring(pstart,pstart+length))
	 */
	private boolean isEqualToSequRegion(String pattern, int sstart, int pstart, int length) {
		assert(sstart >= 0 && sstart <= sequenceLength && pstart >= 0 && pstart <= pattern.length() && length >=0);
		if(sstart >= sequenceLength || pstart >= pattern.length() || length == 0 || sstart+length-1 >= sequenceLength || pstart+length-1 >= pattern.length()) return false;
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
	 * 
	 * @return the length of the longest common prefix of the interval [i,j]
	 */
	private final int getLCP(int i, int j) {
		assert i>= 0 && j>=0 && i<=sequenceLength && j<=sequenceLength;
		if(j == sequenceLength) return 0; // there is NO LCP between the whole string and the sentinel!
		else if(i < getChildTabUP(j+1) && getChildTabUP(j+1) <= j) return getLcpTab(getChildTabUP(j+1));
		else return getLcpTab(getChildTabDown(i));
	}

	/**
	 * Quicksort for an array of suffices.
	 * 
	 * @param from first index to be included in the sorting
	 * @param to last index to be included in the sorting
	 */
	protected void quicksortSuffices(int from, int to) {
		assert(from>=0 && to>=0 && from<this.sequenceLength && to<this.sequenceLength && from < to) : "Suffix index 'from' has to be smaller than 'to'!";
	    int i = from;
	    int j = to;
	    //int pivot = Math.round((from+to+0.0f)/2);
	    int pivot = (from+to)/2;
	    //int pivot = from + (to-from) * (int) Math.random();
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
	private void swap(int a, int b) {
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
	 * 
	 * @return -1|0|1 like compareTo for Strings does
	 */
	private int compareSuffix(int a, int b) {
		assert a>=0 && b>=0 && a<this.getSufTabLength() && b<this.getSufTabLength();
		int idxA = getSufTab(a);
		int idxB = getSufTab(b);
		int end = Math.min(sequenceLength-idxA, sequenceLength-idxB);
		for (int i = 0; i < end; i++) {
			int first = this.sequence[idxA+i];
			int second = this.sequence[idxB+i];
			if(first == second) continue;
			else if(first > second) return 1;
			else if(first < second) return -1;
		}
		if(sequenceLength-idxA > sequenceLength-idxB) return -1;
		else if(sequenceLength-idxA < sequenceLength-idxB) return 1;
		else return 0;
	}

	/**
	 * Returns the suffix starting at position 'i'.
	 * 
	 * @param i the index of the first character of the substring
	 * 
	 * @return the suffix starting at position 'i'
	 */
	protected String getSubstring(int i) {
		return new String(this.sequence, i, sequenceLength - i);
	}

	/**
	 * Returns the substring starting at position 'i' and ending at position 'j-1'.
	 * 
	 * @param i the start position
	 * @param j the stop position (the index of the first character NOT to be included)
	 * 
	 * @return the substring from position 'i' (inclusive) to position j (exclusive)
	 */
	private String getSubstring(int i, int j) {
		return new String(this.sequence, i, Math.min(sequenceLength - i, j));
	}

	/**
	 * Returns the sequence this index is constructed on.
	 * 
	 * @return the sequence this index is constructed on
	 */
	public String getSequence() {
		return new String(this.sequence);
	}
	
	public HashMap<Character, Integer> getStatistics(){
		HashMap<Character, Integer> result = new HashMap<Character, Integer>();
		
		int[] intervals = getChildIntervals(0, sequenceLength);
		// intervals have format: [DOWN,UP]_x, latz child interval of root is termination charater which does NOT appear in the original sequence!
		for(int i=0; i<intervals.length-3; i+=2){
			result.put(sequence[getSufTab(intervals[i])], intervals[i+1] - intervals[i] + 1);
			//System.err.println(sequence[getSufTab(intervals[i])] + "Interval: start: " + intervals[i] + " end: " + intervals[i+1]);
		}		
		return result;
	}
	
	/**
	 * Returns the contig this index was created on.
	 * 
	 * @return the contig this index was created on
	 */
	public SimpleContigImpl[] getContig(){
		return new SimpleContigImpl[]{new SimpleContigImpl(this.name, this.sequence)};
	}
	
	/**
	 * Returns the length of the sequence this index is constucted on.
	 * 
	 * @return the length of the sequence this index is constucted on
	 */
	public int getSequenceLength(){
		return this.sequenceLength;
	}

	/**
	 * Returns the positions of matches of string 'searchString' and its reverse complement! in the index.
	 * 
	 * @param searchString the string to query the index with
	 * 
	 * @return a list containing the position of all matches of 'searchString' in the index
	 */
	public ObjectArrayList findHitPositions(String searchString) {
			if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
			int[] positionsFW = findMatchPositions(searchString.toUpperCase());
			int[] positionsRev = findMatchPositions(SeqTools.revcompDNA(searchString.toUpperCase().toCharArray()));
			
			ObjectArrayList result = new ObjectArrayList();
			for(int position : positionsFW){
				result.add(new IndexHitImpl(new SimpleContigImpl(this.name, this.sequence), position, true));
			}
			for(int position : positionsRev){
				result.add(new IndexHitImpl(new SimpleContigImpl(this.name, this.sequence), position, false));
			}
			return result;
			
	//		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
	//		int[] positionsFW = findMatchPositions(searchString.toUpperCase());
	//		Integer[] result = new Integer[positionsFW.length];
	//		//System.arraycopy(positions, 0, result, 0, positions.length);
	//		for(int i=0; i<positionsFW.length; i++) result[i] = positionsFW[i];
	//		return result;
	}
	
	/**
	 * Returns the positions of matches of string 'searchString' on the forward strand in the index.
	 * 
	 * @param searchString the string to query the index with
	 * 
	 * @return a list containing the position of all matches of 'searchString' on the forward strand in the index
	 */
	public ObjectArrayList findForwardHits(String searchString){
		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
		int[] positionsFW = findMatchPositions(searchString.toUpperCase());
		
		ObjectArrayList result = new ObjectArrayList();
		for(int position : positionsFW){
			result.add(new IndexHitImpl(new SimpleContigImpl(this.name, this.sequence), position, true));
		}
		
		return result;
	}
	
	/**
	 * Returns the positions of matches of string 'searchString' on the reverse strand in the index.
	 * 
	 * @param searchString the string to query the index with
	 * 
	 * @return a list containing the position of all matches of 'searchString' on the reverse strand in the index
	 */
	public ObjectArrayList findReverseHits(String searchString){
		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
		int[] positionsRev = findMatchPositions(SeqTools.revcompDNA(searchString.toUpperCase().toCharArray()));
		
		ObjectArrayList result = new ObjectArrayList();
		for(int position : positionsRev){
			result.add(new IndexHitImpl(new SimpleContigImpl(this.name, this.sequence), position, false));
		}
		return result;
	}

	/**
	 * Returns the number of matches of string 'searchString' and its reverse complement! in the index.
	 * 
	 * @param searchString the string to query the index with
	 * 
	 * @return the number of matches of string 'searchString' in the index
	 */
	public int findHitCount(String searchString) {
		if(searchString == null || searchString.length() < 1) throw new IllegalArgumentException("Invalid search string!");
		return findMatchPositions(searchString.toUpperCase()).length + findMatchPositions(SeqTools.revcompDNA(searchString.toUpperCase().toCharArray())).length;
	}

	/**
	 * Serializes the current index to file 'filename'.
	 * 
	 * @param filename the filename to serialize the index to
	 */
	public final void serialize(File filename) {
		try{
			//ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Returns the name of the sequence of this ESA.
	 * 
	 * @return the name of the sequence of this ESA
	 */
	public final String getName(){
		return this.name;
	}
	
	/**
	 * Prints some information about this ESA to STDOUT.
	 */
	public void printESA(){
		System.out.println("i\tsuftab\tlcptab\tchildtab\tsequence");
		for(int i=0; i<sequenceLength; i++){
			System.out.println(i + "\t" + getSufTab(i) + "\t" + getLcpTab(i) + "\t" + getChildTab(i) + "\t" + getSubstring(getSufTab(i)) + EnhancedSuffixArray.TERMINATION_SYMBOL);
		}
		System.out.println(sequenceLength + "\t" + getSufTab(sequenceLength) + "\t" + getLcpTab(sequenceLength) + "\t" + getChildTab(sequenceLength) + "\t" + EnhancedSuffixArray.TERMINATION_SYMBOL);
	}
}