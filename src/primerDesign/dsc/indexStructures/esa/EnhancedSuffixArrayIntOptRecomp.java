/**
 * 
 */
package primerDesign.dsc.indexStructures.esa;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.text.NumberFormat;
import java.util.Arrays;

import primerDesign.util.SimpleTimer;
import cern.colt.map.OpenIntIntHashMap;

/**
 * This class implements the enhanced suffix array as proposed by Kurtz et.al. 2004 using byte tables (memory usage: ~14n).
 * 
 * Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 * 
 * For increased memory performance, the tables 'lcptab' and 'childtab' store byte values and therefore have to recompute some values!
 * 
 * @author Sebastian Fršhler
 *
 */
public class EnhancedSuffixArrayIntOptRecomp extends EnhancedSuffixArray{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int[] suftab;
	private byte[] lcptab;
	private byte[] childtab;
	private OpenIntIntHashMap lcpmap;  // a table to lookup clipping lcp values

	/**
	 * Initializes an enhanced suffix array.
	 * 
	 * @param sequence the sequence to construct the ESA from
	 * @param name the name of the sequence
	 */
	public EnhancedSuffixArrayIntOptRecomp(String sequence, String name){
		super(sequence, name);
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Initializes an enhanced suffix array.
	 * 
	 * @param file the FASTA file to construct the ESA from
	 */
	public EnhancedSuffixArrayIntOptRecomp(File file) throws IOException{
		super(file);
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#createIndex(int, boolean)
	 */
	@Override
	public void createIndex() {
		this.suftab = new int[sequenceLength+1];
		this.lcptab = new byte[sequenceLength+1];
		//this.lcpmap = new IntArrayList();
		this.lcpmap = new OpenIntIntHashMap();
		this.childtab = new byte[sequenceLength+1];
		
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
	 * Recomuptes the up value at position 'i'.
	 * 
	 * @param i the index to recompute the (clipped) value for
	 * 
	 * @return the up value at position 'i'
	 */
	private int recompChildTabUp(int i){
		int candidate = -1;
		int min = Integer.MAX_VALUE;
		int reference = getLcpTab(i);
		for(int j=i-1; j>=0; j--){
			min = Math.min(min, getLcpTab(j));
			if(getLcpTab(j) > reference && min >= getLcpTab(j)) candidate = j;
			else if(getLcpTab(j) <= reference) break;
		}
		return candidate;
	}
	
	/**
	 * Recomputes the down value at position 'i'.
	 * 
	 * @param i the index to recompute the (clipped) value for
	 * 
	 * @return the down value at position 'i'
	 */
	private int recompChildTabDown(int i){
		int candidate = -1;
		int min = Integer.MAX_VALUE;
		int reference = getLcpTab(i);
		for(int j=i+1; j<=sequenceLength; j++){
			if(getLcpTab(j) > reference && min > getLcpTab(j)) candidate = j;
			else if(getLcpTab(j) <= reference) break;
			min = Math.min(min, getLcpTab(j));
		}
		return candidate;
	}
	
	/**
	 * Recomputes the next value at position 'i' .
	 * 
	 * @param i the index to recompute the (clipped) value for
	 * 
	 * @return the next value at position 'i'
	 */
	private int recompChildTabNext(int i){
		int candidate = -1;
		int min = Integer.MAX_VALUE;
		int reference = getLcpTab(i);
		for(int j=i+1; j<=sequenceLength; j++){
			if(getLcpTab(j) == reference && min > reference) candidate = j;
			else if(getLcpTab(j) < reference) break;
			min = Math.min(min, getLcpTab(j));
		}
		return candidate;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsDownIndex(int)
	 */
	@Override
	protected boolean containsDownIndex(int i) {
		assert(i>=0 && i<this.childtab.length);
		return getChildTab(i) < Byte.MAX_VALUE && getChildTab(i) != -1 ? getLcpTab(getChildTab(i)) > getLcpTab(i) : recompChildTabDown(i) != -1;
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsNextIndex(int)
	 */
	@Override
	protected boolean containsNextIndex(int i) {
		assert(i>=0 && i<this.childtab.length);
		return i != sequenceLength - 1 && getChildTab(i) < Byte.MAX_VALUE && getChildTab(i) != -1 ? getLcpTab(getChildTab(i)) == getLcpTab(i) && getChildTab(i) > i : recompChildTabNext(i) != -1 && recompChildTabNext(i) > i; // if no child value can be found during ChildTabNext-recomputation, value '-1' is returned to signal this
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsUpIndex(int)
	 */
	@Override
	protected boolean containsUpIndex(int i) {
		assert(i>=0 && i<this.childtab.length);
		if(i == sequenceLength) return true;
		else return getLcpTab(i) > getLcpTab(i+1);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTab(int)
	 */
	@Override
	protected int getChildTab(int i) {
		assert(i>=0 && i<this.childtab.length);
		return (this.childtab[i] < Byte.MAX_VALUE) ? this.childtab[i] - Byte.MIN_VALUE + i : -1; 
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabDown(int)
	 */
	@Override
	protected int getChildTabDown(int i) {
		assert containsDownIndex(i) : "Child table does not contain a down value for index " + i + "!";
		return getChildTab(i) < Byte.MAX_VALUE && getChildTab(i) != -1 ? getChildTab(i) : recompChildTabDown(i);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabUP(int)
	 */
	@Override
	protected int getChildTabUP(int i) {
		assert containsUpIndex(i-1) : "Child table does not contain an up value for index " + i + "!";
		return getChildTab(i-1) < Byte.MAX_VALUE && getChildTab(i-1) != -1 ? getChildTab(i-1) : recompChildTabUp(i);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabNext(int)
	 */
	@Override
	protected int getChildTabNext(int i) {
		assert containsNextIndex(i) : "Child table does not contain a next value for index " + i + "!";
		return getChildTab(i) < Byte.MAX_VALUE && getChildTab(i) != -1 ? getChildTab(i) : recompChildTabNext(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabLength()
	 */
	@Override
	protected int getChildTabLength() {
		return this.childtab.length;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getLcpTabLength()
	 */
	@Override
	protected int getLcpTabLength() {
		return this.lcptab.length;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getSufTabLength()
	 */
	@Override
	protected int getSufTabLength() {
		return this.suftab.length;
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getLcpTab(int)
	 */
	@Override
	protected int getLcpTab(int i) {
		assert(i>=0 && i<this.lcptab.length);
		//return (this.lcptab[i] < Byte.MAX_VALUE) ? this.lcptab[i] - Byte.MIN_VALUE : getMap(this.lcpmap, i, 0, this.lcpmap.size()-1);
		return (this.lcptab[i] < Byte.MAX_VALUE) ? this.lcptab[i] - Byte.MIN_VALUE : this.lcpmap.get(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getSufTab(int)
	 */
	@Override
	protected int getSufTab(int i) {
		assert(i>=0 && i<this.suftab.length);
		return this.suftab[i];
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#initChildTable()
	 */
	@Override
	protected void initChildTable() {
		Arrays.fill(this.childtab, Byte.MIN_VALUE);
	}
	
	protected void initLcpTab(){
		Arrays.fill(this.lcptab, Byte.MIN_VALUE);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setChildTab(int, int)
	 */
	@Override
	protected void setChildTab(int i, int j) {
		assert i>=0 && i<this.suftab.length && j>=0 && j<this.suftab.length : "i and j must be in interval [0,sequenceLength]";
		if(j-i + Byte.MIN_VALUE < Byte.MAX_VALUE && j-i >= 0) this.childtab[i] = (byte) (j-i+Byte.MIN_VALUE);
		else{
			this.childtab[i] = Byte.MAX_VALUE;
		}
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setLcpTab(int, int)
	 */
	@Override
	protected void setLcpTab(int i, int j) {
		assert(i>=0 && j>=0 && i<this.lcptab.length && j<=this.sequenceLength);
		if(j + Byte.MIN_VALUE < Byte.MAX_VALUE) this.lcptab[i] = (byte) (j + Byte.MIN_VALUE);
		else{
			this.lcptab[i] = Byte.MAX_VALUE;
			//putLcpMap(i, j);  // since computeLCPTable calls setLcpTable only one per index i, no 'contains' lookup needs to be done in lcpmap!!
			this.lcpmap.put(i, j);
		}
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setSufTab(int, int)
	 */
	@Override
	protected void setSufTab(int i, int j) {
		assert(i>=0 && j>=0 && i<this.suftab.length && j<this.suftab.length);
		this.suftab[i] = j;
	}
	
	public static EnhancedSuffixArrayIntOptRecomp deserialize(String filename){
		EnhancedSuffixArrayIntOptRecomp result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (EnhancedSuffixArrayIntOptRecomp) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IOException{		
		EnhancedSuffixArrayIntOptRecomp.printStatus = false;
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			EnhancedSuffixArrayIntOptRecomp testIndex;
				System.out.println("Reading sequence " + file);
				System.out.println("Reading sequence " + file);
				{					
					System.out.println("Creating EnhancedSuffixArray for first sequence in file " + file);
					
					System.out.println("Constructing EnhancedSuffixArray");
					//sequence = "ACAAACATAT";
					testIndex = new EnhancedSuffixArrayIntOptRecomp(new File(file));
				}
				//System.gc();
				testIndex.createIndex();
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
				testIndex.serialize(new File(file + ".ESAOptIntRecomp.esaidx"));
				System.out.println(" - done in " + timer.getTimeString());
				System.out.println("Index construction done in " + timer.getTotalTimestring());
		}
		else if(mode.equals("query")){
			String query = args[2].toUpperCase();
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayIntOptRecomp index = EnhancedSuffixArrayIntOptRecomp.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for string " + query);
			int result = 0;
			for(int i=0; i<numQueries; i++){
				result = index.findHitCount(query);
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
		}
		else{
			System.out.println("Usage: EnhancedSuffixArrayIntOptRecomp query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: EnhancedSuffixArrayIntOptRecomp create <serializedIndexFile>");
		}
	}
}
