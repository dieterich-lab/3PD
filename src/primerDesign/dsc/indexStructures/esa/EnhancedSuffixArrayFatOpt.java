package primerDesign.dsc.indexStructures.esa;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;

import primerDesign.algo.KoAluruSuffixSort;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;


/**
 * This class implements the enhanced suffix array as proposed by Kurtz et.al. 2004 using integer tables (memory usage: ~20n).
 * 
 * Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 * 
 * For increased query performance, the tables 'lcptab' and 'childtab' store integer values and therefore do not have to
 * recompute any values!
 * 
 * @author Sebastian Fršhler
 *
 */
public class EnhancedSuffixArrayFatOpt extends EnhancedSuffixArray{

	private static final long serialVersionUID = -4436486307664490033L;

	private boolean includesScanRegion;
	private int maxWordSize;
	int[] suftab;
	int[] lcptab;
	int[] childtab;
	
	/**
	 * Initializes an enhanced suffix array of sequence 'sequence'.
	 * 
	 * @param file the file containing the sequence(s) to create the enhance suffix array on
	 * @throws IOException 
	 */
	public EnhancedSuffixArrayFatOpt(File file) throws IOException{
		super(file);
		//if(sequence.length() < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Initializes an enhanced suffix array from a sequence and a sequence name.
	 * 
	 * @param sequence the sequence to construct the ESA from
	 * @param name the name of the sequence
	 */
	public EnhancedSuffixArrayFatOpt(String sequence, String name){
		super(sequence, name);
		//if(sequence.length() < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Initializes an enhanced suffix array from a sequence and a sequence name.
	 * 
	 * @param sequence the sequence to construct the ESA from
	 * @param name the name of the sequence
	 */
	public EnhancedSuffixArrayFatOpt(char[] sequence, String name){
		super(sequence, name);
		//if(sequence.length() < 1 || sequence == null) throw new IllegalArgumentException("Illegal sequence for index creation!");
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/**
	 * Creates an enhanced suffix array. 
	 * 
	 * Parameters have to be provided to conform to the 'DNASequenceIndex' interface which also covers tree based indices
	 * where those parameters make 'much more sense'...
	 * 
	 */
	@Override
	public void createIndex(){

		this.suftab = new int[sequenceLength+1];
		this.lcptab = new int[sequenceLength+1];
		this.childtab = new int[sequenceLength+1];
		
		for(int i=0; i< sequenceLength; i++){
			setSufTab(i, i);
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
		//quicksortSuffices(0, sequenceLength-1);
		{
			KoAluruSuffixSort sorter = new KoAluruSuffixSort();
			this.suftab = sorter.getSuffixArray(this.sequence);
		}
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
	 * Initializes the child table.
	 */
	protected void initChildTable(){
		Arrays.fill(this.childtab, Integer.MIN_VALUE);
	}
	
	/**
	 * Initializes the longest common prefix table.
	 */
	protected void initLcpTab(){
		Arrays.fill(this.lcptab, Integer.MIN_VALUE);
	}
	
	/**
	 * Convenience function to return the entry 'SUF_TAB' for index i.
	 * 
	 * @param i the index
	 * 
	 * @return the entry 'SUF_TAB' for index i
	 */
	@Override
	protected int getSufTab(int i){
		assert(i>=0 && i<this.suftab.length);
		return this.suftab[i];
	}
	
	/**
	 * Convenience function to set the entry 'SUF_TAB' for index i and value j.
	 * 
	 * @param i the index
	 * @param j the value of 'SUF_TAB(i)' to set
	 */
	@Override
	protected void setSufTab(int i, int j){
		assert(i>=0 && j>=0 && i<this.getSufTabLength() && j<this.getSufTabLength());
		this.suftab[i] = j;
	}
	
	/**
	 * Convenience function to return the entry 'LCP_TABLE(i)=j'.
	 * 
	 * @param i the index
	 * 
	 * @return the entry 'LCP_TABLE' for index 'i'
	 */
	@Override
	protected int getLcpTab(int i){
		assert(i>=0 && i<this.lcptab.length);
		return this.lcptab[i];
	}
	
	/**
	 * Convenience function to set the entry 'LCP_TAB(i)=j'.
	 * 
	 * @param i the index
	 * @param j the value of 'LCP_TAB(i)' to set
	 */
	@Override
	protected void setLcpTab(int i, int j){
		assert(i>=0 && j>=0 && i<this.lcptab.length && j<=this.sequenceLength);
		this.lcptab[i] = j;
	}
	
	/**
	 * Convenience function to return the entry 'CHILD_TABLE' for index 'i'.
	 * 
	 * @param i the index
	 * 
	 * @return the entry 'CHILD_TABLE' for index 'i', -1 iff it needs to be recomputed (used in getChildTab{UP|DOWN|NEXT})
	 */
	@Override
	protected int getChildTab(int i){
		assert(i>=0 && i<this.childtab.length);
		return this.childtab[i]; 
	}
	
	/**
	 * Convenience function to set the entry 'CHILD_TAB(i)=j'.
	 * 
	 * @param i the index
	 * @param j the value of 'CHILD_TAB(i)' to set
	 */
	@Override
	protected void setChildTab(int i, int j){
		assert i>=0 && i<this.suftab.length && j>=0 && j<this.suftab.length : "i and j must be in interval [0,sequenceLength]";
		this.childtab[i] = j;
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
	
	/**
	 * Tests whether the child entry for suffix 'i' is an up index.
	 * 
	 * @param i the child to test
	 * 
	 * @return true iff the child entry for suffix 'i' is an up index
	 */
	@Override
	protected boolean containsUpIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		if(i == sequenceLength) return true;
		else return getLcpTab(i) > getLcpTab(i+1);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is a down index.
	 * 
	 * @param i the child to test
	 * 
	 * @return true iff the child entry for suffix 'i' is a down index
	 */
	@Override
	protected boolean containsDownIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		return getLcpTab(getChildTab(i)) > getLcpTab(i);
	}
	
	/**
	 * Tests whether the child entry for suffix 'i' is a next index.
	 * 
	 * @param i the child to test
	 * 
	 * @return true iff the child entry for suffix 'i' is the next index
	 */
	@Override
	protected boolean containsNextIndex(int i){
		assert(i>=0 && i<this.childtab.length);
		return i != sequenceLength - 1 && getChildTab(i) != Integer.MIN_VALUE && getLcpTab(getChildTab(i)) == getLcpTab(i) && getChildTab(i) > i;
	}
	
	/**
	 * Returns the 'up' value stored in childtab(i) - including integrity check.
	 * 
	 * @param i the index
	 * 
	 * @return the 'up' value stored in 'childtab(i)'
	 */
	@Override
	protected int getChildTabUP(int i){
		assert containsUpIndex(i-1) : "Child table does not contain an up value for index " + i + "!";
		return getChildTab(i-1);
	}
	
	/**
	 * Returns the 'down' value stored in 'childtab(i)' - including integrity check.
	 * 
	 * @param i the index
	 * 
	 * @return the 'down' value stored in 'childtab(i)'
	 */
	@Override
	protected int getChildTabDown(int i){
		assert containsDownIndex(i) : "up?: " + containsUpIndex(i)  + " down?: " + containsDownIndex(i) + " next?: " + containsNextIndex(i) + "\t" + printStatus(i) + "Child table does not contain a down value for index " + i + "!";
		return getChildTab(i);
	}
	
	/**
	 * Returns the 'next' value stored in 'childtab(i)' - including integrity check.
	 * 
	 * @param i the index
	 * 
	 * @return the 'down' value stored in 'childtab(i)'
	 */
	@Override
	protected int getChildTabNext(int i){
		assert containsNextIndex(i) : "Child table does not contain a next value for index " + i + "!";
		return getChildTab(i);
	}
	
	private String printStatus(int i){
		return getChildTab(i) + "\t" + getLcpTab(getChildTab(i)) + "\t" + printSuffices(i-1, getChildTab(i) + 1, 15);
	}
	
	private String printSuffices(int from, int to, int length){
		StringBuffer buffer = new StringBuffer();
		buffer.append("\n");
		for(int i=from; i<=to; i++){
			buffer.append(this.getSequence().substring(getSufTab(i), getSufTab(i) + length) + "\n");
		}
		return buffer.toString();
	}
	
	/**
	 * Returns the maximum size of a query string.
	 */
	public int getMaxWordSize() {
		return this.maxWordSize;
	}

	/**
	 * returns whether this index includes the region which was used for primer design of 'just background sequence'
	 */
	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	/**
	 * Deserializes an index of this datatype from file 'filename'.
	 * 
	 * @param filename the filename to deserialize the index from
	 * 
	 * @return an index of this datatype
	 */
	public static EnhancedSuffixArrayFatOpt deserialize(String filename){
		EnhancedSuffixArrayFatOpt result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (EnhancedSuffixArrayFatOpt) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IOException{		
		EnhancedSuffixArrayFatOpt.printStatus = false;
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			EnhancedSuffixArrayFatOpt testIndex;
				System.out.println("Reading sequence " + file);
				{					
					System.out.println("Creating EnhancedSuffixArray for first sequence in file " + file);
					
					System.out.println("Constructing EnhancedSuffixArray");
					//sequence = "ACAAACATAT";
					testIndex = new EnhancedSuffixArrayFatOpt(new File(file));
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
				testIndex.serialize(new File(file + ".ESAFatOpt.esaidx"));
				System.out.println(" - done in " + timer.getTimeString());
				System.out.println("Index construction done in " + timer.getTotalTimestring());
				//testIndex.printArray();
				System.out.println("Statistics:");
				System.out.println("Sequence length: " + testIndex.getSequenceLength());
				HashMap<Character, Integer> stat = testIndex.getStatistics();
				for(Character x : stat.keySet()){
					System.out.println(x + ": " + stat.get(x));
				}
		}
		else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayFatOpt index = EnhancedSuffixArrayFatOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for string " + query);
			int result = 0;
			for(int i=0; i<numQueries; i++){
				result = index.findHitCount(query);
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
			//index.printESA();
		}
		else if(mode.equals("randomquery")){
			int numQueries = Integer.parseInt(args[2]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayFatOpt index = EnhancedSuffixArrayFatOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for random string of length " + 18 + "-" + 30);
			int result = 0;
			for(int i=0; i<numQueries; i++){
					result = index.findHitCount(SeqTools.getRandomPrimerSequence(18, 30));
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
			System.out.print("Creating " + numQueries + " random primers itself took");
			for(int i=0; i<numQueries; i++){
				SeqTools.getRandomPrimerSequence(18, 30);
			}
			System.out.println(" " + timer.getTimeString());
		}
		else{
			System.out.println("Usage: EnhancedSuffixArrayFatOpt query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: EnhancedSuffixArrayFatOpt create <serializedIndexFile>");
		}
	}
}
