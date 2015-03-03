/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.text.NumberFormat;

import primerDesign.algo.KoAluruSuffixSort;
import primerDesign.algo.LinearTimeLCP;
import primerDesign.dsc.nio.MemoryMappedIntFile;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;

/**
 *  This class implements the enhanced suffix array as proposed by Kurtz et.al. 2004 using memory-mapped integer tables.
 * 
 *  Reference: Abouelhoda, Kurtz, Ohlebusch: Replacing suffix trees with enhanced suffix arrays. J. Discr. Algor. 2 (2004) p.53-86
 *  
 *  As of now (2/2008), memory mapped integer files can take up to ~400.000.000 integer values.
 * 
 * @author froehler
 *
 */
public class EnhancedSuffixArrayFatOptMemoryMapped extends primerDesign.dsc.indexStructures.esa.EnhancedSuffixArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String name;
	private MemoryMappedIntFile suftab;
	private MemoryMappedIntFile lcptab;
	private MemoryMappedIntFile childtab;
	private static final int MAX_CONTIG_SIZE = MemoryMappedIntFile.getMaxContigSize();
	protected static boolean printStatus = true;

	protected EnhancedSuffixArrayFatOptMemoryMapped(String sequence, String name){
		super(sequence, name);
		if(sequence.length() > EnhancedSuffixArrayFatOptMemoryMapped.MAX_CONTIG_SIZE) throw new IllegalArgumentException("The maximum size of the sequence in that index is: " + EnhancedSuffixArrayFatOptMemoryMapped.getMaxContigSize());
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	protected EnhancedSuffixArrayFatOptMemoryMapped(char[] sequence, String name){
		super(sequence, name);
		if(sequence.length > EnhancedSuffixArrayFatOptMemoryMapped.MAX_CONTIG_SIZE) throw new IllegalArgumentException("The maximum size of the sequence in that index is: " + EnhancedSuffixArrayFatOptMemoryMapped.getMaxContigSize());
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		if(printStatus){
			System.gc();
			System.out.println("Before all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	protected EnhancedSuffixArrayFatOptMemoryMapped(File file) throws IOException{
		super(file);
		if(this.sequenceLength > EnhancedSuffixArrayFatOptMemoryMapped.MAX_CONTIG_SIZE) throw new IllegalArgumentException("The maximum size of the sequence in that index is: " + EnhancedSuffixArrayFatOptMemoryMapped.getMaxContigSize());
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
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		SimpleTimer timer = new SimpleTimer();
		try{
			this.suftab = new MemoryMappedIntFile(new File(transformName(getName()) + "_suftab") , (sequenceLength+1));
			this.lcptab = new MemoryMappedIntFile(new File(transformName(getName()) + "_lcptab") , (sequenceLength+1));
			this.childtab = new MemoryMappedIntFile(new File(transformName(getName()) + "_childtab") , (sequenceLength+1));
		}
		catch(IOException e){
			e.printStackTrace();
		}
		if(printStatus){
			System.out.print("Initializing memory mapped files");
		}
		for(int i=0; i< sequenceLength; i++){
			setSufTab(i, i);
		}

		this.suftab.put(sequenceLength, Integer.MAX_VALUE);
		if(printStatus){
			System.out.println(" - done in " + timer.getTimeString());
		}
		
		// sort table
		if(printStatus){
			System.gc();
			System.out.print("Before sort: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
			timer.getTimeString();
		}
		{
			int[] result;
			{
				KoAluruSuffixSort sorter = new KoAluruSuffixSort();
				result = sorter.getSuffixArray(this.sequence);
				for(int i=0; i<result.length; i++){
					this.suftab.put(i, result[i]);
				}
			}
			//quicksortSuffices(0, sequenceLength-1);
			if(printStatus){
				System.out.println(" - sorted " + format.format(this.sequenceLength) + "bp in " + timer.getTimeString());
				System.gc();
				System.out.print("Before LCP Table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
				timer.getTimeString();
			}
			{
				int[] lcp = LinearTimeLCP.getLCP(this.sequence, result);
				for(int i=0; i<lcp.length; i++){
					this.lcptab.put(i, lcp[i]);
				}
			}
		}
		//computeLCPTable();
		if(printStatus){
			System.out.println(" - lcp table computed in " + timer.getTimeString());
			System.gc();
			System.out.print("Before child table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
			timer.getTimeString();
		}
		computeChildTable();
		if(printStatus){
			System.out.println(" - child table computed in " + timer.getTimeString());
			System.gc();
			System.out.print("Before bucket table: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
			timer.getTimeString();
		}
		computeBucketTable();
		if(printStatus){
			System.out.println(" - bucket table computed in " + timer.getTimeString() + " size: " + this.bucketTab.size());
			System.gc();
			System.out.println("After all: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#initChildTable()
	 */
	@Override
	protected void initChildTable() {
		for(int i=0; i<this.childtab.length(); i++) this.childtab.put(i, Integer.MIN_VALUE);
	}
	
	protected void initLcpTab(){
		for(int i=0; i<this.lcptab.length(); i++) this.lcptab.put(i, Integer.MIN_VALUE);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsUpIndex(int)
	 */
	@Override
	protected boolean containsUpIndex(int i) {
		assert(i>=0 && i<this.childtab.length());
		if(i == sequenceLength) return true;
		else return getLcpTab(i) > getLcpTab(i+1);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsDownIndex(int)
	 */
	@Override
	protected boolean containsDownIndex(int i) {
		assert(i>=0 && i<this.childtab.length());
		return getLcpTab(getChildTab(i)) > getLcpTab(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#containsNextIndex(int)
	 */
	@Override
	protected boolean containsNextIndex(int i) {
		assert(i>=0 && i<this.childtab.length());
		//return getChildTab(i) != Integer.MIN_VALUE && getLcpTab(getChildTab(i)) == getLcpTab(i);
		return i != sequenceLength - 1 && getChildTab(i) != Integer.MIN_VALUE && getLcpTab(getChildTab(i)) == getLcpTab(i) && getChildTab(i) > i;
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTab(int)
	 */
	@Override
	protected int getChildTab(int i) {
		assert(i>=0 && i<this.childtab.length());
		return this.childtab.get(i); 
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabUP(int)
	 */
	@Override
	protected int getChildTabUP(int i) {
		assert containsUpIndex(i-1) : "Child table does not contain an up value for index " + i + "!";
		return getChildTab(i-1);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabDown(int)
	 */
	@Override
	protected int getChildTabDown(int i) {
		assert containsDownIndex(i) : "up?: " + containsUpIndex(i)  + " down?: " + containsDownIndex(i) + " next?: " + containsNextIndex(i) + "\t" + printStatus(i) + "Child table does not contain a down value for index " + i + "!";
		return getChildTab(i);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabNext(int)
	 */
	@Override
	protected int getChildTabNext(int i) {
		assert containsNextIndex(i) : "Child table does not contain a next value for index " + i + "!";
		return getChildTab(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getChildTabLength()
	 */
	@Override
	protected int getChildTabLength() {
		return this.childtab.length();
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getLcpTab(int)
	 */
	@Override
	protected int getLcpTab(int i) {
		assert(i>=0 && i<this.lcptab.length());
		return this.lcptab.get(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getLcpTabLength()
	 */
	@Override
	protected int getLcpTabLength() {
		return this.lcptab.length();
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getSufTab(int)
	 */
	@Override
	protected int getSufTab(int i) {
		assert(i>=0 && i<this.suftab.length());
		return this.suftab.get(i);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#getSufTabLength()
	 */
	@Override
	protected int getSufTabLength() {
		return this.suftab.length();
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setChildTab(int, int)
	 */
	@Override
	protected void setChildTab(int i, int j) {
		assert i>=0 && i<this.suftab.length() && j>=0 && j<this.suftab.length() : "i and j must be in interval [0,sequenceLength]";
		this.childtab.put(i, j);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setLcpTab(int, int)
	 */
	@Override
	protected void setLcpTab(int i, int j) {
		assert(i>=0 && j>=0 && i<this.lcptab.length() && j<=this.sequenceLength);
		this.lcptab.put(i, j);
	}

	/* (non-Javadoc)
	 * @see primerDesign.Test.esa.EnhancedSuffixArray#setSufTab(int, int)
	 */
	@Override
	protected void setSufTab(int i, int j) {
		assert(i>=0 && j>=0 && i<this.suftab.length() && j<this.suftab.length());
		this.suftab.put(i, j);
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
	
	private static int getMaxContigSize(){
		return EnhancedSuffixArrayFatOptMemoryMapped.MAX_CONTIG_SIZE;
	}
	
	public static EnhancedSuffixArrayFatOptMemoryMapped deserialize(String filename){
		EnhancedSuffixArrayFatOptMemoryMapped result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (EnhancedSuffixArrayFatOptMemoryMapped) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	private String transformName(String name){
		return name.replaceAll("[ /]", "_");
	}
	
	public static void main(String[] args) throws IOException{		
		printStatus = true;
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			EnhancedSuffixArrayFatOptMemoryMapped testIndex;
			System.out.println("Reading first sequence in file " + file);
			{					
				System.out.println("Creating EnhancedSuffixArray for first sequence in file " + file);
				
				System.out.println("Constructing EnhancedSuffixArray");
				//sequence = "ACAAACATAT";
				testIndex = new EnhancedSuffixArrayFatOptMemoryMapped(new File(file));
			}
			System.gc();
			System.out.println("Reading first sequence in file " + file +  ") took " + timer.getTimeString());
					
			//System.gc();
			testIndex.createIndex();
			System.out.println("EnhancedSuffixArray construction complete - " + timer.getTimeString());
			System.out.print("Serializing index");
			testIndex.serialize(new File(file + ".ESAFatOptMemoryMapped.esaidx"));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Index construction done in " + timer.getTotalTimestring());
			System.out.println(testIndex.getName());
		}
		else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayFatOptMemoryMapped index = EnhancedSuffixArrayFatOptMemoryMapped.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for string " + query);
			int result = 0;
			for(int i=0; i<numQueries; i++){
					result = index.findHitCount(query);
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
		}
		else if(mode.equals("randomquery")){
			int numQueries = Integer.parseInt(args[2]);
			System.out.print("Deserializing index");
			EnhancedSuffixArrayFatOptMemoryMapped index = EnhancedSuffixArrayFatOptMemoryMapped.deserialize(file);
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
			System.out.println("Usage: EnhancedSuffixArrayFatOptMemoryMapped query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: EnhancedSuffixArrayFatOptMemoryMapped create <serializedIndexFile>");
		}
	}
}
