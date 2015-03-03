/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.MultiSeqIndex;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;

/**
 * Encapsulates at least one memory mapped enhanced suffix array constructed from a file containing at least one sequence.
 * 
 * @author froehler
 *
 */
public class MultiSeqMemoryMappedESAIndex extends MultiSeqIndex implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private HashMap<String, EnhancedSuffixArrayFatOptMemoryMapped> indices;
	private static boolean printStatus = false;
	
	public MultiSeqMemoryMappedESAIndex(File file){
		super(file);
		//this.indices = new HashMap<String, EnhancedSuffixArrayFatOptMemoryMapped>();
	}

	/**
	 * Dummy implementation to fit the interface :-)
	 */
	public void createIndex() {
		try{
			SlimFastaParser parser = new SlimFastaParser(this.fastaFile);
			SimpleContig currentContig;
			SimpleTimer timer = new SimpleTimer();
			EnhancedSuffixArrayFatOptMemoryMapped currentIndex;
			
			while(parser.hasNextContig()){
				currentContig = parser.parseNextContigIgnoreCase();
				if(MultiSeqMemoryMappedESAIndex.printStatus) System.out.println("Creating index for contig: " + currentContig.getID() + " and size " + currentContig.getSequenceLength());
				currentIndex = new EnhancedSuffixArrayFatOptMemoryMapped(currentContig.getSequence(), this.fastaFile.getPath() + "_" + currentContig.getID().replaceAll("[ \t]+", "_"));
				currentIndex.createIndex();
				//this.indices.put(currentContig.getName(), currentIndex);
				super.putIndex(currentContig.getID(), currentIndex);
				
				if(MultiSeqMemoryMappedESAIndex.printStatus) System.out.println("\tCreated index for contig: " + currentContig.getID() + " and size: " + currentContig.getSequenceLength() + " in " + timer.getTimeString());
			}
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#getSequence()
	 */
	public String getSequence() {
		throw new IllegalStateException("Retrieving the sequence string is not supported by this kind of index!");
	}
	
	public static MultiSeqMemoryMappedESAIndex deserialize(File filename){
		MultiSeqMemoryMappedESAIndex result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (MultiSeqMemoryMappedESAIndex) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getContig()
	 */
	public SimpleContigImpl[] getContig() {
		//SimpleContig[] contigs = new SimpleContig[this.indices.size()];
		SimpleContigImpl[] contigs = new SimpleContigImpl[super.getNumIndices()];
		int i=0;
		DNASequenceIndex currentIndex;
		//for(String contig : this.indices.keySet()){
		for(String contig : super.getIndicesKeyset()){
			//currentIndex = this.indices.get(contig);
			currentIndex = super.getIndex(contig);
			contigs[i] = currentIndex.getContig()[0];
			assert(currentIndex.getContig().length == 1);
		}
		return contigs;
	}
	
	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		printStatus = true;
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			MultiSeqMemoryMappedESAIndex index = new MultiSeqMemoryMappedESAIndex(new File(file));
			
			System.out.println("Creating single indices from file " + file);
			index.createIndex();
			System.out.println("Creating single indices from file " + file + " - done in " + timer.getTimeString());

			System.out.print("Serializing index");
			index.serialize(new File(file + ".ESAFatOptMemoryMapped.esaidx"));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Index construction done in " + timer.getTotalTimestring());
		}
		else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			MultiSeqIndex index = MultiSeqMemoryMappedESAIndex.deserialize(new File(file));
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
			MultiSeqMemoryMappedESAIndex index = MultiSeqMemoryMappedESAIndex.deserialize(new File(file));
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
			System.out.println("Usage: MultiSeqMemoryMappedESAIndex query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: MultiSeqMemoryMappedESAIndex create <serializedIndexFile>");
		}
	}
}
