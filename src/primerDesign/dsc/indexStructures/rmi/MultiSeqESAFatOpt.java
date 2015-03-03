/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.rmi.RemoteException;
import java.text.NumberFormat;
import java.util.HashMap;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.esa.EnhancedSuffixArrayFatOpt;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a multiple-sequence index using the full implementation of the optimized enhanced suffix array.
 * 
 * @author Sebastian Fršhler
 *
 */
//public class MultiSeqESAFatOpt extends UnicastRemoteObject implements DNASequenceIndex, Serializable {
public class MultiSeqESAFatOpt implements DNASequenceIndex, Serializable {

	private static final long serialVersionUID = 6153908947620571392L;

	private HashMap<String, EnhancedSuffixArrayFatOpt> indices;
	private boolean debug = true;
	private File file;
	
	/**
	 * Initializes a new MultiSeqESAFatOpt index.
	 * 
	 * @throws RemoteException
	 */
	public MultiSeqESAFatOpt(File file) throws RemoteException {
		super();
		this.file = file;
	}


	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#createIndex(java.io.File)
	 */
	public void createIndex() {
		this.indices = new HashMap<String, EnhancedSuffixArrayFatOpt>();
		SimpleTimer timer = new SimpleTimer();
		try{
			SlimFastaParser parser = new SlimFastaParser(this.file);
			EnhancedSuffixArrayFatOpt currentIndex;
			
			while(parser.hasNextContig()){
				SimpleContig currentContig = parser.parseNextContigIgnoreCase();
				currentIndex = new EnhancedSuffixArrayFatOpt(currentContig.getSequence(), currentContig.getID());
				currentIndex.createIndex();
				this.indices.put(currentContig.getID(), currentIndex);
				
			if(debug) System.out.println("\tCreated index for contig \"" + currentContig.getID() + "\" done in " + timer.getTimeString());
			}
		}catch(IOException e){
			e.printStackTrace();
		}
	}

	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#findHitCount(java.lang.String)
	 */
	public int findHitCount(String sequence) {
		int count = 0;
		// count matches in all contigs
		for(String contig : this.indices.keySet()){
			count += this.indices.get(contig).findHitPositions(sequence.toUpperCase()).size();
		}
		return count;
	}

	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#findHitPositions(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence) {
		ObjectArrayList result = new ObjectArrayList();
		ObjectArrayList currentList;
		// search matches in all contigs
		for(String contig : this.indices.keySet()){
			// retrieve match positions in current contig
			currentList = this.indices.get(contig).findHitPositions(sequence.toUpperCase());
			for(int i=0; i<currentList.size(); i++){
				result.add(currentList.get(i));
			}
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#serialize(java.io.File)
	 */
	public void serialize(File file) {
		try{
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#deserialize(java.io.File)
	 */
	public static MultiSeqESAFatOpt deserialize(File file) {
		MultiSeqESAFatOpt result = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file));
			result = (MultiSeqESAFatOpt) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public String getName(){
		return this.file.getName();
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getContig()
	 */
	public SimpleContigImpl[] getContig() {
		SimpleContigImpl[] contigs = new SimpleContigImpl[this.indices.size()];
		int i=0;
		EnhancedSuffixArrayFatOpt currentIndex;
		for(String contig : this.indices.keySet()){
			currentIndex = this.indices.get(contig);
			contigs[i] = currentIndex.getContig()[0];
			assert(currentIndex.getContig().length == 1);
		}
		return contigs;
	}
	
	/**
	 * Returns the combined sequence length of all sequences of this index.
	 * 
	 * @return the length of the sequence
	 */
	public int getSequenceLength() {
		int result = 0;
		for(String contig : this.indices.keySet()){
			result += this.indices.get(contig).getContig()[0].getSequenceLength();
		}
		return result;
	}
	
	/**
	 * Returns some statistics of single character frequency in the sequence - over all contigs.
	 * 
	 * @return a hasmap of: character to integer frequency counts
	 */
	public HashMap<Character, Integer> getStatistics(){
		HashMap<Character, Integer> result = new HashMap<Character, Integer>();
		
		for(String contig : this.indices.keySet()){
			HashMap<Character, Integer> temp = this.indices.get(contig).getStatistics();
			for(Character x : temp.keySet()){
				if(result.containsKey(x)){
					result.put(x, result.get(x) + temp.get(x));
				}
				else{
					result.put(x, temp.get(x));
				}
			}
		}
		
		return result;
	}
	
	/**
	 * Prints a textual representation of the current suffix arrays to STDOUT.
	 *
	 */
	public void printArrays(){
		for(String contig : this.indices.keySet()){
			System.out.println("\nIndex for contig: " + contig);
			this.indices.get(contig).printESA();
		}
	}
	
	public static void main(String[] args) throws IOException{
		String mode = args[0];
		File file = new File(args[1]);
		SimpleTimer timer = new SimpleTimer();
		NumberFormat format = NumberFormat.getInstance();
		
		if(mode.equals("create")){
			System.out.println("Creating index for sequence file: " + file);
			MultiSeqESAFatOpt index = new MultiSeqESAFatOpt(file);
			index.createIndex();
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.print("Serializing index");
			index.serialize(new File(file.toString() + ".MultiSeqESAFatOpt.esaidx"));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Total runtime: " + timer.getTotalTimestring());
			
		}else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			MultiSeqESAFatOpt index = MultiSeqESAFatOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.print("Querying index " + format.format(numQueries) + " times for query " + query);
			int matches = 0;
			for(int i=0; i<numQueries; i++){
				matches = index.findHitCount(query);
				matches = index.findHitPositions(query).size();
			}
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Found " + format.format(matches) + " matches to query " + query);
			
//			for(int i=0; i<index.findHitCount(query); i++){
//				System.out.print(((IndexHitImpl)index.findHitPositions(query).get(i)).getPosition() + " ");
//			}
			//index.printArrays();
		}
		else if(mode.equals("randomquery")){
			int numQueries = Integer.parseInt(args[2]);
			System.out.print("Deserializing index");
			MultiSeqESAFatOpt index = MultiSeqESAFatOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + format.format(numQueries) + " times for random string of length " + 18 + "-" + 30);
			int result = 0;
			for(int i=0; i<numQueries; i++){
					result = index.findHitCount(SeqTools.getRandomPrimerSequence(18, 30));
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
			System.out.print("Creating " + format.format(numQueries) + " random primers itself took");
			for(int i=0; i<numQueries; i++){
				SeqTools.getRandomPrimerSequence(18, 30);
			}
			System.out.println(" " + timer.getTimeString());
		}else if(mode.equals("fasta")){
			System.out.print("Deserializing index");
			MultiSeqESAFatOpt index = MultiSeqESAFatOpt.deserialize(file);
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index with sequences from file: " + args[2] + "\n");
			SlimFastaParser parser = new SlimFastaParser(new File(args[2]));
			SimpleContig contig;
			while(parser.hasNextContig()){
				contig = parser.parseNextContigIgnoreCase();
				System.out.println(contig.getID() + " " + index.findHitCount(new String(contig.getSequence())) + " hits");
			}
		}
		else{
			System.out.println("Unknown option: " + args[0]);
			System.exit(1);
		}
	}

}
