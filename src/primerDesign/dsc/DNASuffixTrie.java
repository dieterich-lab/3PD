package primerDesign.dsc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;

import primerDesign.util.Constants;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Creates a suffix tree for a sequence.
 * 
 * @author Sebastian Fršhler
 *
 */
public class DNASuffixTrie extends SuffixTrie implements DNASequenceIndex{

	private static final long serialVersionUID = -3563849880916241991L;
	private boolean includesScanRegion;
	private String sequence;
	private int maxWordSize;

	/**
	 * Initializes a new suffix tree with the specified alphabet.
	 *
	 *@param alphabet the alphabet as an ObjectArrayList of characters
	 */
	public DNASuffixTrie(ObjectArrayList alphabet){
		super(alphabet);
		this.includesScanRegion = false;
		this.sequence = "";
		this.maxWordSize = 0;
	}
	
	/**
	 * Initializes a new suffix tree with the dna alphabet plus a 'masked-basepar' character.
	 *
	 */
	public DNASuffixTrie(){
		super(DNASuffixTrie.initObjectArrayList());
		this.includesScanRegion = false;
		this.sequence = "";
		this.maxWordSize = 0;
	}
	
	private static ObjectArrayList initObjectArrayList(){
		ObjectArrayList list = new ObjectArrayList();
		list.add((char) 'A');
		list.add((char) 'T');
		list.add((char) 'G');
		list.add((char) 'C');
		list.add((char) Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase().charAt(0));
		
		return list;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#createSuffixTree(java.lang.String, int)
	 */
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion){
		if(maxWordSize > sequence.length()) throw new IllegalArgumentException("maxWordSize must be <= sequence length!");
		try{
			this.sequence = sequence.toUpperCase();
			
			this.addSymbols(this.sequence, maxWordSize);
			this.includesScanRegion = includesScanRegion;
			this.maxWordSize = maxWordSize;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchInSuffixTree(java.lang.String)
	 */
	public int searchNbMatchesInIndex(String searchString){
		searchString = searchString.toUpperCase();
		float result = 0;
		if(this.maxWordSize >= searchString.length()){
			try{
				SuffixNode currentNode = this.getRoot();
	
				for(int i=0; i< searchString.length(); i++){
					currentNode = this.getChild(currentNode, searchString.charAt(i));
					// if(currentNode.isTerminal()) 
					result = currentNode.getNumberOfMatches();
				}
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		else{
			result = searchMatchPositionsInIndex(searchString).length;
		}
		return (int) result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchInSuffixTree(java.lang.String)
	 */
	public Integer[] searchMatchPositionsInIndex(String searchString){
		searchString = searchString.toUpperCase();
		ObjectArrayList positions = new ObjectArrayList();
		int stop = this.getSequence().length();
		int match = -1;
		while(match<stop){
			match = this.getSequence().indexOf(searchString, match +1);
			// if primer has no mispriming in index sequence
			if(match == -1) break;
			// else store position
			else positions.add(match);
		}		
		return (Integer[]) positions.toArray(new Integer[positions.size()]);
	}
	
	/**
	 * Returns whether the region to scan for valid primers is included in the region, this index is constructed from.
	 * 
	 * @return true iff sequence region includes primer scan region, false else
	 */
	public boolean includesScanRegion(){
		return this.includesScanRegion;
	}
	
	/**
	 * Returns the sequence this index was created on.
	 * 
	 * @return the sequence this index was created on
	 */
	public String getSequence(){
		return this.sequence;
	}
	
	/**
	 * Returns the maximum word size contained in this index.
	 * 
	 * @return the maximum word size contained in this index
	 */
	public int getMaxWordSize() {
		return this.maxWordSize;
	}
	
	public boolean equals(Object otherTree){
		DNASuffixTrie other = (DNASuffixTrie) otherTree;
		return this.includesScanRegion == other.includesScanRegion && this.maxWordSize == other.maxWordSize 
		&& this.sequence.equals(other.sequence) && this.getAlphabet().equals(other.getAlphabet());
	}
	
	/**
	 * Main class for performance-testing purposes.
	 * 
	 * @throws NoSuchElementException
	 * @throws BioException
	 * @throws IOException
	 */
	public static void main(String[] args) throws NoSuchElementException, BioException, IOException{
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String[] sequences = {"100kb.dna", "200kb.dna", "400kb.dna", "800kb.dna", "../src/yeast-chrXII.dna", "Ppa/Ppacificus.fa", "Cel/Celegans.fa"};
		String path = "/export/Sebastian/PrimerDesign/Testsequenzen/";
		int maxWordLength = 30;
		for(int i=0; i<sequences.length; i++){
			System.out.println("Reading sequence " + sequences[i]);
			BufferedReader reader = new BufferedReader(new FileReader(path + sequences[i]));
			StringBuffer sequence = new StringBuffer();
			String line;
			Pattern fastaStart = Pattern.compile(">.*");
			while((line = reader.readLine()) != null){
				if(!fastaStart.matcher(line).matches()){
					sequence.append(line.trim());
				}
			}
			System.out.println("Reading sequence " + sequences[i] + " (length: " + sequence.length() + ") took " + timer.getTimeString());
			
			for(int j=10; j<= maxWordLength; j += 10){
				System.out.println("Creating Index for sequence " + sequences[i] + " maxWordLength " + j);
				
				System.out.println("Constructing suffix tree");
				DNASequenceIndex testTree = new DNASuffixTrie();
				testTree.createIndex(sequence.toString(), j, true);
				System.out.println("Suffix tree construction complete - " + timer.getTimeString());
				String query = "ACA";
				System.out.println("Querying suffix tree for: " + query);
				System.out.println(testTree.searchNbMatchesInIndex(query) + " matches");
				System.out.println("Queryied suffix tree - " + timer.getTimeString());
				System.out.println();
			}
		}
	}
}
