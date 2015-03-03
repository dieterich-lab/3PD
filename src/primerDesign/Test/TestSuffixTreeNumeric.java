package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SuffixTree;

import primerDesign.util.SimpleTimer;

/**
 * Creates a suffix tree for a sequence.
 * 
 * @author froehler
 *
 */
public class TestSuffixTreeNumeric extends SuffixTree {
	/**
	 * 
	 */
	private static final long serialVersionUID = -3563849880916241991L;
	private boolean includesScanRegion;
	private String sequence;

	/**
	 * Initializes a new suffix tree with the specified alphabet.
	 *
	 */
	public TestSuffixTreeNumeric(){
		super((FiniteAlphabet) AlphabetManager.alphabetForName("DNA"));
		this.includesScanRegion = false;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#createSuffixTree(java.lang.String, int)
	 */
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion){
		try{
			this.sequence = sequence;
			Alphabet alpha = AlphabetManager.alphabetForName("DNA");
			SymbolTokenization charTokenization = alpha.getTokenization("token");
			this.addSymbols(new SimpleSymbolList(charTokenization, sequence), maxWordSize);
			this.includesScanRegion = includesScanRegion;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchInSuffixTree(java.lang.String)
	 */
	public float searchMatchPositionsInIndex(String searchString){
		float result = 0;
		try{
			SimpleSymbolList sequenceString = new SimpleSymbolList(AlphabetManager.alphabetForName("DNA").getTokenization("token"), searchString);
			SuffixNode currentNode = this.getRoot();

			for(int i=1; i<= sequenceString.length(); i++){
				currentNode = this.getChild(currentNode, sequenceString.symbolAt(i));
				//if(currentNode.isTerminal()) 
				result = currentNode.getNumber();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		return result;
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
	 * Main class for testing purposes.
	 * 
	 * @param args args[1] input file
	 * 
	 * @throws NoSuchElementException
	 * @throws BioException
	 * @throws IOException
	 */
	public static void main(String[] args) throws NoSuchElementException, BioException, IOException{
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		System.out.println("Reading sequence");
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		String sequence = "";
		String line;
		while((line = reader.readLine()) != null){
			if(!line.startsWith(">", 0)){
				sequence += line.trim();
			}
		}
		System.out.println("Read sequence - " + timer.getTimeString());
		
		System.out.println("Constructing suffix tree");
		TestSuffixTreeNumeric testTree = new TestSuffixTreeNumeric();
		testTree.createIndex(sequence, 30, true);
		System.out.println("Suffix tree construction complete - " + timer.getTimeString());
		
		String querystring = "GGTCGTGGTCAGTTGTTG";
		
		System.out.println("Querying suffix tree");
		float result = testTree.searchMatchPositionsInIndex(querystring);
		System.out.println("Queried suffix tree - " + timer.getTimeString());
		System.out.println("Result: " + result);
		System.out.println("Queried suffix tree positions - " + timer.getTimeString());
		
		System.out.println("Querying suffix tree again");
		result = testTree.searchMatchPositionsInIndex(querystring);
		System.out.println("Queried suffix tree again - " + timer.getTimeString());
		System.out.println("Result: " + result);
		System.out.println("Queried suffix tree positions - " + timer.getTimeString());
	}
}
