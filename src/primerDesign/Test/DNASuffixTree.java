package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.CharacterTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.util.Constants;
import primerDesign.util.SimpleTimer;

/**
 * Creates a suffix tree for a sequence.
 * 
 * @author froehler
 *
 */
public class DNASuffixTree extends SuffixTreeWithPositions implements DNASequenceIndex{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3563849880916241991L;
	private boolean includesScanRegion;
	private String sequence;
	private int maxWordSize;
	private SymbolTokenization tokenization;

	/**
	 * Initializes a new suffix tree with the specified alphabet.
	 *
	 */
	public DNASuffixTree(){
		super(DNASuffixTree.createExtendedAlphabet());
		this.includesScanRegion = false;
	}
	
	/**
	 * Creates a standard DNA alphabet extended by a 'masked-basepair' character.
	 * 
	 * @return a standard DNA alphabet extended by a 'masked-basepair' character
	 */
	public static FiniteAlphabet createExtendedAlphabet(){
		Set<AtomicSymbol> set = new HashSet<AtomicSymbol>();
		set.add(AlphabetManager.createSymbol("a"));
		set.add(AlphabetManager.createSymbol("t"));
		set.add(AlphabetManager.createSymbol("g"));
		set.add(AlphabetManager.createSymbol("c"));
		set.add(AlphabetManager.createSymbol(Constants.REPETITIVE_ELEMENT_CHARACTER.toLowerCase()));
		
		String alphabet_name = "masked_dna";
		FiniteAlphabet alphabet = (FiniteAlphabet) new SimpleAlphabet(set, alphabet_name);
		AlphabetManager.registerAlphabet(alphabet_name, alphabet);
		
		//FiniteAlphabet alphabet = DNATools.getDNA();
		
		return alphabet;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#createSuffixTree(java.lang.String, int)
	 */
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion){
		try{
			this.sequence = sequence.toLowerCase();
			
			Symbol a = (Symbol) AlphabetManager.createSymbol("a");
			Symbol t = (Symbol) AlphabetManager.createSymbol("t");
			Symbol g = (Symbol) AlphabetManager.createSymbol("g");
			Symbol c = (Symbol) AlphabetManager.createSymbol("c");
			Symbol n = (Symbol) AlphabetManager.createSymbol(Constants.REPETITIVE_ELEMENT_CHARACTER.toLowerCase());
			
			CharacterTokenization tokenization = new CharacterTokenization(this.getAlphabet(), false);
			tokenization.bindSymbol(a, 'a');
			tokenization.bindSymbol(t, 't');
			tokenization.bindSymbol(g, 'g');
			tokenization.bindSymbol(c, 'c');
			tokenization.bindSymbol(n, 'n');
			
			//SymbolTokenization tokenization = this.getAlphabet().getTokenization("token"); 
			((SimpleAlphabet)this.getAlphabet()).putTokenization("token", tokenization);
			
			SimpleSymbolList list = new SimpleSymbolList(tokenization, sequence);
			
			this.addSymbols(list, maxWordSize);
			this.includesScanRegion = includesScanRegion;
			this.maxWordSize = maxWordSize;
			this.tokenization = (SymbolTokenization) tokenization;
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchInSuffixTree(java.lang.String)
	 */
	public int searchNbMatchesInIndex(String searchString){
		float result = 0;
		try{
			SimpleSymbolList sequenceString = new SimpleSymbolList(this.tokenization, searchString);
			SuffixNode currentNode = this.getRoot();

			for(int i=1; i<= sequenceString.length(); i++){
				currentNode = this.getChild(currentNode, sequenceString.symbolAt(i));
				// if(currentNode.isTerminal()) 
				result = currentNode.getNumberOfMatches();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		return (int) result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchInSuffixTree(java.lang.String)
	 */
	public Integer[] searchMatchPositionsInIndex(String searchString){
		Integer[] result = new Integer[0];
		try{
			SimpleSymbolList sequenceString = new SimpleSymbolList(this.tokenization, searchString);
			SuffixNode currentNode = this.getRoot();

			for(int i=1; i<= sequenceString.length(); i++){
				currentNode = this.getChild(currentNode, sequenceString.symbolAt(i));
				//if(currentNode.isTerminal()) 
				result = currentNode.getAllMatchPositions();
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
	 * Returns the maximum word size contained in this index.
	 * 
	 * @return the maximum word size contained in this index
	 */
	public int getMaxWordSize() {
		return this.maxWordSize;
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
		DNASequenceIndex testTree = new DNASuffixTree();
		testTree.createIndex(sequence, 20, true);
		System.out.println("Suffix tree construction complete - " + timer.getTimeString());
		
		String querystring = "GGTCGTGGTCAGTTGTTG";
		
		System.out.println("Querying suffix tree");
		float result = testTree.searchNbMatchesInIndex(querystring);
		System.out.println("Queried suffix tree - " + timer.getTimeString());
		Integer[] positions = testTree.searchMatchPositionsInIndex(querystring);
		System.out.println("Result: " + result);
		for(int i=0; i<positions.length; i++){
			System.out.println("\tPosition: " + (Integer) positions[i]);
		}
		System.out.println("Queried suffix tree positions - " + timer.getTimeString());
		
		System.out.println("Querying suffix tree again");
		result = testTree.searchNbMatchesInIndex(querystring);
		System.out.println("Queried suffix tree again - " + timer.getTimeString());
		positions = testTree.searchMatchPositionsInIndex(querystring);
		System.out.println("Result: " + result);
		for(int i=0; i<positions.length; i++){
			System.out.println("\tPosition: " + (Integer) positions[i]);
		}
		System.out.println("Queried suffix tree positions - " + timer.getTimeString());
	}
}
