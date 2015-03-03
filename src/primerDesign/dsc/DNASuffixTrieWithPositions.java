package primerDesign.dsc;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Pattern;

import primerDesign.util.Constants;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

public class DNASuffixTrieWithPositions extends SuffixTrieWithPositions implements DNASequenceIndex {
	
	private boolean includesScanRegion;
	private String sequence;
	private int maxWordSize;
	
	public DNASuffixTrieWithPositions(){
		super(DNASuffixTrieWithPositions.initObjectArrayList());
		this.includesScanRegion = false;
		this.sequence = "";
		this.maxWordSize = 0;
	}
	
	public DNASuffixTrieWithPositions(ObjectArrayList alphabet) {
		super(alphabet);
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

	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
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

	public Integer[] searchMatchPositionsInIndex(String searchString) {
		searchString = searchString.toUpperCase();
		HashSet<Integer> tempResult = new HashSet<Integer>();
		HashSet<Integer> result = new HashSet<Integer>();
		for(int i=0; i<searchString.length(); i+=maxWordSize){
			try{
				SuffixNode currentNode = this.getRoot();
				
				for(int j=i; j< Math.min(i+maxWordSize, searchString.length()); j++){
					currentNode = this.getChild(currentNode, searchString.charAt(j));
					tempResult = currentNode.getPositions();
				}
				if(i==0) result = tempResult;
				else result = DNASuffixTrieWithPositions.compareMatchSets(result, tempResult, i);
				if(result.isEmpty()) break;
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		return (Integer[]) result.toArray(new Integer[result.size()]);
	}

	public int searchNbMatchesInIndex(String searchString) {
		return searchMatchPositionsInIndex(searchString.toUpperCase()).length;
	}
	
	public static HashSet<Integer> compareMatchSets(HashSet<Integer> first, HashSet<Integer> second, int offset){
		HashSet<Integer> commonElements = new HashSet<Integer>();
		Iterator iter = first.iterator();
		while(iter.hasNext()){
			int position = (Integer) iter.next();
			if(second.contains(position + offset)) commonElements.add(position);
		}
		return commonElements;
	}
	
	public boolean includesScanRegion(){
		return this.includesScanRegion;
	}
	
	public String getSequence(){
		return this.sequence;
	}

	public int getMaxWordSize() {
		return this.maxWordSize;
	}
	
	public boolean equals(Object otherTree){
		DNASuffixTrieWithPositions other = (DNASuffixTrieWithPositions) otherTree;
		return this.includesScanRegion == other.includesScanRegion && this.maxWordSize == other.maxWordSize 
		&& this.sequence.equals(other.sequence) && this.getAlphabet().equals(other.getAlphabet());
	}
	
	public static void main(String[] args) throws IOException{
		//String[] sequences = {"100kb.dna", "200kb.dna", "400kb.dna", "800kb.dna", "../src/yeast-chrXII.dna", "Cel/Celegans.fa", "Ppa/Ppacificus.fa"};
		String[] sequences = {"../src/yeast-chrXII.dna"};
		String path = "/export/Sebastian/PrimerDesign/Testsequenzen/";
		SimpleTimer timer = new SimpleTimer();
		DNASequenceIndex testIndex;
		BufferedReader reader;
		String sequence;
		String line;
		Pattern fastaStart = Pattern.compile(">.*");
		for(int i=0; i<sequences.length; i++){
			System.out.println("Reading sequence " + sequences[i]);
			{
				reader = new BufferedReader(new FileReader(path + sequences[i]));
				StringBuffer seq = new StringBuffer();
				while((line = reader.readLine()) != null){
					if(!fastaStart.matcher(line).matches()){
						seq.append(line.trim());
					}
				}
				sequence = seq.toString();
			}
			System.out.println("Reading sequence " + sequences[i] + " (length: " + sequence.length() + ") took " + timer.getTimeString());
			
			System.out.println("Creating DNASuffixTrieWithPositions for sequence " + sequences[i]);
			
			System.out.println("Constructing DNASuffixTrieWithPositions");
			Runtime runtime = Runtime.getRuntime();
			NumberFormat format = NumberFormat.getInstance();
			System.gc();
			System.out.println("Before DNASuffixTrieWithPositions creation: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
			
			testIndex = new DNASuffixTrieWithPositions();
			testIndex.createIndex(sequence, 10, true);
			System.out.println("After DNASuffixTrieWithPositions creation: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
			System.out.println("DNASuffixTrieWithPositions construction complete - " + timer.getTotalTimestring());
				String query = "ATG";
				System.out.println("Querying DNASuffixTrieWithPositions for: " + query);
				System.out.println(testIndex.searchNbMatchesInIndex(query) + " matches");
				System.out.println("Queried DNASuffixTrieWithPositions - " + timer.getTimeString());
				System.out.println();
		}
	}
}
