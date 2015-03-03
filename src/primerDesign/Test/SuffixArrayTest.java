package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

public class SuffixArrayTest implements DNASequenceIndex {
	private ObjectArrayList array;
	private int minWordLength;
	private int maxWordLength;
	private boolean includesScanRegion;
	private String sequence;
	
	/**
	 * Creates a suffix array of the sequence 'sequence'.
	 * 
	 * All suffices in the array are at least 'minWordLength' and at most 'maxWordLength' long!
	 * 
	 * @param sequence
	 * @param maxWordLength
	 * @param minWordLength
	 */
	public void createSuffixArray(String sequence, int maxWordLength, int minWordLength){
		if(maxWordLength > 0 && maxWordLength >= minWordLength) this.maxWordLength = maxWordLength;
		else throw new IllegalArgumentException("MaxWordLength has to be > 0 and >= MinWordLength!");
		if(minWordLength > 0 && minWordLength <= maxWordLength) this.minWordLength = minWordLength;
		else throw new IllegalArgumentException("MinWordLength has to be > 0 and <= MaxWordLength!");
		if(sequence.length() < 1) throw new IllegalArgumentException("Sequence length has to be > 0!");
		
		this.array = new ObjectArrayList();
		
		for(int i=0; i<= sequence.length()-minWordLength; i++){
			for(int j=minWordLength; j<= maxWordLength;j++){
				if(i+j <= sequence.length()){
					SuffixArrayElement currentElement = new SuffixArrayElement(sequence.substring(i, i+j).toUpperCase(), i);
					int position = array.indexOf(currentElement, true);
					if(position >= 0){
						((SuffixArrayElement) array.get(position)).addPosition(i);
					}
					else array.add(currentElement);
				}
			}
		}
		array.sort();
	}
		
	/**
	 * Scans a suffix array for the occurrence of suffix 'suffix'.
	 * 
	 * @param suffix the suffix to search in the suffix array
	 */
	public Integer[] findSuffix(String suffix){
		SuffixArrayElement searchElement = new SuffixArrayElement(suffix.toUpperCase(), Integer.MAX_VALUE);

		int position = array.binarySearch(searchElement);
		if(position >= 0){
			SuffixArrayElement matchingElement = (SuffixArrayElement) array.get(position);
			return matchingElement.getPositions();
		}
		else return new Integer[0];
	}
	
	/**
	 * @return the maxWordLength
	 */
	public int getMaxWordLength() {
		return maxWordLength;
	}

	/**
	 * @return the minWordLength
	 */
	public int getMinWordLength() {
		return minWordLength;
	}
	
	private class SuffixArrayElement implements Comparable{
		private String sequence;
		private HashSet positions;
		
		public SuffixArrayElement(String suffix, int position){
			if(suffix.length() > 0) this.sequence = suffix;
			else throw new IllegalArgumentException("Sequence length has to be > 0!");
			
			this.positions = new HashSet<Integer>();
			this.addPosition(position);
		}
		
		/**
		 * Adds a match position to a SuffixArrayElement.
		 * 
		 * @param position to position to add to an element
		 */
		public void addPosition(int position){
			if(position >= 0) this.positions.add(position);
			else throw new IllegalArgumentException("Position has to be positive!");
		}
		
		/**
		 * @return the match positions of this element
		 */
		public Integer[] getPositions() {
			return (Integer[]) this.positions.toArray(new Integer[this.positions.size()]);
		}
		/**
		 * @return the sequence
		 */
		public String getSequence() {
			return sequence;
		}
		
		public int compareTo(Object other) {
			return this.sequence.compareTo(((SuffixArrayElement) other).sequence);
		}
		
		public boolean equals(Object other){
			return this.sequence.equals(((SuffixArrayElement) other).sequence);
		}
	}
	
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
		this.createSuffixArray(sequence, maxWordSize, 0);
		this.includesScanRegion = includesScanRegion;
		this.sequence = sequence;
	}

	public int getMaxWordSize() {
		return this.maxWordLength;
	}

	public String getSequence() {
		return this.sequence;
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	public Integer[] searchMatchPositionsInIndex(String searchString) {
		return this.findSuffix(searchString);
	}

	public int searchNbMatchesInIndex(String searchString) {
		return this.findSuffix(searchString).length;
	}
	
	public static void main(String[] args) throws IOException{		
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		StringBuffer sequence = new StringBuffer();
		String line;
		SimpleTimer timer = new SimpleTimer();
		System.out.print("Reading sequence - ");
		while((line = reader.readLine()) != null){
			if(!line.matches(">.*")) sequence.append(line.trim());
		}
		System.out.println("done " + timer.getTimeString());
		
		System.out.print("Constructing suffix array - ");
		SuffixArrayTest suffixArray = new SuffixArrayTest();
		suffixArray.createSuffixArray(sequence.toString(), 10, 1);
		System.out.println("done " + timer.getTimeString());
		
		String queryString = "ATG";
		Integer[] matchPositions = suffixArray.findSuffix(queryString);
		System.out.print("Querying suffix array for string " + queryString + " - ");
//		for(int i=0; i<matchPositions.length; i++){
//			System.out.print(matchPositions[i] + " ");
//		}
		System.out.println("done " + timer.getTimeString());
	}
}
