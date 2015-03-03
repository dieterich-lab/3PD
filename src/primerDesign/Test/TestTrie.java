package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;

import weka.core.FastVector;

import com.sun.org.apache.xml.internal.utils.Trie;

public class TestTrie {
	private Trie trie;
	
	/**
	 * Creates a trie of substrings of sequence 'sequence' of length [from,to]
	 * 
	 * @param sequence the sequence to build the trie from
	 * @param from the minimum length of a substring contained in this trie
	 * @param to the maximum length of a substring contained in this trie
	 */
	public TestTrie(String sequence, int from, int to){
		trie = new Trie();
		for(int j=from; j<=to; j++){
			for(int i=0; i<sequence.length()-j+1; i++){
				String substring = sequence.substring(i, i+j);
				if(trie.get(substring) == null){
					FastVector positions = new FastVector();
					positions.addElement(i);
					trie.put(substring, positions);
				}
				else{
					FastVector position = (FastVector) trie.get(substring);
					position.addElement(i);
					trie.put(substring, position);
				}
			}
		}
	}
	
	/**
	 * Finds all occurrences of querystring 'querystring' in the trie.
	 * 
	 * @param querystring the querystring to be searched in the trie
	 * @return all occurrences of 'querystring' in the trie
	 */
	public FastVector findPositionsInTrie(String querystring){
		FastVector positions = (FastVector) trie.get(querystring);
		
		return positions;
	}
	
	/**
	 * Main class for testing the trie.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		try{
			long time = System.currentTimeMillis();
			System.out.println("Reading sequence:");
			BufferedReader reader = new BufferedReader(new FileReader(args[0]));
			String sequence = "";
			String line;
			while((line = reader.readLine()) != null){
				if(!line.startsWith(">", 0)){
					sequence += line.trim();
				}
			}
			System.out.println("Reading sequence took: " + (System.currentTimeMillis() - time) + "ms");
			
			time = System.currentTimeMillis();
			System.out.println("Constructing trie of input sequence");
			TestTrie trie = new TestTrie(sequence, 20, 30);
			time = System.currentTimeMillis();
			System.out.println("Constructing trie of input sequence took " + (System.currentTimeMillis() - time));
			
			String querystring = "TTTTTTTTCGGTTTGCGGG";
			time = System.currentTimeMillis();
			System.out.println("Scanning sequence for query string " + querystring);
			FastVector positions = trie.findPositionsInTrie(querystring);
			System.out.println("Scanning sequence for query string " + querystring + " took " + (System.currentTimeMillis() - time));
			
			for(int i=0; i<positions.size(); i++)
				System.out.println("Found string " + querystring + " at position " + positions.elementAt(i));
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}
