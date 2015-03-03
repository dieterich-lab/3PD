package primerDesign.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 * Implements some sequence transformation tools that are more efficient than the corresponding biojava tools.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SeqTools {
	
	private static Random random = new Random();
	private static char DEFAULT_MASK_CHAR = 'N';

	/**
	 * Computres and returns a random dna sequence.
	 * 
	 * @param minLength the minimum length of the random dna sequence to be returned
	 * @param maxLength tha maximum length of the random dna sequence to be returned
	 * 
	 * @return a random dna sequence of length [minLength,maxLength]
	 */
	public static String getRandomPrimerSequence(int minLength, int maxLength) {
		int stop = random.nextInt(maxLength-minLength+1);
		
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<=minLength + stop; i++){
			buffer.append(getRandomBase());
		}
		return buffer.toString();
	}
	
	/**
	 * Returns a random DNA base.
	 * 
	 * @return a random DNA base
	 */
	public static String getRandomBase(){
		int base = random.nextInt(4);
		
		String result = null;
		
		switch(base){
			case 0 : result = "A"; break;
			case 1 : result = "T"; break;
			case 2 : result = "G"; break;
			case 3 : result = "C"; break;
		}
		return result;
	}
	
	/**
	 * Returns a random DNA base character.
	 * 
	 * @return a random DNA base character
	 */
	public static char getRandomBaseChar(){
		int base = random.nextInt(4);
		
		char result;
		
		switch(base){
			case 0 : result = 'A'; break;
			case 1 : result = 'T'; break;
			case 2 : result = 'G'; break;
			case 3 : result = 'C'; break;
			default : throw new IllegalStateException("unhaldled case!");
		}
		return result;
	}
	
	/**
	 * Computes the reverse complement of a sequence.
	 * 
	 * @param sequence the sequence for which the reverse complement is to be computed
	 * 
	 * @return the reverse complement of the sequence (the other dna strand in 5'->3' direction)
	 */
	public static String revcompDNA(char[] sequence){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<sequence.length; i++){
			switch(sequence[i]){
			case 'A' : buffer.append('T'); break;
			case 'T' : buffer.append('A'); break;
			case 'G' : buffer.append('C'); break;
			case 'C' : buffer.append('G'); break;
			case 'N' : buffer.append('N'); break;
			case 'U' : buffer.append('A'); break;
			case 'M' : buffer.append('K'); break;
			case 'R' : buffer.append('Y'); break;
			case 'W' : buffer.append('W'); break;
			case 'S' : buffer.append('S'); break;
			case 'Y' : buffer.append('R'); break;
			case 'K' : buffer.append('M'); break;
			case 'V' : buffer.append('B'); break;
			case 'H' : buffer.append('D'); break;
			case 'D' : buffer.append('H'); break;
			case 'B' : buffer.append('V'); break;
			case 'X' : buffer.append('X'); break;

			case 'a' : buffer.append('t'); break;
			case 't' : buffer.append('a'); break;
			case 'g' : buffer.append('c'); break;
			case 'c' : buffer.append('g'); break;
			case 'n' : buffer.append('n'); break;
			case 'u' : buffer.append('a'); break;
			case 'm' : buffer.append('k'); break;
			case 'r' : buffer.append('y'); break;
			case 'w' : buffer.append('w'); break;
			case 's' : buffer.append('s'); break;
			case 'y' : buffer.append('r'); break;
			case 'k' : buffer.append('m'); break;
			case 'v' : buffer.append('b'); break;
			case 'h' : buffer.append('d'); break;
			case 'd' : buffer.append('h'); break;
			case 'b' : buffer.append('v'); break;
			case 'x' : buffer.append('x'); break;
			default : throw new IllegalArgumentException("Unsupported character detected!: " + sequence[i]);
			}
		}
		return buffer.reverse().toString();
	}
	
	/**
	 * Computes the reverse complement of a sequence.
	 * 
	 * @param sequence the sequence for which the reverse complement is to be computed
	 * 
	 * @return the reverse complement of the sequence (the other dna strand in 5'->3' direction)
	 */
	public static byte[] revcompDNAByte(byte[] sequence){
		byte[] result = new byte[sequence.length];
		final int length = sequence.length;
		for(int i=0; i<sequence.length; i++){
			switch(sequence[i]){
			case 'A' : result[length - i - 1] = 'T'; break;
			case 'T' : result[length - i - 1] = 'A'; break;
			case 'G' : result[length - i - 1] = 'C'; break;
			case 'C' : result[length - i - 1] = 'G'; break;
			case 'N' : result[length - i - 1] ='N'; break;
			case 'U' : result[length - i - 1] = 'A'; break;
			case 'M' : result[length - i - 1] = 'K'; break;
			case 'R' : result[length - i - 1] ='Y'; break;
			case 'W' : result[length - i - 1] ='W'; break;
			case 'S' : result[length - i - 1] ='S'; break;
			case 'Y' : result[length - i - 1] ='R'; break;
			case 'K' : result[length - i - 1] ='M'; break;
			case 'V' : result[length - i - 1] ='B'; break;
			case 'H' : result[length - i - 1] ='D'; break;
			case 'D' : result[length - i - 1] ='H'; break;
			case 'B' : result[length - i - 1] ='V'; break;
			case 'X' : result[length - i - 1] ='X'; break;

			case 'a' : result[length - i - 1] = 't'; break;
			case 't' : result[length - i - 1] ='a'; break;
			case 'g' : result[length - i - 1] ='c'; break;
			case 'c' : result[length - i - 1] ='g'; break;
			case 'n' : result[length - i - 1] ='n'; break;
			case 'u' : result[length - i - 1] ='a'; break;
			case 'm' : result[length - i - 1] ='k'; break;
			case 'r' : result[length - i - 1] ='y'; break;
			case 'w' : result[length - i - 1] ='w'; break;
			case 's' : result[length - i - 1] ='s'; break;
			case 'y' : result[length - i - 1] = 'r'; break;
			case 'k' : result[length - i - 1] ='m'; break;
			case 'v' : result[length - i - 1] ='b'; break;
			case 'h' : result[length - i - 1] ='d'; break;
			case 'd' : result[length - i - 1] ='h'; break;
			case 'b' : result[length - i - 1] ='v'; break;
			case 'x' : result[length - i - 1] ='x'; break;
			default : throw new IllegalArgumentException("Unsupported character detected!: " + sequence[i]);
			}
		}
		return result;
	}
	
	/**
	 * Computes the reverse complement of a sequence.
	 * 
	 * @param sequence the sequence for which the reverse complement is to be computed
	 * 
	 * @return the reverse complement of the sequence (the other dna strand in 5'->3' direction)
	 */
	public static String revcompDNA(byte[] sequence){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<sequence.length; i++){
			switch(sequence[i]){
			case 'A' : buffer.append('T'); break;
			case 'T' : buffer.append('A'); break;
			case 'G' : buffer.append('C'); break;
			case 'C' : buffer.append('G'); break;
			case 'N' : buffer.append('N'); break;
			case 'U' : buffer.append('A'); break;
			case 'M' : buffer.append('K'); break;
			case 'R' : buffer.append('Y'); break;
			case 'W' : buffer.append('W'); break;
			case 'S' : buffer.append('S'); break;
			case 'Y' : buffer.append('R'); break;
			case 'K' : buffer.append('M'); break;
			case 'V' : buffer.append('B'); break;
			case 'H' : buffer.append('D'); break;
			case 'D' : buffer.append('H'); break;
			case 'B' : buffer.append('V'); break;
			case 'X' : buffer.append('X'); break;

			case 'a' : buffer.append('t'); break;
			case 't' : buffer.append('a'); break;
			case 'g' : buffer.append('c'); break;
			case 'c' : buffer.append('g'); break;
			case 'n' : buffer.append('n'); break;
			case 'u' : buffer.append('a'); break;
			case 'm' : buffer.append('k'); break;
			case 'r' : buffer.append('y'); break;
			case 'w' : buffer.append('w'); break;
			case 's' : buffer.append('s'); break;
			case 'y' : buffer.append('r'); break;
			case 'k' : buffer.append('m'); break;
			case 'v' : buffer.append('b'); break;
			case 'h' : buffer.append('d'); break;
			case 'd' : buffer.append('h'); break;
			case 'b' : buffer.append('v'); break;
			case 'x' : buffer.append('x'); break;
			default : throw new IllegalArgumentException("Unsupported character detected!: " + sequence[i]);
			}
		}
		return buffer.reverse().toString();
	}
	
	/**
	 * Computes the reverse complement of a DNA base.
	 * 
	 * @param x the DNA base for which the reverse complement is to be computed
	 * 
	 * @return the reverse complement of the DNA base
	 */
	public static char revcompChar(char x){
		switch(x){
			case 'a' : return 't';
			case 't' : return 'a';
			case 'g' : return 'c';
			case 'c' : return 'g';
			case 'n' : return 'n';
			case 'u' : return 'a'; 
			case 'm' : return 'k'; 
			case 'r' : return 'y'; 
			case 'w' : return 'w'; 
			case 's' : return 's'; 
			case 'y' : return 'r'; 
			case 'k' : return 'm'; 
			case 'v' : return 'b'; 
			case 'h' : return 'd'; 
			case 'd' : return 'h'; 
			case 'b' : return 'v'; 
			case 'x' : return 'x'; 
			
			case 'A' : return 'T';
			case 'T' : return 'A';
			case 'G' : return 'C';
			case 'C' : return 'G';
			case 'N' : return 'N';
			case 'U' : return 'A'; 
			case 'M' : return 'K'; 
			case 'R' : return 'Y'; 
			case 'W' : return 'W'; 
			case 'S' : return 'S'; 
			case 'Y' : return 'R'; 
			case 'K' : return 'M'; 
			case 'V' : return 'B'; 
			case 'H' : return 'D'; 
			case 'D' : return 'H'; 
			case 'B' : return 'V'; 
			case 'X' : return 'X'; 
			default : throw new IllegalArgumentException("Unsupported character detected!: " + x);
		}
	}
	
	/**
	 * Computes the complement of a sequence.
	 * 
	 * @param sequence the sequence for which the complement is to be computed
	 * 
	 * @return the complement of the sequence (the other dna strand in 3'->5' direction)
	 */
	public static String complementDNA(char[] sequence){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<sequence.length; i++){
			switch(sequence[i]){
			case 'A' : buffer.append('T'); break;
			case 'T' : buffer.append('A'); break;
			case 'G' : buffer.append('C'); break;
			case 'C' : buffer.append('G'); break;
			case 'N' : buffer.append('N'); break;
			case 'U' : buffer.append('A'); break;
			case 'M' : buffer.append('K'); break;
			case 'R' : buffer.append('Y'); break;
			case 'W' : buffer.append('W'); break;
			case 'S' : buffer.append('S'); break;
			case 'Y' : buffer.append('R'); break;
			case 'K' : buffer.append('M'); break;
			case 'V' : buffer.append('B'); break;
			case 'H' : buffer.append('D'); break;
			case 'D' : buffer.append('H'); break;
			case 'B' : buffer.append('V'); break;
			case 'X' : buffer.append('X'); break;

			case 'a' : buffer.append('t'); break;
			case 't' : buffer.append('a'); break;
			case 'g' : buffer.append('c'); break;
			case 'c' : buffer.append('g'); break;
			case 'n' : buffer.append('n'); break;
			case 'u' : buffer.append('a'); break;
			case 'm' : buffer.append('k'); break;
			case 'r' : buffer.append('y'); break;
			case 'w' : buffer.append('w'); break;
			case 's' : buffer.append('s'); break;
			case 'y' : buffer.append('r'); break;
			case 'k' : buffer.append('m'); break;
			case 'v' : buffer.append('b'); break;
			case 'h' : buffer.append('d'); break;
			case 'd' : buffer.append('h'); break;
			case 'b' : buffer.append('v'); break;
			case 'x' : buffer.append('x'); break;
			default : throw new IllegalArgumentException("Unsupported character detected!: " + sequence[i]);
			}
		}
		return buffer.toString();
	}
	
	/**
	 * Computes the complement of a DNA base.
	 * 
	 * @param x the DNA base for which the complement is to be computed
	 * 
	 * @return the complement of the DNA base
	 */
	public static char complementChar(char x){
		return revcompChar(x);
	}
	
	/**
	 * Computes the reverse of a sequence.
	 * 
	 * @param sequence the sequence for which the reverse is to be computed
	 * @return the reverse of the sequence (the same strand in 3'->5' direction)
	 */
	public static String reverseDNA(String sequence){
		return (new StringBuilder(sequence)).reverse().toString();
	}
	
	/**
	 * Reads in a file in FASTA format.
	 * 
	 * @param fileName the filename to be parsed
	 * 
	 * @return a String of the FASTA sequence
	 * 
	 * @throws IOException
	 */
	public static String readFastaFile(String fileName) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		StringBuffer backgroundSeq = new StringBuffer();
		String line;
		Pattern fastaStart = Pattern.compile(">.*");
		while((line = reader.readLine()) != null){
			if(!fastaStart.matcher(line).matches()){
				//background += line.trim();
				backgroundSeq.append(line.trim());
			}
			else{
				// prohibit direct concatenation of different sequences - pseudo-primers could be found
				backgroundSeq.append(Constants.REPETITIVE_ELEMENT_CHARACTER);
			}
		}
		return backgroundSeq.toString();
	}
	
	/**
	 * Returns all sequence contained in file 'file' concatenated together.
	 * 
	 * @param file the file in FASTA format
	 * 
	 * @return a String containing all sequence contained in file 'file' concatenated together
	 * 
	 * @throws FileNotFoundException
	 * @throws NoSuchElementException
	 * @throws BioException
	 */
	public static String getConcatenatedFastaDNASequence(File file) throws FileNotFoundException, NoSuchElementException, BioException{
		RichSequenceIterator iter = RichSequence.IOTools.readFastaDNA(new BufferedReader(new FileReader(file)), null);
		String sequence = "";
		while(iter.hasNext()){
			if(!sequence.equals("")) sequence += Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase();
			sequence += iter.nextSequence().seqString().toUpperCase();
		}
		return sequence;
	}
	
	/**
	 * Returns the common characters of two sequences (in forward direction!).
	 * 
	 * Comparison is started from the end of both sequences and continued until either sequence is completely traversed.
	 * This method is designed for comparing primer ends to a match region in a dna sequence!
	 * 
	 * @param first the first sequence to be compared
	 * @param second the second sequence to be compared
	 * 
	 * @return the common characters of thw two sequences
	 */
	public static String getCommonChars(char[] first, char[] second){
		StringBuffer buffy = new StringBuffer();
		int firstLength = first.length;
		int secondLength = second.length;
		int minLength = Math.min(firstLength, secondLength);
		for(int i=1; i<=minLength; i++){
			if(first[firstLength - i] == second[secondLength - i]){
				buffy.append(first[firstLength - i]);
			}
		}
		return SeqTools.reverseDNA(buffy.toString());
	}
	
	/**
	 * Masks a given sequence by thresholding the sequence quality values.
	 * 
	 * @param sequence the sequence
	 * @param qualities the quality values of the sequence
	 * @param threshold the threshold for basepair masking
	 * 
	 * @return the masked sequence (in upper case!)
	 */
	public static char[] maskSequence(char[] sequence, byte[] qualities, byte threshold){
		if(sequence.length != qualities.length) throw new IllegalArgumentException("Sequence and quality values array must have same length!");
		if(threshold < 0 || threshold > 97) throw new IllegalArgumentException("Quality threshold must be 0<x<=97!");

		for(int i=0; i<sequence.length; i++){
			if(qualities[i] < threshold) sequence[i] = DEFAULT_MASK_CHAR;
		}
		return sequence;
	}
	
	public static void main(String[] args){
		int instances = Integer.parseInt(args[0]);
		PrimerSearchParameters searchParams = new PrimerSearchParameters();
		for(int i=0; i<instances; i++) System.out.println(SeqTools.getRandomPrimerSequence(searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH()));
	}
}
