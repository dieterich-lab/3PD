/**
 * 
 */
package primerDesign.util;

import java.io.Serializable;
import java.util.Arrays;

/**
 * An implementation of the SimpleContig interface.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SimpleContigImpl implements Serializable, SimpleContig{

	private static final long serialVersionUID = 1L;
	private String id;
	private String comments;
	private char[] sequence;

	/**
	 * Initializes a simple contig.
	 * 
	 * @param id the ID of the contig
	 */
	public SimpleContigImpl(String id) {
		this.id = id;
	}
	
	/**
	 * Initializes a simple contig.
	 * 
	 * @param id the ID of the contig
	 * @param comments the header of the contig
	 */
	public SimpleContigImpl(String id, String comments){
		this.id = id;
		this.comments = comments;
	}
	
	/**
	 * Initializes a simple contig.
	 * 
	 * @param id the ID of the contig
	 * @param sequence the sequence of the contig
	 */
	public SimpleContigImpl(String id, char[] sequence) {
		this.id = id;
		this.sequence =sequence;
	}
	
	/**
	 * Initializes a simple contig.
	 * 
	 * @param id the ID of the contig
	 * @param comments the header of the contig
	 * @param sequence the sequence of the contig
	 */
	public SimpleContigImpl(String id, String comments, char[] sequence){
		this.id = id;
		this.comments = comments;
		this.sequence = sequence;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getName()
	 */
	public String getID(){
		return this.id;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#setName(java.lang.String)
	 */
	public void setID(String name){
		this.id = name;
	}
	
	/*
	 * (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getComments()
	 */
	public String getComments(){
		return this.comments;
	}
	
	/*
	 * (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#setComments(java.lang.String)
	 */
	public void setComments(String comments){
		this.comments = comments;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getSequence()
	 */
	public char[] getSequence(){
		return this.sequence;
	}
		
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#setSequence(char[])
	 */
	public void setSequence(char[] sequence){
		this.sequence = sequence;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getSubsequence(int, int)
	 */
	public String getSubsequence(int start, int end){
		char[] temp = new char[end-start];
		for(int i=start; i<end; i++){
			temp[i-start] = this.sequence[i];
		}
		return new String(temp);
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getSequenceLength()
	 */
	public int getSequenceLength(){
		return this.sequence.length;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#toFastaString()
	 */
	public String toFastaString(){
		return ">" + this.id + "\n" + new String(this.sequence) + "\n";
	}
	
	/**
	 * Returns true iff both SimpleContigs have identical names PLUS identical sequences
	 */
	public boolean equals(Object o){
		SimpleContigImpl other = (SimpleContigImpl) o;
		return this.id.equals(other.getID()) && Arrays.equals(this.sequence,other.sequence);
	}
}
