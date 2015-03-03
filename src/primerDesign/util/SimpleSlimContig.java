/**
 * 
 */
package primerDesign.util;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Implements a simple slim contig.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SimpleSlimContig implements SimpleContig, Serializable{
	private static final long serialVersionUID = 1L;
	private String id;
	private String comments;
	// private char[] sequence;
	private byte[] sequence;

	/**
	 * Initializes a slim contig.
	 * 
	 * @param id the ID of a contig
	 */
	public SimpleSlimContig(String id) {
		this.id = id;
	}
	
	/**
	 * Initializes a slim contig.
	 * 
	 * @param id the ID of a contig
	 * @param sequence the sequence of a contig
	 */
	public SimpleSlimContig(String id, char[] sequence) {
		this.id = id;
		this.sequence = new byte[sequence.length];
		for(int i=0; i<sequence.length; i++) this.sequence[i] = (byte) sequence[i];
		//this.sequence = sequence;
	}
	
	/**
	 * Initializes a slim contig.
	 * 
	 * @param id the ID of a contig
	 * @param comments the header of a contig
	 * @param sequence the sequence of a contig
	 */
	public SimpleSlimContig(String id, String comments, char[] sequence) {
		this.id = id;
		this.comments = comments;
		this.sequence = new byte[sequence.length];
		for(int i=0; i<sequence.length; i++) this.sequence[i] = (byte) sequence[i];
		//this.sequence = sequence;
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
	public void setID(String id){
		this.id = id;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getComments()
	 */
	@Override
	public String getComments() {
		return this.comments;
	}

	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#setComments(java.lang.String)
	 */
	@Override
	public void setComments(String comments) {
		this.comments = comments;		
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getSequence()
	 */
	public char[] getSequence(){
		char[] result = new char[this.sequence.length];
		for(int i=0; i<this.sequence.length; i++) result[i] = (char) this.sequence[i];
		return result;
		//return this.sequence;
	}
	
	/**
	 * Returns the sequence of a slim contig.
	 * 
	 * @return the sequence of a slim contig
	 */
	public byte[] getSlimSequence(){
		return this.sequence;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#setSequence(char[])
	 */
	public void setSequence(char[] sequence){
		this.sequence = new byte[sequence.length];
		for(int i=0; i<sequence.length; i++) this.sequence[i] = (byte) sequence[i];
		//this.sequence = sequence;
	}
	
	/**
	 * Sets the sequence of a slim contig.
	 * 
	 * @param sequence the sequence to be set 
	 */
	public void setSlimSequence(byte[] sequence){
		this.sequence = sequence;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.util.SimpleContig#getSubsequence(int, int)
	 */
	public String getSubsequence(int start, int end){
		char[] temp = new char[end-start];
		for(int i=start; i<end; i++){
			temp[i-start] = (char) this.sequence[i];
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
		SimpleSlimContig other = (SimpleSlimContig) o;
		return this.id.equals(other.getID()) && Arrays.equals(this.sequence,other.sequence);
	}
}
