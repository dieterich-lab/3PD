/**
 * 
 */
package primerDesign.util;

/**
 * Implements a simple contig consisting of name, comments (optional) and sequence (optional).
 * 
 * Implementations should conform to http://www.ncbi.nlm.nih.gov/blast/fasta.shtml.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface SimpleContig {

	/**
	 * Returns the contig ID.
	 * 
	 * @return the ID of the contig
	 */
	public abstract String getID();

	/**
	 * Sets a contig's ID.
	 * 
	 * @param id the contig ID to be set
	 */
	public abstract void setID(String id);
	
	/**
	 * Returns the header of a FASTA sequence.
	 * 
	 * @return the header of a FASTA sequence
	 */
	public abstract String getComments();
	
	/**
	 * Sets the header of a FASTA sequence.
	 * 
	 * @param comments the header to be set
	 */
	public abstract void setComments(String comments);

	/**
	 * Returns a char[] representation of the sequenc ein FASTA format.
	 * 
	 * @return a char[] representation of the sequenc ein FASTA format
	 */
	public abstract char[] getSequence();

	/**
	 * Sets the FASTA sequence of a contig.
	 * 
	 * @param sequence the FASTA sequence
	 */
	public abstract void setSequence(char[] sequence);

	/**
	 * Returns the substring of the sequence of this contig starting at 'start' and extending until 'end-1'.
	 * 
	 * The length of the returned subsequence is end-start!
	 * 
	 * @param start the start - inclusive
	 * @param end the end - exclusive
	 * 
	 * @return a subsequence ranging from start to end-1
	 */
	public abstract String getSubsequence(int start, int end);

	/**
	 * Returns the length of a sequence.
	 * 
	 * @return the length of a sequence
	 */
	public abstract int getSequenceLength();

	/**
	 * Returns a String representation of the sequenc ein FASTA format.
	 * 
	 * @return a String representation of the sequenc ein FASTA format
	 */
	public abstract String toFastaString();
}