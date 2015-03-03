/**
 * 
 */
package primerDesign.dsc.indexStructures;

import primerDesign.util.SimpleContig;

/**
 * Specifies the properties of a hit in an index structure.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface IndexHit {

	/**
	 * Returns the name of the contig this hit is located on.
	 * 
	 * @return the name of the contig
	 */
	public abstract String getContigName();
	
	/**
	 * Returns the contig this hit is located on.
	 * 
	 * @return the contig this hit is located on
	 */
	public abstract SimpleContig getContig();
	
	/**
	 * Returns the sequence of the contig this hit is locate on.
	 * 
	 * @return the sequence of the contig this hit is locate on
	 */
	public abstract char[] getContigSequence();
	
	/**
	 * Returns the sequence - starting at hit position - of length 'length'.
	 * 
	 * If this hit is a reverse hit, the sequence of the reverse strand (in 5'->3' direction) - starting at hit position will be returned,
	 * else the sequence of the forward strand (downstream) of hit position will be returned - starting at hit position.
	 * 
	 * @return the sequence - starting at hit position - of length 'length'.
	 */
	public abstract char[] getSequence(int length);
	
	/**
	 * Returns the length of the contig this hit is located on.
	 * 
	 * @return the length of the contig this hit is located on
	 */
	public abstract int getContigLength();

	/**
	 * Returns the position of the hit on contig 'contig' .
	 * 
	 * @return the position of the hit
	 */
	public abstract int getPosition();
	
	/**
	 * Re-sets the position of a hit in a contig.
	 * 
	 * This can be used, e.g., if index hits are used as seeds for further filtering.
	 *
	 * @param position the position of the hit
	 */
	public abstract void setPosition(int position);
	
	/**
	 * Returns whether this hit was found on the forward strand of the dna - if false: found on the reverse strand.
	 *
	 */
	public abstract boolean isForwardHit();
	
	public abstract boolean equals(Object other);
}