/**
 * 
 */
package primerDesign.dsc.indexStructures;

import java.io.Serializable;

import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContig;

/**
 * Encapsulates a match position of a query in an index.
 * 
 * The name and the position of the match are supposed to be immutable!
 * 
 * @author Sebastian Fršhler
 *
 */
public class IndexHitImpl implements IndexHit, Serializable {

	private static final long serialVersionUID = 1L;
	private final SimpleContig contig; //the contig the match is located on
	private int position; // the position of the match on contig 'contig'
	private boolean isForwardHit;
	
	/**
	 * Creates a new hit in an index structure.
	 * 
	 * @param contig the contig the hit is located on
	 * @param position the position of the match on contig 'contig'
	 */
	public IndexHitImpl(SimpleContig contig, int position, boolean isForwardHit){
		if(contig == null || position < 0) throw new IllegalArgumentException("Invalid values for String or Position!");
		this.contig = contig;
		this.position = position;
		this.isForwardHit = isForwardHit;
	}
	
	/* (non-Javadoc)
	 * @see esaRmi.IndexHit#getContig()
	 */
	public String getContigName(){
		return this.contig.getID();
	}
	
	/* (non-Javadoc)
	 * @see esaRmi.IndexHit#getPosition()
	 */
	public int getPosition(){
		return this.position;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.IndexHit#getContig()
	 */
	public SimpleContig getContig() {
		return this.contig;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.IndexHit#getContigSequence()
	 */
	public char[] getContigSequence() {
		return this.contig.getSequence();
	}
	
	/*
	 * @see primerDesign.dsc.indexStructures.IndexHit#getSequence(int)
	 */
	public char[] getSequence(int length){
		if(length < 0) throw new IllegalArgumentException();
		char[] result = new char[length];
		char[] temp = this.getContigSequence();
		
		System.arraycopy(temp, this.position, result, 0, length);
		
		if(!this.isForwardHit){
			result = SeqTools.revcompDNA(result).toCharArray();
		}
		return result;
	}
	
	/**
	 * Returns the subsequence of the contig the hit was found on.
	 * 
	 * In case this is a hit on the reverse strand, the reverse complemented automatically.
	 * 
	 * @param start the start of the sequence to be returned (inclusive)
	 * @param end the end of the sequence to be returned (inclusive)
	 * 
	 * @return the subsequence from start to end (both inclusive) of the contig the hit was found on
	 */
	public char[] getSubSequence(int start, int end){
		if(start > end) throw new IllegalArgumentException();
		int length = end - start + 1;
		char[] result = new char[length];
		
		System.arraycopy(this.getContigSequence(), start, result, 0, length);
		
		if(!this.isForwardHit){
			result = SeqTools.revcompDNA(result).toCharArray();
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.IndexHit#setPosition()
	 */
	public void setPosition(int position) {
		if(position < 0) throw new IllegalArgumentException();
		this.position = position;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.IndexHit#getContigLength()
	 */
	public int getContigLength() {
		return this.contig.getSequenceLength();
	}
	
	public boolean equals(Object o){
		IndexHitImpl other = (IndexHitImpl) o;
		return this.position == other.position && this.contig.equals(other.contig) && this.isForwardHit == other.isForwardHit;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.IndexHit#isForwardHit()
	 */
	public boolean isForwardHit() {
		return this.isForwardHit;
	}
	
	public String toString(){
		return this.contig.getID() + " " + this.position + " " + this.isForwardHit;
	}
}
