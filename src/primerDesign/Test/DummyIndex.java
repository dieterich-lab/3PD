/**
 * 
 */
package primerDesign.Test;

import primerDesign.dsc.DNASequenceIndex;

/**
 * Implements a dummy index for testing purposes.
 * 
 * This dummy index always returns one match iff sequence includes scan region and no match else.
 * 
 * @author froehler
 *
 */
public class DummyIndex implements DNASequenceIndex {
	private String sequence;
	private int maxWordSize;
	private boolean includesScanRegion;

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#createIndex(java.lang.String, int, boolean)
	 */
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
		this.sequence = sequence;
		this.maxWordSize = maxWordSize;
		this.includesScanRegion = includesScanRegion;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#getMaxWordSize()
	 */
	public int getMaxWordSize() {
		return this.maxWordSize;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#getSequence()
	 */
	public String getSequence() {
		return this.sequence;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#includesScanRegion()
	 */
	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchMatchPositionsInIndex(java.lang.String)
	 */
	public Integer[] searchMatchPositionsInIndex(String searchString) {
		if(this.includesScanRegion) return new Integer[]{-1};
		else return new Integer[0];
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchNbMatchesInIndex(java.lang.String)
	 */
	public int searchNbMatchesInIndex(String searchString) {
		if(this.includesScanRegion) return 1;
		else return 0;
	}
}
