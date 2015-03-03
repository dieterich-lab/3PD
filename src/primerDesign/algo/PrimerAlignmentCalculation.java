/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.PrimerAlignmentScores;

/**
 * Specifies a generic primer alignment.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface PrimerAlignmentCalculation {

	/**
	 * Returns the sa (self alignment) score of a  primer.
	 * 
	 * @return the sa score of a primer
	 */
	public abstract int getSAScore();

	/**
	 * Returns the sea (self end alignment) score of a primer.
	 * 
	 * @return the sea score of a primer
	 */
	public abstract int getSEAScore();

	/**
	 * Returns the pa (pair alignment) score of two primers.
	 * 
	 * @return the pa score of two primers
	 */
	public abstract int getPAScore();

	/**
	 * Returns the pea (pair end alignment) score of two primers.
	 * 
	 * @return the pea score of two primers
	 */
	public abstract int getPEAScore();

	/**
	 * This specifies the 'primer self alignments'.
	 *
	 *@param primer the primer for which the self alignment is to be computed
	 */
	public abstract PrimerAlignmentScores computeSelfAlignment(String primer);

	/**
	 * This specifies the primer pair alignment ('primer1-primer2 alignments').
	 *
	 *@param primer1 the first primer for which the pair alignment is to be computed
	 *@param primer2 the second primer for which the pair alignment is to be computed
	 */
	public abstract PrimerAlignmentScores computePairAlignment(String primer1, String primer2);

}