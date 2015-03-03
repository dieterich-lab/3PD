/**
 * 
 */
package primerDesign.dsc.indexStructures.primerMisprimingCheck;

import primerDesign.dsc.Primer;
import primerDesign.util.PrimerSearchParameters;

/**
 * Specifies a primer mispriming check used by 3PD.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface PrimerMisprimingCheck{
	/**
	 * Checks whether a given primer has misprimings in a defined background sequence.
	 * 
	 * @param primer the primer sequence to check for misprimings
	 * 
	 * @return true iff there are more than ONE priming in the given background sequence
	 */
	public boolean hasMisprimings(Primer primer);
	
	/**
	 * Sets the 3PD primer search parameters.
	 * 
	 * @param prams the 3PD primer search parameters
	 */
	public void setPrimerSearchParams(PrimerSearchParameters prams);
}
