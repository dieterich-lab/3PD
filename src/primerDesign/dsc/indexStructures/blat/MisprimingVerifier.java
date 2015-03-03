/**
 * 
 */
package primerDesign.dsc.indexStructures.blat;

import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.PrimerSearchParameters;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a primer mispriming verified
 * 
 * @author Sebastian Fršhler
 *
 */
public class MisprimingVerifier {
	/**
	 * Checks for unsafe misprimings of a primer.
	 * 
	 * An unsafe mispriming is an additional priming of a primer that is sufficiently close to a restriction site such that a false-positive amplicon can be generated.
	 * 
	 * @param hits the list of primings for the primer (as found in the index structure)
	 * @param index the restriction sites index
	 * @param params the primer search parameters used by 3PD
	 * 
	 * @return true iff more than one priming exist AND at least one additional priming is sufficiently close to a restriction site
	 */
	public static boolean hasUnsafeMisprimings(ObjectArrayList hits, RestrictionSitesIndex index, PrimerSearchParameters params){
		final int size = hits.size();
		boolean realPrimingFound = false;
		IndexHitImpl hit;
		if(size <= 1) return true; // only the real priming was found
		else{
			for(int i=0; i<size;i++){
				hit = (IndexHitImpl) hits.getQuick(i);
					if(index.isSufficientlyClose(hit.getPosition(), params.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH(), params.getEnzymeName(), params.getTargetOrganism(), hit.getContigName())){
						if(!realPrimingFound) realPrimingFound = true;
						else return true;
					}
			}
			return false;
		}
	}
}
