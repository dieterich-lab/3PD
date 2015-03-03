/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.RestrictionSite;
import cern.colt.list.ObjectArrayList;

/**
 * @author Sebastian Fršhler
 *
 */
public class RestrictionFragmentSizeFilter {
	/**
	 * Excludes restriction sites with very short three prime fragments.
	 * 
	 * In order to prevent ligation products of fragments A and B of the form: A_xB_y which could
	 * generate multiple amplicons of different lengths, restriction sites with
	 * small three prime fragments are pruned such that false-positive amplicons are required to have at least size > threshold.
	 * 
	 * @param restrictionSites a list of restriction sites to be checked for short fragments
	 * @param minSize the minimum size of an acceptable restriction fragment
	 * 
	 * @return the same list without short fragments
	 */
	public static ObjectArrayList excludeSitesMinLength(ObjectArrayList restrictionSites, int minSize){
		return RestrictionFragmentSizeFilter.excludeSitesMinMaxLength(restrictionSites, minSize, Integer.MAX_VALUE);
	}
	
	/**
	 * Excludes restriction sites with long three prime fragments.
	 * 
	 * In order to prevent very long ligation products, restriction sites with
	 * large three prime fragments are pruned.
	 * 
	 * @param restrictionSites a list of restriction sites to be checked for short fragments
	 * @param maxSize the maximum size of an acceptable restriction fragment
	 * 
	 * @return the same list without long fragments
	 */
	public static ObjectArrayList excludeSitesMaxLength(ObjectArrayList restrictionSites, int maxSize){
		return RestrictionFragmentSizeFilter.excludeSitesMinMaxLength(restrictionSites, 0, maxSize);
	}
	
	/**
	 * Excludes restriction sites with very short and very long three prime fragments.
	 * 
	 * In order to prevent ligation products of fragments A and B of the form: A_xB_y which could
	 * generate multiple amplicons of different lengths, restriction sites with
	 * small three prime fragments are pruned such that false-positive amplicons are required to have at least size > threshold.
	 * 
	 * @param restrictionSites a list of restriction sites to be checked for short fragments
	 * @param minSize the minimum size of an acceptable restriction fragment
	 * 
	 * @return the same list without short fragments
	 */
	public static ObjectArrayList excludeSitesMinMaxLength(ObjectArrayList restrictionSites, int minSize, int maxSize){
		restrictionSites.quickSortFromTo(0, restrictionSites.size()-1, new RestrictionSitePositionComparator());
		
		int lastPosition = -1;
		RestrictionSite currentSite;
		for(int i=0; i<restrictionSites.size(); i++){
			currentSite = (RestrictionSite) restrictionSites.get(i);
			if(lastPosition >= 0 && lastPosition <= currentSite.getPosition() && (currentSite.getPosition() - lastPosition < minSize || currentSite.getPosition() - lastPosition > maxSize)){
				currentSite.getSequenceRegion().removeRestrictionSite(currentSite);
				restrictionSites.remove(restrictionSites.indexOf(currentSite, true));
				i--;
			}
			lastPosition = currentSite.getPosition();
		}
		return restrictionSites;
	}
}
