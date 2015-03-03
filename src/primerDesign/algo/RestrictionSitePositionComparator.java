/**
 * 
 */
package primerDesign.algo;

import java.util.Comparator;

import primerDesign.dsc.RestrictionSite;

/**
 * @author Sebastian Fršhler
 *
 */
public class RestrictionSitePositionComparator implements Comparator {
	public int compare(Object o1, Object o2) {
		RestrictionSite first = (RestrictionSite) o1;
		RestrictionSite second = (RestrictionSite) o2;
		if(first.getPosition() < second.getPosition()) return -1;
		else if(first.getPosition() > second.getPosition()) return 1;
		else if(first.getPosition() == second.getPosition()) return 0;
		else throw new IllegalStateException("Unhandled case");
	}
}
