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
public class PrimerSetArrayComparator implements Comparator{
	public int compare(Object o1, Object o2) {
		RestrictionSite first = (RestrictionSite) o1;
		RestrictionSite second = (RestrictionSite) o2;
		if(first.getNumberOfValidPrimerPairs() < second.getNumberOfValidPrimerPairs()) return -1;
		else if(first.getNumberOfValidPrimerPairs() > second.getNumberOfValidPrimerPairs()) return 1;
		else if(first.getNumberOfValidPrimerPairs() == second.getNumberOfValidPrimerPairs()) return 0;
		else throw new IllegalStateException("Unhandled case");
	}
}
