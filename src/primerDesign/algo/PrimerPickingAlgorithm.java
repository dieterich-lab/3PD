/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.PrimerSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.PrimerSearchParameters;

/**
 * This interface specifies the picking procedure of the best primer pairs from a set of restriction sites with valid 
 * forward, reverse and hybridization probe primers.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface PrimerPickingAlgorithm {
	/**
	 * This method picks a set of most-homogenous primer pairs from a list of restriction sites with associated valid primers.
	 * 
	 * Each restriction site contains valid forward, reverse and hybridization probe primers. The goal is to find the most homogenous
	 * set of valid primer pairs w.r.t. the search parameters. This corresponds to a path throuhg the list of restriction sites from left to right
	 * containing exactly one primer pair per restriction site. The overall set is supposed to be most-homogenous w.r.t the search params
	 * and to be valid.
	 * 
	 * @param optimalSites
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return the optimal primer sets according to the primer search parameters specified
	 */
	public PrimerSet pickBestPrimerSet(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams);
}
