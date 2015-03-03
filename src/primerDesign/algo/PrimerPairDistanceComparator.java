/**
 * 
 */
package primerDesign.algo;

import java.util.Comparator;

import primerDesign.dsc.PrimerPair;

/**
 * Compares two primer pairs, used to order primer pairs w.r.t. increasing distance to optimal primer pair.
 * 
 * @author Sebastian fršhler
 *
 */
public class PrimerPairDistanceComparator implements Comparator{
	
	private static double epsilon = 1E-03;

	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(Object arg0, Object arg1) {
		PrimerPair first = (PrimerPair) arg0;
		PrimerPair second = (PrimerPair) arg1;
		
		if(Math.abs(first.getDistanceToOptimalPrimerPair() - second.getDistanceToOptimalPrimerPair()) < epsilon) return 0;
		else if(first.getDistanceToOptimalPrimerPair() < second.getDistanceToOptimalPrimerPair()) return -1;
		else return 1;
	}
}
