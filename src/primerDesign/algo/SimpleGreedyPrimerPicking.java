/**
 * 
 */
package primerDesign.algo;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import primerDesign.dsc.PrimerSet;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;

/**
 * This class implements a simple greedy primer picking strategy for forward primers only.
 * 
 * Starting from an initial primer at one specific restriction site, for the remaining restriction sites
 * the most homogenous primer is iteratively added to the current set.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SimpleGreedyPrimerPicking implements PrimerPickingAlgorithm {
	
	/**
	 * Picks the best primer set - non-exhaustively enumerates all complete paths through the primer matrix.
	 * 
	 * The primer matrix consists of the most-homogenously distributed restriction sites as the columns and
	 * the valid primers for each restriction as the rows for each column.
	 * The best primer set has to contain one "best" primer for each restriction site and therefore is
	 * equivalent to a complete path through the matrix traversing each column.
	 * 
	 * @param optimalSites an array of restriction sites with associated valid primers
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return the best primer set consisting of the most homogenuous primers according to the constraint imposed
	 */
	public PrimerSet pickBestPrimerSet(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams){
		// enumerate and score all possible primer sets and pick best primer set 
		double bestScore = Double.MAX_VALUE;
		PrimerSet bestPrimerSet = null;
		Arrays.sort(optimalSites, new RSSArrayComparator());
		RestrictionSite firstSite = optimalSites[0];
		// create and score a primer set for each valid restriction site of the first optimal restriction site
		// iteratively add the best primer of the subsequent restriction sites
		for(int i=0; i<firstSite.getNumberOfValidUpstreamPrimers(); i++){
			PrimerSet primerSet = new PrimerSet(searchParams);
			primerSet.addForwardPrimer(firstSite.getUpstreamPrimer(i));
			for(int j=1; j < optimalSites.length; j++){
				primerSet.addBestPrimer(optimalSites[j].getValidUpstreamPrimers(), PrimerTypes.forwardPrimer);
			}
			if(primerSet.getScore() < bestScore){
				bestScore = primerSet.getScore();
				bestPrimerSet = primerSet;
			}
		}
		return bestPrimerSet;
	}
	
	/**
	 * Picks the best primer set - using threading.
	 * 
	 * @param optimalSites an array of restriction sites with associated valid primers
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return the best primer set consisting of the most homogenuous primers according to the constraint imposed
	 */
	public PrimerSet pickBestPrimerSetThreaded(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams){		
		Set<Thread> threadNames = new HashSet<Thread>();

		// create and start threads
		for(int i=0; i< Constants.MAX_NUM_PICKING_THREADS; i++){
			PrimerPickingThread thread = new PrimerPickingThread(optimalSites, i, searchParams);
			thread.setName("PrimerPickingThread" + i);
			threadNames.add(thread);
			thread.start();
		}
		// join threads
		Iterator iter = threadNames.iterator();
		while(iter.hasNext()){
			Thread currentThread = (Thread) iter.next();
			try {
				currentThread.join();
			} catch (InterruptedException e) {
				System.err.println(currentThread.getName() + ": Error joining thread!");
				e.printStackTrace();
			}
		}
		
		return PrimerPickingThread.getBestPrimerSet();
	}
	
	/**
	 * Picks the best primer set - using threading.
	 * 
	 * @param optimalSites an array of restriction sites with associated valid primers
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return the best primer set consisting of the most homogenuous primers according to the constraint imposed
	 */
	public PrimerSet pickBestPrimerSetThreaded_test(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams){		
		Set<Thread> threadNames = new HashSet<Thread>();
		Vector<PrimerSet> sets = new Vector<PrimerSet>();

		// create and start threads
		for(int i=0; i< Constants.MAX_NUM_PICKING_THREADS; i++){
			PrimerPickingThread_test thread = new PrimerPickingThread_test(optimalSites, i, searchParams);
			thread.setName("PrimerPickingThread" + i);
			threadNames.add(thread);
			thread.start();
		}
		// join threads
		Iterator iter = threadNames.iterator();
		while(iter.hasNext()){
			Thread currentThread = (Thread) iter.next();
			try {
				currentThread.join();
			} catch (InterruptedException e) {
				System.err.println(currentThread.getName() + ": Error joining thread!");
				e.printStackTrace();
			}
		}
		// retrieve results
		iter = threadNames.iterator();
		while(iter.hasNext()){
			PrimerPickingThread_test currentThread = (PrimerPickingThread_test) iter.next();
			sets.add(currentThread.getBestPrimerSet());
		}
		// sort best primer set from each thread in ascending order (w.r.t score)
		Collections.sort(sets);
		
		return sets.elementAt(0);
	}
}
