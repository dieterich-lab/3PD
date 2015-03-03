package primerDesign.algo;

import primerDesign.dsc.PrimerSet;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;

/**
 * Implements a threaded primer picking.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPickingThread extends Thread {

	private RestrictionSite[] optimalSites;
	private static PrimerSet bestPrimerSet;
	private static int numParallelThreads = Constants.MAX_NUM_PICKING_THREADS;
	private static double bestScore = Double.MAX_VALUE;
	private int start;
	private PrimerSearchParameters searchParams;
	
	public PrimerPickingThread(RestrictionSite[] optimalsites, int start, PrimerSearchParameters searchParams){
		this.optimalSites = optimalsites;
		this.start = start;
		this.searchParams = searchParams;
		this.bestPrimerSet = new PrimerSet(searchParams);
	}
	
	public void run() {
		RestrictionSite firstSite = this.optimalSites[0];
		for (int i = this.start; i < firstSite.getNumberOfValidUpstreamPrimers(); i+=PrimerPickingThread.numParallelThreads) {
			// create primer set
			PrimerSet primerSet = new PrimerSet(this.searchParams);
			primerSet.addForwardPrimer(firstSite.getUpstreamPrimer(i));
			// add "best" primer from each optimal restriction site
			for(int j=1; j < optimalSites.length; j++){
				primerSet.addBestPrimer(optimalSites[j].getValidUpstreamPrimers(), PrimerTypes.forwardPrimer);
			}
			// score primer set: store it in result set
			try{
				update(primerSet);
			}catch(Exception e){
				e.printStackTrace();
				System.err.println(Thread.currentThread().getName() + ": Error writing to result set!");
			}
		}
	}
	
	private synchronized void update(PrimerSet primerSet){
		if(primerSet.getScore() < PrimerPickingThread.bestScore){
			PrimerPickingThread.bestPrimerSet = primerSet;
			PrimerPickingThread.bestScore = primerSet.getScore();
		}
	}
	
	public static PrimerSet getBestPrimerSet(){
		return PrimerPickingThread.bestPrimerSet;
	}

}
