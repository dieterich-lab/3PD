/**
 * 
 */
package primerDesign.algo;

import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerPairPickingStatistics;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.Constants;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.MyExtendedMath;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * This class implements a simple greedy primer set picking algorithm.
 * 
 *  This is done by first enumerating each candidate set at each optimal restriction site,
 *  thereafter the sets are ranked by their score per restriction site.
 *  Finally, a 'best path' is greedily searched throught the matrix of candidate primer sets
 *  by starting with one specific primer set at one specific retriction site and then iteratively
 *  adding the primer set with the most-similar score for the remaining restriction sites.
 *  
 *  The best path through this matrix is returned as the best set of primer pairs.
 *  This set of primer pairs has the most-homogenous score-of-single primer pairs of all enumerated sets.
 * 
 * 
 * @author Sebastian Fr�hler
 *
 */
public class SimpleGreedyPrimerPairPicking{ // implements PrimerPickingAlgorithm {
	private boolean printTimingStatusInfo = false;
	private boolean printDebugInfo = true;
	private SimpleTimer timer = new SimpleTimer();
	private boolean doStat = true;
	PrimerPairPickingStatistics stat = null;

	/**
	 * @param printTimingStatusInfo the printTimingStatusInfo to set
	 */
	public void setPrintTimingStatusInfo(boolean printTimingStatusInfo) {
		this.printTimingStatusInfo = printTimingStatusInfo;
	}

	/**
	 * @param printDebugInfo the printDebugInfo to set
	 */
	public void setPrintDebugInfo(boolean printDebugInfo) {
		this.printDebugInfo = printDebugInfo;
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerPickingAlgorithm#pickBestPrimerSet(primerDesign.dsc.RestrictionSite[])
	 */
	public PrimerPairSet pickBestPrimerSet(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams) {
		if(optimalSites.length < 2) throw new IllegalArgumentException("The list of optimal restriction sites must contain >= 2 elements!");
		
		// enumerate valid primer pairs at each restriction site
		if(printTimingStatusInfo) System.out.print("\nEnumerating valid primer pairs for each restriction site");
		
		this.doStat = searchParams.isComputePickingStatistics();
		
		for(int i=0; i<optimalSites.length; i++){
			if(!optimalSites[i].hasEnumeratedPrimerPairs()){
				optimalSites[i].enumeratePrimerPairs();
				//optimalSites[i].sortPrimerPairs();
				if(printDebugInfo) System.err.print("\nsite: " + i + " (valid primer pairs " + optimalSites[i].getNumberOfValidPrimerPairs() + ")");
			}
		}
		if(printDebugInfo) System.err.println();
		for(int i=0; i<optimalSites.length; i++){
			if(optimalSites[i].getNumberOfValidPrimerPairs() == 0){
				PrimerSearch search = searchParams.getPrimerSearch();
			
				while(optimalSites[i].getNumberOfValidPrimerPairs() == 0){
					if(printDebugInfo) System.err.print("Refining sequence region: " + optimalSites[i].getPosition() + " (" + optimalSites[i].getNumberOfValidPrimerPairs() + " valid primer pairs)");
					optimalSites = search.refineSeqRegion(optimalSites, i, searchParams);
					optimalSites[i].enumeratePrimerPairs();
					if(printDebugInfo) System.err.println(" -> " + optimalSites[i].getPosition() + " (" + optimalSites[i].getNumberOfValidPrimerPairs() + " valid primer pairs)");
				}
				if(optimalSites[i].getNumberOfValidPrimerPairs() == 0) throw new IllegalArgumentException("No valid! primer pairs exist for at least one interval: " + i + "!");
			}
		}
		if(Constants.doSortedSitesScreening){
			if(Constants.doGreedySitesScreening){
				// sort optimal restriction sites by increasing primer set number, start enumerating pair sets with site with lowest number of primer pairs
				// -> enumerate and check only the minimum number of primer pair sets!
				Arrays.sort(optimalSites, new PrimerSetArrayComparator());
			}
			else{			
				// sort optimal restriction sites by decreasing primer set number, start enumerating pair sets with site with highest number of primer pairs
				// -> enumerate and check the maximum number of primer pair sets!
				Arrays.sort(optimalSites, new PrimerSetArrayComparatorRev());
			}
		}
		
		if(printTimingStatusInfo) System.out.println(" - done in " + timer.getTimeString());
		
		if(printTimingStatusInfo) System.out.print("Sorting candidates matrix");
		for(int i=0; i<optimalSites.length; i++) optimalSites[i].sortPrimerPairs();
		if(printTimingStatusInfo) System.out.println(" - done in " + timer.getTimeString());
		
		// pick best set of primer pairs
		if(printTimingStatusInfo) System.out.print("Picking best primer pair");
		
//		PrimerPairSet bestPrimerPairs = new PrimerPairSet();
//		int startElements = optimalSites[0].getNumberOfValidPrimerPairs();
//		double currentBestHomogenityScore = Double.MAX_VALUE;
//		PrimerPairSet currentBestSet;
//		
//		int pa_max = Integer.MAX_VALUE;
//		int pea_max = Integer.MAX_VALUE;
//		
//		for(int i=0; i<startElements; i++){
//			currentBestSet = new PrimerPairSet(optimalSites[0].getPrimerPair(i));
//			for(int j=1; j<optimalSites.length; j++){
//				currentBestSet.addBestScoringPrimerPair(optimalSites[j].getPrimerPairs());
//			}
//			if(currentBestSet.getHomogenityScore() < currentBestHomogenityScore){
//				currentBestSet.computeMaxAlignmentScore();
//					if(currentBestSet.getMaxPairAlignScore() <= pa_max && currentBestSet.getMaxPairAlignEndScore() <= pea_max){
//						bestPrimerPairs = currentBestSet;
//						currentBestHomogenityScore = currentBestSet.getHomogenityScore();
//						pa_max = currentBestSet.getMaxPairAlignScore();
//						pea_max = currentBestSet.getMaxPairAlignEndScore();
//					}
//				else if(Constants.PRINT_DEBUG_LOG) System.err.println("PA: " + currentBestSet.getMaxPairAlignScore() + " PEA: " + currentBestSet.getMaxPairAlignEndScore());
//			}
//		}
		Vector<PrimerPairSet> threadResults = new Vector<PrimerPairSet>();
		
		int elements = optimalSites[0].getNumberOfValidPrimerPairs();
		
		double increment = 1.0 / Constants.MAX_NUM_PICKING_THREADS;
		int start = 0;
		int end = Math.min(MyExtendedMath.round(elements * increment), elements - 1);
		ObjectArrayList threads = new ObjectArrayList();
		int[] emptyBestPairCount = new int[optimalSites.length];
		
		for(int i=1; i<=Constants.MAX_NUM_PICKING_THREADS; i++){
			threads.add(new PairSetEvalThread(start, end, optimalSites, threadResults, emptyBestPairCount, this.doStat, searchParams));
			start = end +1;
			if(start >= elements) break;
			assert(i == Constants.MAX_NUM_PICKING_THREADS || start < elements);
			end = Math.min(MyExtendedMath.round(start + elements * increment), elements-1);
		}
		
//		PairSetEvalThread thread1 = new PairSetEvalThread(0, Math.round(elements*0.25f), optimalSites, threadResults);
//		PairSetEvalThread thread2 = new PairSetEvalThread(Math.round(elements*0.25f)+1, Math.round(elements*0.5f), optimalSites, threadResults);
//		PairSetEvalThread thread3 = new PairSetEvalThread(Math.round(elements*0.5f)+1, Math.round(elements*0.75f), optimalSites, threadResults);
//		PairSetEvalThread thread4 = new PairSetEvalThread(Math.round(elements*0.75f)+1, elements-1, optimalSites, threadResults);
		
		for(int i=0; i<threads.size(); i++){
			((PairSetEvalThread)threads.get(i)).start();
		}
//		thread1.start();
//		thread2.start();
//		thread3.start();
//		thread4.start();
		
		try {
			for(int i=0; i<threads.size(); i++){
				((PairSetEvalThread)threads.get(i)).join();
			}
//			thread1.join();
//			thread2.join();
//			thread3.join();
//			thread4.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		PrimerPairSet bestPrimerPairs = null;
		PrimerPairSet currentPairSet;
		double currentBestHomogenityScore = Double.MAX_VALUE;
		double currentdOptPrimerPairSet = Double.MAX_VALUE;
		
//		int pa_max = searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE();
//		int pea_max = searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE();
		int pa_max = Integer.MAX_VALUE;
		int pea_max = Integer.MAX_VALUE;
		
		if(this.doStat){
			this.stat = searchParams.getPickingStat();
			
			for(int i=0; i<threads.size(); i++){
				stat.combineStats(((PairSetEvalThread) threads.get(i)).getStat());
			}
		}
		
		for(int i=0; i<threadResults.size(); i++){
			currentPairSet = threadResults.get(i);
			
			if(currentPairSet == null) continue;
			else if(currentPairSet.getHomogenityScore() * searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentPairSet.getAvgDistOptPrimerPair() * searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)currentPairSet.getMaxPairAlignScore())/searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)currentPairSet.getMaxPairAlignEndScore())/searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParams.getPAIR_END_ALIGNMENT_WEIGHT() < currentBestHomogenityScore * searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentdOptPrimerPairSet * searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)pa_max)/searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)pea_max)/searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParams.getPAIR_END_ALIGNMENT_WEIGHT()){
			//if(currentPairSet.getHomogenityScore() < currentBestHomogenityScore && currentPairSet.getMaxPairAlignScore() <= pa_max && currentPairSet.getMaxPairAlignEndScore() <= pea_max && currentPairSet.getAvgDistOptPrimerPair() < currentdOptPrimerPairSet){
				bestPrimerPairs = currentPairSet;
				currentBestHomogenityScore = currentPairSet.getHomogenityScore();
				currentdOptPrimerPairSet = currentPairSet.getAvgDistOptPrimerPair();
				pa_max = currentPairSet.getMaxPairAlignScore();
				pea_max = currentPairSet.getMaxPairAlignEndScore();
			}
		}

		if(bestPrimerPairs == null){
			
			int refinePos = -1;
			int temp = 0;
			for(int i=0; i<optimalSites.length; i++){
				if(emptyBestPairCount[i] > temp){
					temp = emptyBestPairCount[i];
					refinePos = i;
				}
			}
			if(refinePos < 0) refinePos = new Random().nextInt(optimalSites.length);
			if(refinePos >= 0){
				if(printDebugInfo){ 
					System.err.print("Refining sequence region: " + refinePos + " RSS: " + optimalSites[refinePos].getPosition());
				}
						
				optimalSites = searchParams.getPrimerSearch().refineSeqRegion(optimalSites, refinePos, searchParams);
				if(printDebugInfo){
					System.err.println(" -> " + optimalSites[refinePos].getPosition());
				}
				 
				bestPrimerPairs = pickBestPrimerSet(optimalSites, searchParams);
			}
		}
		
		if(printTimingStatusInfo) System.out.println(" - done in " + timer.getTimeString());
		
		return bestPrimerPairs;
	}
	
	public PrimerPairPickingStatistics getStat(){
		return this.stat;
	}
	
	/**
	 * This class encapsulates a set of tasks to be processed by a thread.
	 * 
	 * From a set of lists of valid primer pairs, a thread is supposed to search through a subset of all possible combinations.
	 * 
	 * @author Sebastian Fr�hler
	 *
	 */
	class PairSetEvalThread extends Thread{
		private int start;
		private int end;
		private Vector<PrimerPairSet> result;
		private RestrictionSite[] optimalSites;
		private int[] emptyBestPairCount;
		private boolean doStat;
		private PrimerPairPickingStatistics stat;
		private PrimerSearchParameters searchParams;
		
		public PairSetEvalThread(int start, int end, RestrictionSite[] optimalSites, Vector<PrimerPairSet> result, int[] emptyBestPairCount, boolean doStat, PrimerSearchParameters searchParams){
			this.start = start;
			this.end = end;
			this.optimalSites = optimalSites;
			this.result = result;
			this.emptyBestPairCount = emptyBestPairCount;
			this.doStat = doStat;
			if(this.doStat) this.stat = new PrimerPairPickingStatistics();
			this.searchParams = searchParams;
		}
		
		public void run(){
			PrimerPairSet bestPrimerPairs = null;
			
			double currentBestHomogenityScore = Double.MAX_VALUE;
			double currentdOptPrimerPairSet = Double.MAX_VALUE;
			PrimerPairSet currentBestSet;
			PrimerPairSet lastBestSet = null;
			int[] emptyBestPairCount = new int[optimalSites.length];
			int temp = -1;
			
			// this is required to prevent init of best pair with pair having undesired alignment values but low dOpt and homogenity
			int pa_max = searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() + 1;
			int pea_max = searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() + 1;
			final int pa = searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() + 1;
			final int pea = searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() + 1;
			
			PrimerPairSetAlignments alignments = new PrimerPairSetAlignments(this.optimalSites.length);
			for(int i=0; i<this.optimalSites.length; i++){
				alignments.addScanRegion(optimalSites[i].getForwardScanSequence(), optimalSites[i].getProbeScanSequence(), this.searchParams);
			}
			
			PrimerPair candidate;
			PrimerPairIterator iter;
			PrimerPair bestPair;
			double minScoreDifference;
			int sa_max;
			int sea_max;
			PrimerAlignmentScores scores;
			
			for(int i=this.start; i<=this.end; i++){
				currentBestSet = new PrimerPairSet(this.optimalSites[0].getPrimerPair(i), this.searchParams); 
//				switch(this.optimalSites[0].getPrimerPair(i).getHashCode()){
//					case -1981834031 : {
//						System.out.println();;
//						break;
//					}
//					case 839705632 : ; break;
//					case 689692681 : ; break;
//					case 678751907 : ; break;
//					case 1709977903 : ; break;
//					case -1180945243 : ; break;
//				}
				if(this.doStat) this.stat.incPairSetsScanned();
				try{
					for(int j=1; j<this.optimalSites.length; j++){
						if(j != currentBestSet.size()) break; // in case no best pair can be added at one position
						temp = j;
//						try{
//							assert(j== currentBestSet.size());
//						}
//						catch(AssertionError e){
//							e.printStackTrace();
//							System.err.println(currentBestSet.getHashCode());
////							System.exit(1);
//						}
						if(Constants.doSortedDOPTScreen){
							// <--- devel test
							//
							assert(currentBestSet.size() == j);
							
							iter = new PrimerPairIterator(this.optimalSites[j].getPrimerPairs(), currentBestSet);
							bestPair = null;
							//int bestPairIndex = -1;
							minScoreDifference = Integer.MAX_VALUE;
							sa_max = -1;
							sea_max = -1;
							scores = null;
							
//							try{
//								assert(j== currentBestSet.size());
//							}
//							catch(AssertionError e){
//								e.printStackTrace();
//								System.err.println(currentBestSet.getHashCode());
////								System.exit(1);
//							}
							
							while(iter.hasNext()){
								candidate = iter.getNext();
								
//								try{
//									assert(j== currentBestSet.size());
//								}
//								catch(AssertionError e){
//									e.printStackTrace();
//									System.err.println(currentBestSet.getHashCode());
////									System.exit(1);
//								}
								
								if(candidate == null) break;
								if(candidate.isNotAcceptablePair()) continue;
								
//								try{
//									assert(j== currentBestSet.size());
//								}
//								catch(AssertionError e){
//									e.printStackTrace();
//									System.err.println(currentBestSet.getHashCode());
////									System.exit(1);
//								}
								
								if(stat != null) stat.incPairsScanned();
								
								if(Math.abs(candidate.getAverageTMPrimers() - currentBestSet.getAverageTMPrimers()) > searchParams.getMAX_PRIMER_TM_DIFFERENCE()){
									if(stat != null)stat.incPairPrimerPrimerTMreject();
									continue;
								}
								else if(stat != null) stat.incPairPrimerPrimerTMaccept();
								
								if(searchParams.isPickTaqManProbe() && Math.abs(candidate.getAverageTMProbe() - currentBestSet.getAverageTMProbes()) > searchParams.getMAX_PRIMER_TM_DIFFERENCE()){ 
									if(stat != null)stat.incPairPrimerProbeTMreject();
									continue;
								}
								else if(stat != null) stat.incPairPrimerProbeTMaccept();
								
								if(candidate.getMaxPairElementsAlignmentScore().getPairScore() > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || candidate.getMaxPairElementsAlignmentScore().getPairEndScore() > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
									if(stat != null) stat.incIntraPairAlignmentReject();
									if(printDebugInfo) System.err.println("Add pair Thread: " + this.getId() + " - invalid intra pair alignment: PA: " + candidate.getMaxPairElementsAlignmentScore().getPairScore() + " PEA: " + candidate.getMaxPairElementsAlignmentScore().getPairEndScore());
									continue;
								}
								
								if(stat != null) stat.incIntraPairAlignmentAccept();
								
								//assert(currentBestSet.getPrimerPairs().size() == j);
								if(searchParams.isCheckInterPairAlignments()){
									scores = PrimerPairSet.checkInterPairAlignments(candidate, currentBestSet, alignments, currentBestSet.size(), searchParams);
									sa_max = scores.getPairScore();
									sea_max = scores.getPairEndScore();
								}
								
								// -> adapted for inter-pair DP-based cross-hybridization check
	//							if(sa_max > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || sea_max > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
	//								continue;
	//							}
								if(sa_max > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || sea_max > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
									if(stat != null) stat.incInterPairAlignmentReject();
									if(printDebugInfo) System.err.println("Add pair Thread: " + this.getId() + "  - invalid inter-pair set alignment: PA: " + sa_max + " PEA: " + sea_max);
									continue;
								}
								
								if(stat != null) stat.incInterPairAlignmentAccept();
								
	//							if(Math.abs(currentPair.getDistanceToOptimalPrimerPair() - getAvgDistOptPrimerPair()) < minScoreDifference && sa_max <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && sea_max <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
	//								bestPair = currentPair;
	//								bestPairIndex = i;
	//								minScoreDifference = Math.abs(currentPair.getDistanceToOptimalPrimerPair() - getAvgDistOptPrimerPair());
	//							}
								assert(!searchParams.isCheckInterPairAlignments() || (sa_max >= 0 && sea_max >= 0 && candidate.getDistanceToOptimalPrimerPair() >= 0));
								if(Math.abs(candidate.getDistanceToOptimalPrimerPair() - currentBestSet.getAvgDistOptPrimerPair()) < minScoreDifference && sa_max <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && sea_max <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
									bestPair = candidate;
									//bestPairIndex = this.optimalSites[j].getPrimerPairs().indexOf(candidate, true);
									minScoreDifference = Math.abs(candidate.getDistanceToOptimalPrimerPair() - currentBestSet.getAvgDistOptPrimerPair());
									if(stat != null) stat.incPairsAccepted();
									if(printDebugInfo) System.err.println("Add pair Thread: " + this.getId() + ": Accept candidate pair " + i + ": sa " + sa_max + " sea " + sea_max + " dOpt " + candidate.getDistanceToOptimalPrimerPair() + " minDiff " + minScoreDifference + " PairSetSize: " + currentBestSet.getNumPrimerPairs() + " new  size: " + (currentBestSet.getNumPrimerPairs()+1) + "/" + searchParams.getNumPrimers());
									
//									try{
//										assert(j== currentBestSet.size());
//									}
//									catch(AssertionError e){
//										e.printStackTrace();
//										System.err.println(currentBestSet.getHashCode());
////										System.exit(1);
//									}
									
									currentBestSet.addPrimerPair(bestPair, sa_max, sea_max);
									
//									try{
//										assert(j + 1 == currentBestSet.size());
//									}
//									catch(AssertionError e){
//										e.printStackTrace();
//										System.err.println(currentBestSet.getHashCode());
////										System.exit(1);
//									}
									
									break;
								}
								else{
									if(stat != null) stat.incPairsRejected();
									if(printDebugInfo) System.err.println("Add pair Thread: " + this.getId() + ": Reject candidate pair " + i + " (NO improvement!): sa " + sa_max + " sea " + sea_max + " dOpt " + candidate.getDistanceToOptimalPrimerPair() + " avgDist " + currentBestSet.getAvgDistOptPrimerPair() + " PairSetSize: " + currentBestSet.getNumPrimerPairs() + "/" + searchParams.getNumPrimers());
								}
							}
	//						if(bestPair == null) throw new EmptyResultSetException();
	//						else{
	//							currentBestSet.addPrimerPair(bestPair);
	//						}
							
						
							// devel test --->
							///
						}
						else{					
							/// OLD
							//currentBestSet.addBestScoringPrimerPair(this.optimalSites[j].getPrimerPairs(), alignments, j, searchParams, true, this.getId(), this.stat);
							/// NEW
							//currentBestSet.addBestScoringPrimerPair(this.optimalSites[j].getPrimerPairs(), alignments, j, searchParams, searchParams.isCheckInterPairAlignments(), this.getId(), this.stat);
							assert(currentBestSet.size() == j);
							currentBestSet.addBestScoringPrimerPair(this.optimalSites[j].getPrimerPairs(), alignments, currentBestSet.size(), searchParams, searchParams.isCheckInterPairAlignments(), this.getId(), this.stat);
						}
					}
				}catch(EmptyResultSetException e){
					
					//
					if(!Constants.doEarlyMMScan){
						boolean restart = false;
						for(int j=0; j<currentBestSet.size(); j++){
							if(j>0 && currentBestSet.getPrimerPair(j).containsUnacceptablePrimers()) restart = true;
						}
						if(restart) i--;
					}
					//
					emptyBestPairCount[temp]++;
					if(this.doStat) this.stat.incEmptyBestPair();
					continue;
				}
				
				if(!Constants.doEarlyMMScan){	
					/*
					 * primers may be included in primer pairs list despite not being tested for misprimings, 
					 * in these cases, a 'late (lazy) mispriming check' needs to be done RIGHT NOW!
					 */
					this.searchParams.getSearchStat().incPrimerIndexCount();
					boolean rescan = false;
					int mismatchAt = -1;
//						System.err.println("Checking primer pairs of size: " + currentBestSet.size());
					for(int j=0; j<currentBestSet.size(); j++){
						if(currentBestSet.getPrimerPair(j).containsUnacceptablePrimers()){
							rescan = true;
							mismatchAt = j;
							synchronized(this.searchParams.getSearchStat()){
								this.searchParams.getSearchStat().incMispriming();
							}
							break;
//								System.err.println("<REJECT " + this.getId() + ">");
//								System.err.println(currentBestSet.getPrimerPair(j).toFormattedString());
//								System.err.println("</REJECT " + this.getId() + ">");
						}else{
							synchronized(this.searchParams.getSearchStat()){
								this.searchParams.getSearchStat().incNoMispriming();
							}
						}
//							else{
//								System.err.println("<ACCEPT " + this.getId() + ">");
//								System.err.println(currentBestSet.getPrimerPair(j).toFormattedString());
//								System.err.println("</ACCEPT " + this.getId() + ">");
//							}
						//lastBestSet = currentBestSet;
					}
					if(rescan && ((lastBestSet != null) ? !lastBestSet.equals(currentBestSet) : true)){
						if(mismatchAt > 0){
							i--;
							lastBestSet = currentBestSet;
//							System.err.println("pair set mispriming prune: " + this.getId() + " i: " + i + "/" + this.end)
							continue;
						}
						else if(mismatchAt == 0){
							lastBestSet = currentBestSet;
							continue;
						}
						else throw new IllegalStateException("Invalid state!");
					}
					else if(lastBestSet != null && lastBestSet.equals(currentBestSet)) break;
				}

				//if(currentBestSet.size() == optimalSites.length) System.err.println("Distance:\t" + currentBestSet.getAvgDistOptPrimerPair() + "\tHomogenity\t" + currentBestSet.getHomogenityScore());
				if(currentBestSet.size() == optimalSites.length && currentBestSet.getHomogenityScore() * this.searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentBestSet.getAvgDistOptPrimerPair() * this.searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)currentBestSet.getMaxPairAlignScore())/this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * this.searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)currentBestSet.getMaxPairAlignEndScore())/this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT() < currentBestHomogenityScore * this.searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentdOptPrimerPairSet * this.searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)pa_max)/this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * this.searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)pea_max)/this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT()){
				//if(currentBestSet.getHomogenityScore() < currentBestHomogenityScore && currentBestSet.getAvgDistOptPrimerPair() < currentdOptPrimerPairSet){
					//currentBestSet.computeMaxAlignmentScore();
						// Inter pair alignments are efficiently checked during incremental pair set extension
						//currentBestSet.checkInterPairAlignments(alignments);
					
					if(stat != null) this.stat.incPairSetsScannedCheck();
					//if(currentBestSet.getMaxPairAlignScore() <= pa_max && currentBestSet.getMaxPairAlignEndScore() <= pea_max && currentBestSet.getMaxPairAlignScore() + currentBestSet.getMaxPairAlignEndScore() < pa_max + pea_max){
					//if((currentBestSet.getMaxPairAlignScore() < pa_max && currentBestSet.getMaxPairAlignEndScore() <= pea_max) || (currentBestSet.getMaxPairAlignScore() <= pa_max && currentBestSet.getMaxPairAlignEndScore() < pea_max)){
					if(currentBestSet.getMaxPairAlignScore() <= pa && currentBestSet.getMaxPairAlignEndScore() <= pea){ // && (currentBestSet.getMaxPairAlignScore() <= pa_max || currentBestSet.getMaxPairAlignEndScore() <= pea_max)){
					//					if(currentBestSet.getMaxPairAlignScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && currentBestSet.getMaxPairAlignEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() &&
//							currentBestSet.getMaxPairAlignScore() <= pa_max && currentBestSet.getMaxPairAlignEndScore() <= pea_max){
						
						bestPrimerPairs = currentBestSet;
						currentBestHomogenityScore = currentBestSet.getHomogenityScore();
						currentdOptPrimerPairSet = currentBestSet.getAvgDistOptPrimerPair();
						pa_max = currentBestSet.getMaxPairAlignScore();
						pea_max = currentBestSet.getMaxPairAlignEndScore();
						if(this.doStat) this.stat.incAcceptedPairSets();
						if(printDebugInfo) System.err.println("Picking Thread " + this.getId() + ": Accepting pair set - PA: " + pa_max + " PEA: " + pea_max);
					}
					else{
						if(this.doStat) this.stat.incNoImprovementOfAlignmentScore();
					}
//					}
					//else if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Rejecting pair set - PA: " + currentBestSet.getMaxPairAlignScore() + " PEA: " + currentBestSet.getMaxPairAlignEndScore());
				}
				else if(currentBestSet.size() < optimalSites.length){
					if(this.doStat) this.stat.incInvalidPairSize();
					if(printDebugInfo) System.err.println("Picking Thread " + this.getId() + ": Reject: Ivalid pair set of size: " + currentBestSet.size());
				}
				else{
					if(this.doStat) this.stat.incNoImprovementInWeightedSum();
					if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Reject: NO improvement of score");					
				}
			}
			synchronized (this.result) {
				this.result.add(bestPrimerPairs);
			}
			synchronized (this.emptyBestPairCount) {
				for(int i=0; i<emptyBestPairCount.length; i++){
					this.emptyBestPairCount[i] += emptyBestPairCount[i];
				}
			}
		}
		
		public PrimerPairPickingStatistics getStat(){
			return this.stat;
		}
	}
}
