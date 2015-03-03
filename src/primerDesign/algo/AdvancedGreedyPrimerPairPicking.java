package primerDesign.algo;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;

import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.Constants;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.MyExtendedMath;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/*
 * THIS IS A DEADEND CLASS - NOT TO BE FURTHER IMPLEMENTED!!!
 */
public class AdvancedGreedyPrimerPairPicking{
	private static boolean printTimingStatusInfo = true;
	private static boolean printDebugInfo = true;
	private static boolean runDebugAssertions = true;
	private SimpleTimer timer = new SimpleTimer();
	
	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerPickingAlgorithm#pickBestPrimerSet(primerDesign.dsc.RestrictionSite[])
	 */
	public PrimerPairSet pickBestPrimerSet(RestrictionSite[] optimalSites, PrimerSearchParameters searchParams) {
		if(optimalSites.length < 2) throw new IllegalArgumentException("The list of optimal restriction sites must contain >= 2 elements!");
		
		// enumerate valid primer pairs at each restriction site
		if(printTimingStatusInfo) System.out.print("Enumerating valid primer pairs for each restriction site");
		
		for(int i=0; i<optimalSites.length; i++){
			try{
				if(!optimalSites[i].hasEnumeratedPrimerPairs()){
					optimalSites[i].enumeratePrimerPairs();
				}
				assert(optimalSites[i].hasEnumeratedPrimerPairs() || optimalSites[i].getNumberOfValidPrimerPairs() == 0);
			}
			catch(EmptyResultSetException e){
				// do nothing, just refine sequence region with NO valid primer pairs
			}
			if(optimalSites[i].getNumberOfValidPrimerPairs() < 1){
				PrimerSearch search = searchParams.getPrimerSearch();
				while(optimalSites[i].getNumberOfValidPrimerPairs() < 1){
					optimalSites = search.refineSeqRegion(optimalSites, i, searchParams);
					try{
						optimalSites[i].enumeratePrimerPairs();
					}
					catch(EmptyResultSetException e){
						
					}
					// quit loop when there are no more RSSs to scan (implicit! in refineSeqRegion)
				}
			}
			assert(optimalSites[i].hasEnumeratedPrimerPairs());
			assert(optimalSites[i].getNumberOfValidPrimerPairs() >= 1);
			
			// sort pairs w.r.t. dist to optimal pair
			optimalSites[i].sortPrimerPairs();
			if(runDebugAssertions){
				// assert pairs are properly sorted
				for(int k=1; k<optimalSites[i].getNumberOfValidPrimerPairs(); k++){
					assert(optimalSites[i].getPrimerPair(k-1).getDistanceToOptimalPrimerPair() <= optimalSites[i].getPrimerPair(k).getDistanceToOptimalPrimerPair());
				}
			}
		}

		// sort optimal restriction sites by increasing primer set number, start enumerating pair sets with site with lowest number of primer pairs
		// -> enumerate and check only the minimum number of primer pair sets!
		Arrays.sort(optimalSites, new PrimerSetArrayComparator());
		if(runDebugAssertions){
			// assert sites are sorted in creasing order w.r.t # valid pairs
			for(int j=1; j<optimalSites.length; j++){
				assert(optimalSites[j-1].getNumberOfValidPrimerPairs() <= optimalSites[j].getNumberOfValidPrimerPairs());
			}
		}
		
		if(printTimingStatusInfo) System.out.println(" - done in " + timer.getTimeString());
		
		// pick best set of primer pairs
		if(printTimingStatusInfo) System.out.print("Picking best primer pair");
		
		Vector<PrimerPairSet> threadResults = new Vector<PrimerPairSet>();
		
		int elements = optimalSites[0].getNumberOfValidPrimerPairs();
		
		double increment = 1.0 / Constants.MAX_NUM_PICKING_THREADS;
		int start = 0;
		int end = Math.min(MyExtendedMath.round(elements * increment), elements - 1);
		ObjectArrayList threads = new ObjectArrayList();
		
		for(int i=1; i<=Constants.MAX_NUM_PICKING_THREADS; i++){
			threads.add(new PairSetEvalThread(start, end, optimalSites, threadResults, searchParams));
			if(printDebugInfo) System.err.println("\nInit thread: " + start + "-" + end + "/" + (elements-1));
			start = end +1;
			//assert(i == Constants.MAX_NUM_PICKING_THREADS || start < elements);
			if(start >= elements) break;
			else end = Math.min(MyExtendedMath.round(start + elements * increment), elements-1);
		}
		assert(end == elements - 1);
		
		for(int i=0; i<threads.size(); i++){
			((PairSetEvalThread)threads.get(i)).start();
		}
		
		try {
			for(int i=0; i<threads.size(); i++){
				((PairSetEvalThread)threads.get(i)).join();
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		PrimerPairSet bestPrimerPairs = null;
		PrimerPairSet currentPairSet;
		double currentBestHomogenityScore = Double.MAX_VALUE;
		double currentdOptPrimerPairSet = Double.MAX_VALUE;
		
		int pa_max = Integer.MAX_VALUE;
		int pea_max = Integer.MAX_VALUE;
		
		for(int i=0; i<threadResults.size(); i++){
			currentPairSet = threadResults.get(i);
			
			if(currentPairSet == null) continue;
			
			if(currentPairSet != null && currentPairSet.getHomogenityScore() * searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentPairSet.getAvgDistOptPrimerPair() * searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)currentPairSet.getMaxPairAlignScore())/searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)currentPairSet.getMaxPairAlignEndScore())/searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParams.getPAIR_END_ALIGNMENT_WEIGHT() < currentBestHomogenityScore * searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentdOptPrimerPairSet * searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)pa_max)/searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)pea_max)/searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParams.getPAIR_END_ALIGNMENT_WEIGHT()){
			//if(currentPairSet.getHomogenityScore() < currentBestHomogenityScore && currentPairSet.getMaxPairAlignScore() <= pa_max && currentPairSet.getMaxPairAlignEndScore() <= pea_max && currentPairSet.getAvgDistOptPrimerPair() < currentdOptPrimerPairSet){
				bestPrimerPairs = currentPairSet;
				currentBestHomogenityScore = currentPairSet.getHomogenityScore();
				currentdOptPrimerPairSet = currentPairSet.getAvgDistOptPrimerPair();
				pa_max = currentPairSet.getMaxPairAlignScore();
				pea_max = currentPairSet.getMaxPairAlignEndScore();
			}
		}
		
		if(printTimingStatusInfo) System.out.println(" - done in " + timer.getTimeString());
		assert(bestPrimerPairs == null || bestPrimerPairs.getNumPrimerPairs() == optimalSites.length);
		
		if(bestPrimerPairs == null){
			PrimerSearch search = searchParams.getPrimerSearch();
			assert(search != null);
			Iterator<RestrictionSite> iter = optimalSites[0].getSequenceRegion().getCurrentRestrictionSitesIterator();
			if(iter == null) iter = optimalSites[0].getSequenceRegion().getRestrictionSitesIterator(true);
			
			while(iter.hasNext()){
				if(printTimingStatusInfo) System.out.println("Refine region: 0, using site " + (optimalSites[0].getSequenceRegion().getIteratorPosition() + 1) + "/" + optimalSites[0].getSequenceRegion().getNumRestrictionSites());
				RestrictionSite before = optimalSites[0];
				System.out.print("RSS position: " + before.getPosition() + " -> ");
				optimalSites = search.refineSeqRegion(optimalSites, 0, searchParams);
				System.out.println(optimalSites[0].getPosition());
				assert(!before.equals(optimalSites[0]) && before.getPosition() != optimalSites[0].getPosition());
					
				try{
					bestPrimerPairs = pickBestPrimerSet(optimalSites, searchParams);
					if(bestPrimerPairs != null && bestPrimerPairs.getNumPrimerPairs() == optimalSites.length) break;
				}
				catch(EmptyResultSetException e){
					// do nothing, just discard the current restriction site and test the next one...until all sites have been checked OR a valid primer pair set has been found!
				}
			}
			assert(bestPrimerPairs == null || bestPrimerPairs.getNumPrimerPairs() == optimalSites.length);
			if(bestPrimerPairs == null || bestPrimerPairs.getNumPrimerPairs() < optimalSites.length) throw new EmptyResultSetException("No valid primer pairs can be picked - check your selection parameters/ change enzyme!");
		}
		
		return bestPrimerPairs;
	}
	
	class PairSetEvalThread extends Thread{
		private int start;
		private int end;
		private Vector<PrimerPairSet> result;
		private RestrictionSite[] optimalSites;
		private PrimerSearchParameters searchParams;
		
		private static  final int LOW = 0;  // the next best pair below the current position (list is sorted by distance to optimal primer pair!)
		private static final int HIGH = 1; // the next best pair above the current position (list is sorted by distance to optimal primer pair!)
		
		public PairSetEvalThread(int start, int end, RestrictionSite[] optimalSites, Vector<PrimerPairSet> result, PrimerSearchParameters searchParams){
			this.start = start;
			this.end = end;
			this.optimalSites = new RestrictionSite[optimalSites.length];
			for(int i=0; i<optimalSites.length; i++) this.optimalSites[i] = optimalSites[i];  // copy references because original optimalSites might change between threads during execution time
			this.result = result;
			this.searchParams = searchParams;
		}
		
		public void run(){
			
			System.err.println("Thread " + this.getId() + ": Thread started - range: " + this.start + "-" + this.end);
			PrimerPairSet bestPrimerPairs = null;
			
			double currentBestHomogenityScore = Double.MAX_VALUE;
			double currentdOptPrimerPairSet = Double.MAX_VALUE;
			PrimerPairSet currentBestSet;
			
			int pa_max = Integer.MAX_VALUE;
			int pea_max = Integer.MAX_VALUE;
			
			int[][] bestPos = new int[this.optimalSites.length][2];  // stores the boundary positions of the primer pair which was just added in the list of valid primer pairs foreach optimalSite
			int pos;
			int refineCandidatePos;
			PrimerPair nextBestPair;
			PrimerSearch search = this.searchParams.getPrimerSearch();
			PrimerPairSet currentBestSetBeforeRefinement;
			int[][] bestPosBeforeRefinement = new int[this.optimalSites.length][2];
			
			PrimerPairSetAlignments alignments = new PrimerPairSetAlignments(this.optimalSites.length);
			for(int i=0; i<this.optimalSites.length; i++){
				alignments.addScanRegion(optimalSites[i].getForwardScanSequence(), optimalSites[i].getProbeScanSequence(), this.searchParams);
			}
			
			for(int i=this.start; i<=this.end; i++){
				currentBestSet = new PrimerPairSet(this.optimalSites[0].getPrimerPair(i));
				bestPos[0][LOW] = i;
				bestPos[0][HIGH] = i;
				refineCandidatePos = 0;
				if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Init candidate set with pair: " + i + " having sa: " + this.optimalSites[0].getPrimerPair(i).getMaxPairElementsAlignmentScore().getPairScore() + " and sea: " + this.optimalSites[0].getPrimerPair(i).getMaxPairElementsAlignmentScore().getPairEndScore());
				
				for(int j=1; j<this.optimalSites.length; j++){
					try{
						pos = currentBestSet.addBestScoringPrimerPair(this.optimalSites[j].getPrimerPairs(), alignments, j, this.searchParams, true);
						bestPos[j][LOW] = pos;
						bestPos[j][HIGH] = pos;
						refineCandidatePos = j;
						if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Added best scoring primer pair " + pos + "/" + this.optimalSites[j].getNumberOfValidPrimerPairs() + " at position: " + j + " total scores sa: " + currentBestSet.getMaxPairAlignScore() + " (last pair: " + currentBestSet.getPrimerPair(currentBestSet.getNumPrimerPairs()-1).getMaxPairElementsAlignmentScore().getPairScore() + ") sea: " + currentBestSet.getMaxPairAlignEndScore() + " (last pair: " + currentBestSet.getPrimerPair(currentBestSet.getNumPrimerPairs()-1).getMaxPairElementsAlignmentScore().getPairEndScore() + ")");
					}
					catch(EmptyResultSetException e){
						if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Extending set at position " + j + " failed!");
						currentBestSetBeforeRefinement = new PrimerPairSet(currentBestSet);
						for(int l=0; l<bestPos.length; l++){
							bestPosBeforeRefinement[l][LOW] = bestPos[l][LOW];
							bestPosBeforeRefinement[l][HIGH] = bestPos[l][HIGH];
						}
						// "go to next position that can be refined"
						// while current position to refine is valid (w.r.t index) but has NO more pairs to refine with
						while(refineCandidatePos > 0 && refineCandidatePos < this.optimalSites.length && (bestPos[refineCandidatePos][LOW] < 0 && bestPos[refineCandidatePos][HIGH] >= optimalSites[refineCandidatePos].getNumberOfValidPrimerPairs())){
							if(printDebugInfo) System.err.println("Thread " + this.getId() + ": No more primer pairs to refine position: " + refineCandidatePos + " low: " + bestPos[refineCandidatePos][LOW] + " high: " + bestPos[refineCandidatePos][HIGH] + " (" + optimalSites[refineCandidatePos].getNumberOfValidPrimerPairs() + ")");
							bestPos[refineCandidatePos][LOW] = -1;
							bestPos[refineCandidatePos][HIGH] = -1;
							refineCandidatePos--;
							currentBestSet.deleteLastPrimerPair();
							//alignments.deleteLastRegion();
						}
						// "quit extending primer pair set if there is no position to refine left"
						if(refineCandidatePos <= 0){
							if(printDebugInfo) System.err.println("\t" + "Thread " + this.getId() + ": No valid pairs in first extension step!");
							break; // quit extending this current first primer pair, check next one
						}
						// "else refine position"
						else{
							assert(currentBestSet.getNumPrimerPairs() >= 1);
							if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Refining positon " + refineCandidatePos + ", checking bestPos list - low: " + bestPos[refineCandidatePos][LOW] + " high: " + bestPos[refineCandidatePos][HIGH]);
							// in case no valid best scoring pair can be found, refine list j-1, take 'next-best' pair, rescan list i: if list j has no valid pairs to extend, refine list j-1 and so on...
							nextBestPair = getNextBestPrimerPair(optimalSites, refineCandidatePos, bestPos, currentBestSet, alignments, searchParams);
							while(nextBestPair == null && refineCandidatePos > 1 && refineCandidatePos < this.optimalSites.length){
								if(printDebugInfo) System.err.println("Thread " + this.getId() + ": No valid pair to extend set at position " + refineCandidatePos + ", checking previous position ");
								bestPos[refineCandidatePos][LOW] = -1;
								bestPos[refineCandidatePos][HIGH] = -1;
								refineCandidatePos--;
								currentBestSet.deleteLastPrimerPair();
								//alignments.deleteLastRegion();
								assert(currentBestSet.getNumPrimerPairs() >= 1);
								nextBestPair = getNextBestPrimerPair(optimalSites, refineCandidatePos, bestPos, currentBestSet, alignments, searchParams);
								if(printDebugInfo && nextBestPair == null) System.err.println("Thread " + this.getId() + ": Refining position " + refineCandidatePos + " failed");
							}
							if(nextBestPair != null){
								assert(!currentBestSet.getPrimerPair(refineCandidatePos).toString().equals(nextBestPair.toString()));
								currentBestSet.setPrimerPair(refineCandidatePos, nextBestPair);
								//alignments.addScanRegion(nextBestPair.getForwardPrimer().getRestrictionSite().getForwardScanSequence(), nextBestPair.getHybridizationProbe().getRestrictionSite().getProbeScanSequence(), this.searchParams);
								j = refineCandidatePos;
								if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Refined (substituted) next best pair at position " + refineCandidatePos + " with pair: " + optimalSites[refineCandidatePos].getPrimerPairs().indexOf(nextBestPair, false) + ", continuing at position: " + (j+1));
								
							}else{
								assert(refineCandidatePos >= 1);
								// if all remaining 'next-best' combinations have been scanned, refine sequence region, list j is from, then resume list eval at original position (before entering this namespace)
								RestrictionSite oldSite = optimalSites[j];
								try{
									do{
										optimalSites = search.refineSeqRegion(optimalSites, j, searchParams);
										if(!optimalSites[j].hasEnumeratedPrimerPairs()){
											optimalSites[j].enumeratePrimerPairs();
										}
										assert(optimalSites[j].hasEnumeratedPrimerPairs());
										optimalSites[j].sortPrimerPairs();
										if(runDebugAssertions){
											// assert pairs are properly sorted
											for(int k=1; k<optimalSites[j].getNumberOfValidPrimerPairs(); k++){
												assert(optimalSites[j].getPrimerPair(k-1).getDistanceToOptimalPrimerPair() <= optimalSites[j].getPrimerPair(k).getDistanceToOptimalPrimerPair());
											}
										}
										if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Refining optimal site " + j);
									}while(optimalSites[j].getNumberOfValidPrimerPairs() == 0);
								}
								catch(EmptyResultSetException ex){
									break;
								}
								bestPos[j][LOW] = -1;
								bestPos[j][HIGH] = -1;
								currentBestSet = currentBestSetBeforeRefinement;
								bestPos = bestPosBeforeRefinement;
								refineCandidatePos = currentBestSetBeforeRefinement.getNumPrimerPairs() - 1;
								// update alignment column j and all subsequent (dependent) columns
								for(int k=this.optimalSites.length - 1; k>=j; k--){
									alignments.deleteLastRegion();
								}
								for(int k=j; k<this.optimalSites.length; k++){
									alignments.addScanRegion(this.optimalSites[k].getForwardScanSequence(), this.optimalSites[k].getProbeScanSequence(), searchParams);
								}
								assert(!optimalSites[j].equals(oldSite));
								if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Refined sequence region: " + j + " continuing at position: " + j);
								j--;
								assert(j == refineCandidatePos);
							}
						}
					}
				}
				
				// if not every optimal site has a picked primer pair, skip checking this set, evaluate next candidate set
				if(currentBestSet.getNumPrimerPairs() < this.optimalSites.length){
					if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Rejecting incomplete pair set: " + currentBestSet.getNumPrimerPairs() + "/" + this.optimalSites.length);
					continue;
				}

				if(currentBestSet.getHomogenityScore() < currentBestHomogenityScore || currentBestSet.getAvgDistOptPrimerPair() < currentdOptPrimerPairSet){
				//if(currentBestSet.getHomogenityScore() < currentBestHomogenityScore && currentBestSet.getAvgDistOptPrimerPair() < currentdOptPrimerPairSet){
					currentBestSet.checkInterPairAlignments(alignments); //computeMaxAlignmentScore();
					if(currentBestSet.getMaxPairAlignScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && currentBestSet.getMaxPairAlignEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
						if(currentBestSet.getHomogenityScore() * this.searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentBestSet.getAvgDistOptPrimerPair() * this.searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)currentBestSet.getMaxPairAlignScore())/this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * this.searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)currentBestSet.getMaxPairAlignEndScore())/this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT() < currentBestHomogenityScore * this.searchParams.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentdOptPrimerPairSet * this.searchParams.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)pa_max)/this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * this.searchParams.getPAIR_ALIGNMENT_WEIGHT() + ((double)pea_max)/this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT()){
							bestPrimerPairs = currentBestSet;
							currentBestHomogenityScore = currentBestSet.getHomogenityScore();
							currentdOptPrimerPairSet = currentBestSet.getAvgDistOptPrimerPair();
							pa_max = currentBestSet.getMaxPairAlignScore();
							pea_max = currentBestSet.getMaxPairAlignEndScore();
							if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Accepted best pair with start " + i + " having homogenity: " + currentBestHomogenityScore + " and dOpt: " + currentdOptPrimerPairSet + " pa: " + pa_max + " pea: " + pea_max);
						}
						else if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Reject candidate pair d.t. no improvement in weighted sum");
					}
					else if(printDebugInfo) System.err.println("Thread " + this.getId() + ": Reject d.t. alignment values: PA: " + currentBestSet.getMaxPairAlignScore() + " PEA: " + currentBestSet.getMaxPairAlignEndScore());
				}
			}
			synchronized (this.result) {
				this.result.add(bestPrimerPairs);
			}
		}
		
		private PrimerPair getNextBestPrimerPair(RestrictionSite[] optimalSites, int refineCandidatePos, int[][] bestPos, PrimerPairSet currentBestSet, PrimerPairSetAlignments alignments, PrimerSearchParameters searchParams){
			boolean found = false;
			int low = bestPos[refineCandidatePos][LOW];
			int high = bestPos[refineCandidatePos][HIGH];
			int candidates = optimalSites[refineCandidatePos].getNumberOfValidPrimerPairs();
			PrimerAlignmentScores scores = new PrimerAlignmentScores(0,0);
			PrimerPair result = null;
			PrimerPair candidateLow;
			PrimerPair candidateHigh;
			
			--low;
			++high;
			

			while(!found && low >= 0 && high < candidates){
				candidateLow = optimalSites[refineCandidatePos].getPrimerPair(low);
				candidateHigh = optimalSites[refineCandidatePos].getPrimerPair(high);
				
				if(Math.abs(currentBestSet.getAvgDistOptPrimerPair() - candidateLow.getDistanceToOptimalPrimerPair()) < Math.abs(currentBestSet.getAvgDistOptPrimerPair() - candidateHigh.getDistanceToOptimalPrimerPair())){
					// first check candidate 'low'

					scores = PrimerPairSet.checkInterPairAlignments(candidateLow, currentBestSet, alignments, refineCandidatePos, searchParams);

					if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
						result = candidateLow;
						bestPos[refineCandidatePos][LOW] = low;
						found = true;
					}
					else{
						//scores = PrimerPairSet.checkInterPairAlignments(candidateHigh, currentBestSet, alignments, refineCandidatePos, searchParams);
						if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
							result = candidateHigh;
							bestPos[refineCandidatePos][HIGH] = high;
							found = true;
						}
					}
				}
				else{
					// first check candidate 'high'
					scores = PrimerPairSet.checkInterPairAlignments(candidateHigh, currentBestSet, alignments, refineCandidatePos, searchParams);
					if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
						result = candidateHigh;
						bestPos[refineCandidatePos][HIGH] = high;
						found = true;
					}
					else{
						scores = PrimerPairSet.checkInterPairAlignments(candidateLow, currentBestSet, alignments, refineCandidatePos, searchParams);
						if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
							result = candidateLow;
							bestPos[refineCandidatePos][LOW] = low;
							found = true;
						}
					}
				}
				low--;
				high++;
			}
			while(!found && low < 0 && high < candidates){
				candidateHigh = optimalSites[refineCandidatePos].getPrimerPair(high);
				scores = PrimerPairSet.checkInterPairAlignments(candidateHigh, currentBestSet, alignments, refineCandidatePos, searchParams);
				if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
					result = candidateHigh;
					bestPos[refineCandidatePos][HIGH] = high;
					found = true;
				}
				high++;
			}
			while(!found && low >= 0 && high >= candidates){
				candidateLow = optimalSites[refineCandidatePos].getPrimerPair(low);
				scores = PrimerPairSet.checkInterPairAlignments(candidateLow, currentBestSet, alignments, refineCandidatePos, searchParams);
				if(scores.getPairScore() <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && scores.getPairEndScore() <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
					result = candidateLow;
					bestPos[refineCandidatePos][LOW] = low;
					found = true;
				}
				low--;
			}
			return result;
		}
	}
}
