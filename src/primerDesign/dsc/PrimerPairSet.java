/**
 * 
 */
package primerDesign.dsc;

import primerDesign.algo.PrimerPairSetAlignments;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.PrimerSearchParameters;
import cern.colt.list.ObjectArrayList;

/**
 * This class encapsulates a set of primer pairs.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPairSet{
	private ObjectArrayList primerPairs;
	private double sumDistanceToOptimalPrimerPair;  // the sum of the scores of the individual primer pairs
	private double homogenityScore; // the score representing the homogenity of primers in this primer pair set -> the StdDev of the score
	private double averageTMPrimers;
	private double averageTMProbes;
	private boolean wasChanged;
	private int pa_max;
	private int pea_max;
	private PrimerSearchParameters searchParams;
	
	private static boolean printDebugInfo = false;
	
	/***
	 * Initializes a new empty primer pair set.
	 * 
	 * @param searchParams the 3PD search prameters
	 */
	public PrimerPairSet(PrimerSearchParameters searchParams){
		this.primerPairs = new ObjectArrayList();
		this.sumDistanceToOptimalPrimerPair = -1;
		this.averageTMPrimers = -1;
		this.averageTMProbes = -1;
		this.wasChanged = true;
		this.pa_max = -1;
		this.pea_max = -1;
		this.searchParams = searchParams;
	}
	
	/**
	 * Initializes a new primer pair set, starting with primer pair 'pair'.
	 * 
	 * @param pair the primer pair to be the initial element of this primer pair set
	 * @param searchParams the 3PD search parameters
	 */
	public PrimerPairSet(PrimerPair pair, PrimerSearchParameters searchParams){
		this.primerPairs = new ObjectArrayList();
		this.primerPairs.add(pair);

		double score = pair.getDistanceToOptimalPrimerPair();
		this.sumDistanceToOptimalPrimerPair = score;
		this.averageTMPrimers = pair.getAverageTMPrimers();
		this.averageTMProbes = pair.getAverageTMProbe();
		
		this.wasChanged = true;
		this.pa_max = -1;
		this.pea_max = -1;
		this.searchParams = searchParams;
	}
	
	/**
	 * Initializes a primer pair set.
	 * 
	 * @param pairSet the pair set to initialize with
	 * @param searchParams the 3PD search parameters
	 */
	public PrimerPairSet(PrimerPairSet pairSet, PrimerSearchParameters searchParams){
		this.primerPairs = new ObjectArrayList();
		for(int i=0; i<pairSet.getNumPrimerPairs(); i++){
			this.primerPairs.add(pairSet.getPrimerPair(i));
		}
		
		this.sumDistanceToOptimalPrimerPair = pairSet.getAvgDistOptPrimerPair() * pairSet.getNumPrimerPairs();
		this.homogenityScore = pairSet.getHomogenityScore();
		this.averageTMPrimers = pairSet.averageTMPrimers;
		this.averageTMProbes = pairSet.getAverageTMProbes();
		this.pa_max = pairSet.getMaxPairAlignScore();
		this.pea_max = pairSet.getMaxPairAlignEndScore();
		this.wasChanged = false;
		this.searchParams = searchParams;
	}
	
	/**
	 * Adds the best scoring primer pair to the current primer pair set.
	 * 
	 * The best scoring primer pair is the pair with the closest (most-similar) score w.r.t the current primer pair set
	 * 
	 * @param pairs
	 * @param alignments the alignments
	 * @param position the position
	 * @param searchParams the 3PD search parameters
	 * @param checkInterPairAlignments whether to check inter pair alignments as well
	 * 
	 * @return the position of the best primer pair which was added in the list of enumerated valid primer pairs
	 */
	public int addBestScoringPrimerPair(ObjectArrayList pairs, PrimerPairSetAlignments alignments, int position, PrimerSearchParameters searchParams, boolean checkInterPairAlignments, long threadID, PrimerPairPickingStatistics stat){
		if(pairs.size() == 0) throw new IllegalArgumentException("Empty primer pair list detected!");
		
		PrimerPair bestPair = null;
		int bestPairIndex = -1;
		PrimerPair currentPair;
		double minScoreDifference = Integer.MAX_VALUE;
		int sa_max = -1;
		int sea_max = -1;
		PrimerAlignmentScores scores;
		
		for(int i=0; i<pairs.size(); i++){
			currentPair = (PrimerPair) pairs.getQuick(i);
			if(currentPair.isNotAcceptablePair()) continue;
			
			if(stat != null) stat.incPairsScanned();
			assert(currentPair.getMaxPairElementsAlignmentScore().getPairScore() >= 0 && currentPair.getMaxPairElementsAlignmentScore().getPairEndScore() >= 0);
			
			//Test for primer pair sets TM compatibility (primers and probes separated)
			if(Math.abs(currentPair.getAverageTMPrimers() - this.getAverageTMPrimers()) > searchParams.getMAX_PRIMER_TM_DIFFERENCE()){
				stat.incPairPrimerPrimerTMreject();
				continue;
			}
			else stat.incPairPrimerPrimerTMaccept();
			if(searchParams.isPickTaqManProbe() && Math.abs(currentPair.getAverageTMProbe() - this.getAverageTMProbes()) > searchParams.getMAX_PRIMER_TM_DIFFERENCE()){ 
				stat.incPairPrimerProbeTMreject();
				continue;
			}
			else stat.incPairPrimerProbeTMaccept();
			
			if(currentPair.getMaxPairElementsAlignmentScore().getPairScore() > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || currentPair.getMaxPairElementsAlignmentScore().getPairEndScore() > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
				if(stat != null) stat.incIntraPairAlignmentReject();
				if(printDebugInfo) System.err.println("Add pair Thread: " + threadID + " - invalid intra pair alignment: PA: " + currentPair.getMaxPairElementsAlignmentScore().getPairScore() + " PEA: " + currentPair.getMaxPairElementsAlignmentScore().getPairEndScore());
				continue;
			}
			if(stat != null) stat.incIntraPairAlignmentAccept();
			
			assert(this.primerPairs.size() == position);
			if(checkInterPairAlignments){
				//scores = checkInterPairAlignments(currentPair, this, alignments, position, searchParams);
				scores = checkInterPairAlignments(currentPair, this, alignments, this.primerPairs.size(), searchParams);
				sa_max = scores.getPairScore();
				sea_max = scores.getPairEndScore();
			}
			
			// -> adapted for inter-pair DP-based cross-hybridization check
//			if(sa_max > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || sea_max > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
//				continue;
//			}
			if(sa_max > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || sea_max > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
				if(stat != null) stat.incInterPairAlignmentReject();
				if(printDebugInfo) System.err.println("Add pair Thread: " + threadID + "  - invalid inter-pair set alignment: PA: " + sa_max + " PEA: " + sea_max);
				continue;
			}
			if(stat != null) stat.incInterPairAlignmentAccept();
			
//			if(Math.abs(currentPair.getDistanceToOptimalPrimerPair() - getAvgDistOptPrimerPair()) < minScoreDifference && sa_max <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && sea_max <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
//				bestPair = currentPair;
//				bestPairIndex = i;
//				minScoreDifference = Math.abs(currentPair.getDistanceToOptimalPrimerPair() - getAvgDistOptPrimerPair());
//			}
			assert(!checkInterPairAlignments || (sa_max >= 0 && sea_max >= 0 && currentPair.getDistanceToOptimalPrimerPair() >= 0));
			if(Math.abs(currentPair.getDistanceToOptimalPrimerPair() - this.getAvgDistOptPrimerPair()) < minScoreDifference && sa_max <= searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && sea_max <= searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
				bestPair = currentPair; 
				bestPairIndex = i;
				minScoreDifference = Math.abs(currentPair.getDistanceToOptimalPrimerPair() - this.getAvgDistOptPrimerPair());
				if(stat != null) stat.incPairsAccepted();
				if(printDebugInfo) System.err.println("Add pair Thread: " + threadID + ": Accept candidate pair " + i + ": sa " + sa_max + " sea " + sea_max + " dOpt " + currentPair.getDistanceToOptimalPrimerPair() + " minDiff " + minScoreDifference + " PairSetSize: " + this.primerPairs.size() + " new  size: " + (this.primerPairs.size()+1) + "/" + searchParams.getNumPrimers());
			}
			else{
				if(stat != null) stat.incPairsRejected();
				if(printDebugInfo) System.err.println("Add pair Thread: " + threadID + ": Reject candidate pair " + i + " (NO improvement!): sa " + sa_max + " sea " + sea_max + " dOpt " + currentPair.getDistanceToOptimalPrimerPair() + " avgDist " + this.getAvgDistOptPrimerPair() + " PairSetSize: " + this.primerPairs.size() + "/" + searchParams.getNumPrimers());
			}
		}
		if(bestPair == null) throw new EmptyResultSetException();
		this.primerPairs.add(bestPair);
		this.averageTMPrimers += bestPair.getAverageTMPrimers();
		this.averageTMProbes += bestPair.getAverageTMProbe();
		if(this.pa_max < sa_max) this.pa_max = sa_max;
		if(this.pea_max < sea_max) this.pea_max = sea_max;
		this.sumDistanceToOptimalPrimerPair += bestPair.getDistanceToOptimalPrimerPair();
		this.wasChanged = true;
		
		if(printDebugInfo){
			System.err.println("Add pair Thread: " + threadID + " Adding pair: " + bestPairIndex + " new pairSetSize: " + this.primerPairs.size() + "/" + searchParams.getNumPrimers());
			if(this.primerPairs.size() == searchParams.getNumPrimers()) System.err.println(this.toString());
		}
		
		return bestPairIndex;
	}
	
	/**
	 * Adds another primer pair to the set.
	 * 
	 * @param pair another primer pair to add
	 * @param saMax the maximum self-alignment value of 'pair'
	 * @param seaMax the maximum self-end alignment value of 'pair'
	 */
	public void addPrimerPair(PrimerPair pair, int saMax, int seaMax){
		this.primerPairs.add(pair);
		this.averageTMPrimers += pair.getAverageTMPrimers();
		this.averageTMProbes += pair.getAverageTMProbe();
		if(this.pa_max < saMax) this.pa_max = saMax;
		if(this.pea_max < seaMax) this.pea_max = seaMax;
		this.sumDistanceToOptimalPrimerPair += pair.getDistanceToOptimalPrimerPair();
		this.wasChanged = true;
	}
	
	/**
	 * Checks the pairwise alignments of a primer pair 'pair' to a list of primer pairs 'pairSet'.
	 * 
	 * The evaluation stops in case one alignment value exceed the pair thresholds as specified in 'searchParams'.
	 * Therefore this is just a screen if all all (pairwise and intra-pair) alignment values are <= the threshold specified
	 * and NO exhaustive alignment between each pair (in this specific case!).
	 * 
	 * @param pair the pair to pairwise align against a list of pairs
	 * @param pairSet the list of primer pair to align the pair against
	 * @param searchParams the parameters containing max thresholds for alignment values
	 * 
	 * @return the scores of the pairwise alignment OR the first score exceeding the specified threshold!
	 */
	public static PrimerAlignmentScores checkInterPairAlignments(PrimerPair pair, PrimerPairSet pairSet, PrimerPairSetAlignments alignments, int position, PrimerSearchParameters searchParams){
//		PrimerPair other;
//		PrimerAlignmentScores scores;
//		int sa_max = -1;
//		int sea_max = -1;
//		// compute pairwise alignment to all primer pairs already included in set (only: FW1-FW2, FW1-Probe2, FW2-Probe1)!
//		for(int j=0; j<pairSet.getNumPrimerPairs(); j++){
//			other = pairSet.getPrimerPair(j);
//			scores = other.getAlignmentValues(pair);
//			
//			// check inter-pair alignment scores
//			if(scores.getPairScore() > sa_max){
//				sa_max = scores.getPairScore();
//			}
//			if(scores.getPairEndScore() > sea_max){
//				sea_max = scores.getPairEndScore();
//			}
//			
//			// check intra-pair alignment scores of already included pairs
//			if(other.getMaxPairElementsAlignmentScore().getPairScore() > sa_max){
//				sa_max = other.getMaxPairElementsAlignmentScore().getPairScore();
//			}
//			if(other.getMaxPairElementsAlignmentScore().getPairEndScore() > sea_max){
//				sea_max = other.getMaxPairElementsAlignmentScore().getPairEndScore();
//			}
//			if(sa_max > searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || sea_max > searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
//				break;
//			}
//		}
//		if(pair.getMaxPairElementsAlignmentScore().getPairScore() > sa_max){
//			sa_max = pair.getMaxPairElementsAlignmentScore().getPairScore();
//		}
//		if(pair.getMaxPairElementsAlignmentScore().getPairEndScore() > sea_max){
//			sea_max = pair.getMaxPairElementsAlignmentScore().getPairEndScore();
//		}
//		
//		return new PrimerAlignmentScores(sa_max, sea_max);
		
		assert(position <= pairSet.getNumPrimerPairs());
		
		PrimerPair other = null;
		Primer primer1;
		Primer primer2;
		int fw1Fw2 = -1;
		int fw1Fw2Max = -1;
		int fw1Fw2MaxEnd = -1;
		int fw1P2 = -1;
		int fw1P2Max = -1;
		int fw2P1 = -1;
		int fw2P1Max = -1;
		PrimerAlignmentScores scores = null;
		for(int i=0; i<position; i++){

			other = pairSet.getPrimerPair(i);

			primer1 = other.getForwardPrimer();
			primer2 = pair.getForwardPrimer();

			scores = alignments.getAlignment(i, primer1.getPositionInScanSequence(), primer1.getLength(), position, primer2.getPositionInScanSequence(), primer2.getLength(), AlignmentType.forward1Forward2);
			fw1Fw2 = scores.getPairScore();

			if(searchParams.isPickTaqManProbe()){
				primer2 = pair.getHybridizationProbe();
				fw1P2 = alignments.getAlignment(i, primer1.getPositionInScanSequence(), primer1.getLength(), position, primer2.getPositionInScanSequence(), primer2.getLength(), AlignmentType.forward1Probe2).getPairScore();
				
				primer1 = other.getHybridizationProbe();
				primer2 = pair.getForwardPrimer();
				fw2P1 = alignments.getAlignment(i, primer2.getPositionInScanSequence(), primer2.getLength(), position, primer1.getPositionInScanSequence(), primer1.getLength(), AlignmentType.forward2Probe1).getPairScore();
			}
			if(fw1Fw2 > fw1Fw2Max) fw1Fw2Max = fw1Fw2;
			if(scores.getPairEndScore() > fw1Fw2MaxEnd) fw1Fw2MaxEnd = scores.getPairEndScore();
			if(fw1P2 > fw1P2Max) fw1P2Max = fw1P2;
			if(fw2P1 > fw2P1Max) fw2P1Max = fw2P1;
		}
		
		return new PrimerAlignmentScores(Math.max(fw1Fw2Max, Math.max(fw1P2Max, fw2P1Max)), fw1Fw2MaxEnd);
		//return new PrimerAlignmentScores(0,0);
	}
	
	/**
	 * Checks inter pair alignments.
	 * 
	 * @param alignments the pair alignments
	 */
	public void checkInterPairAlignments(PrimerPairSetAlignments alignments){
		PrimerAlignmentScores scores;
		PrimerPair pair1;
		PrimerPair pair2;
		Primer primer1;
		Primer primer2;
		int fw1Fw2 = -1;
		int fw1Fw2Max = -1;
		int fw1Fw2MaxEnd = -1;
		int fw1P2 = -1;
		int fw1P2Max = -1;
		int fw2P1 = -1;
		int fw2P1Max = -1;
		
		for(int i=0; i<this.primerPairs.size(); i++){
			pair1 = (PrimerPair) this.primerPairs.getQuick(i);
			for(int j=i+1; j<this.primerPairs.size(); j++){
				pair2 = (PrimerPair) this.primerPairs.getQuick(j);
				
				primer1 = pair1.getForwardPrimer();
				primer2 = pair2.getForwardPrimer();

				scores = alignments.getAlignment(i, primer1.getPositionInScanSequence(), primer1.getLength(), j, primer2.getPositionInScanSequence(), primer2.getLength(), AlignmentType.forward1Forward2);
				fw1Fw2 = scores.getPairScore();
				
				
				primer2 = pair2.getHybridizationProbe();
				fw1P2 = alignments.getAlignment(i, primer1.getPositionInScanSequence(), primer1.getLength(), j, primer2.getPositionInScanSequence(), primer2.getLength(), AlignmentType.forward1Probe2).getPairScore();
				
				primer1 = pair1.getHybridizationProbe();
				primer2 = pair2.getForwardPrimer();
				fw2P1 = alignments.getAlignment(i, primer2.getPositionInScanSequence(), primer2.getLength(), j, primer1.getPositionInScanSequence(), primer1.getLength(), AlignmentType.forward2Probe1).getPairScore();

				if(fw1Fw2 > fw1Fw2Max) fw1Fw2Max = fw1Fw2;
				if(scores.getPairEndScore() > fw1Fw2MaxEnd) fw1Fw2MaxEnd = scores.getPairEndScore();
				if(fw1P2 > fw1P2Max) fw1P2Max = fw1P2;
				if(fw2P1 > fw2P1Max) fw2P1Max = fw2P1;
			}
		}
		this.pa_max = Math.max(fw1Fw2Max, Math.max(fw1P2Max, fw2P1Max));
		this.pea_max = fw1Fw2MaxEnd;
	}
	
	/**
	 * Sets a primer pair.
	 * 
	 * @param index the position of the pair within this set
	 * @param primerPair the primer pair to set
	 */
	public void setPrimerPair(int index, PrimerPair primerPair){
		if(index < this.primerPairs.size()){
			this.averageTMPrimers -= ((PrimerPair)this.primerPairs.getQuick(index)).getAverageTMPrimers();
			this.averageTMProbes -= ((PrimerPair)this.primerPairs.getQuick(index)).getAverageTMProbe();
		}
		this.averageTMPrimers += primerPair.getAverageTMPrimers();
		this.averageTMProbes += primerPair.getAverageTMProbe();
		this.primerPairs.set(index, primerPair);
	}
	
	/**
	 * Deletes the last primer pair of this set.
	 * 
	 */
	public void deleteLastPrimerPair(){
		this.averageTMPrimers -= ((PrimerPair)this.primerPairs.getQuick(this.primerPairs.size() - 1)).getAverageTMPrimers();
		this.averageTMProbes -= ((PrimerPair)this.primerPairs.getQuick(this.primerPairs.size() - 1)).getAverageTMProbe();
		this.primerPairs.remove(this.primerPairs.size() - 1);
	}
	
	/**
	 * Returns all primer pairs of this set.
	 * 
	 * @return a reference to all primer pairs of this set
	 */
	public ObjectArrayList getPrimerPairs(){
		return this.primerPairs;
	}
	
	/**
	 * Returns the number of primer pairs in this set.
	 * 
	 * @return the number of primer pairs in this set
	 */
	public int getNumPrimerPairs(){
		return this.primerPairs.size();
	}
	
	/**
	 * Returns a specific primer pair.
	 * 
	 * @param index the index of the pair to return to returns
	 * @return the primer pair at position 'index' in this set
	 */
	public PrimerPair getPrimerPair(int index){
		if(index<0 || index >= this.primerPairs.size()) throw new IndexOutOfBoundsException();
		return (PrimerPair) this.primerPairs.getQuick(index);
	}
	
	/**
	 * Updates scores of this primer pair set.
	 */
	private void updateScores(){
		this.sumDistanceToOptimalPrimerPair = 0;
		for(int i=0; i<this.primerPairs.size(); i++){
			this.sumDistanceToOptimalPrimerPair += ((PrimerPair)this.primerPairs.getQuick(i)).getDistanceToOptimalPrimerPair();
		}
	}
	
	/**
	 * Returns the average distance to the virtual optimal primer pair.
	 * 
	 * @return the average distance to the virtual optimal primer pair
	 */
	public double getAvgDistOptPrimerPair(){
		return this.sumDistanceToOptimalPrimerPair/this.primerPairs.size();
	}
	
	/**
	 * Comoutes the homogenity of this primer pair set.
	 * 
	 * The homogenity is defined here to be the standard deviation of the scores of all primer pairs in this set.
	 *
	 */
	public void computeHomogenityScore(){
		this.homogenityScore = 0;
		this.updateScores();
		double temp;
		for(int i=0; i<this.primerPairs.size(); i++){
			temp = ((PrimerPair)this.primerPairs.getQuick(i)).getDistanceToOptimalPrimerPair() - this.getAvgDistOptPrimerPair();
			this.homogenityScore += (temp*temp);
		}
		this.homogenityScore = Math.sqrt(this.homogenityScore/ this.primerPairs.size());
		
		this.wasChanged = false;
	}
	
	/**
	 * Computes the maximum alignment score of all pairs in the current primer pair set.
	 * 
	 * The following combinations need to be checked, depending on the type of experiment:
	 * Experiment: Only forward1-forward2 and forward{1|2}-probe primer alignments are computed since reverse primers do not co-occur in the same reaction tube with a probe.
	 * Control: Additionally, within each pair, the following alignments need to be checked: forward-reverse, forward-probe, reverse-probe
	 *
	 */
	public void computeMaxAlignmentScore(){
		PrimerAlignmentScores scores;
		this.pa_max = -1;
		this.pea_max = -1;
		
		for(int i=0; i<this.primerPairs.size()-1; i++){
			scores = ((PrimerPair)this.primerPairs.getQuick(i)).getMaxPairElementsAlignmentScore();
			if(this.pa_max < scores.getPairScore()) this.pa_max = scores.getPairScore();
			if(this.pea_max < scores.getPairEndScore()) this.pea_max = scores.getPairEndScore();
			
			for(int j=i+1; j<this.primerPairs.size(); j++){
				scores = ((PrimerPair)this.primerPairs.getQuick(i)).getAlignmentValues((PrimerPair)this.primerPairs.getQuick(j));
				if(this.pa_max < scores.getPairScore()) this.pa_max = scores.getPairScore();
				if(this.pea_max < scores.getPairEndScore()) this.pea_max = scores.getPairEndScore();
			}
		}
		scores = ((PrimerPair)this.primerPairs.getQuick(this.primerPairs.size()-1)).getMaxPairElementsAlignmentScore();
		if(this.pa_max < scores.getPairScore()) this.pa_max = scores.getPairScore();
		if(this.pea_max < scores.getPairEndScore()) this.pea_max = scores.getPairEndScore();
	}
	
	public int getMaxPairAlignScore(){
		return this.pa_max;
	}
	
	public int getMaxPairAlignEndScore(){
		return this.pea_max;
	}
	
	public double getHomogenityScore(){
		if(this.wasChanged)	this.computeHomogenityScore();
		return this.homogenityScore;
	}
	
	public int size(){
		return this.primerPairs.size();
	}
	
	public double getAverageTMPrimers(){
		return this.averageTMPrimers / this.primerPairs.size();
	}
	
	public double getAverageTMProbes(){
		return this.averageTMProbes / this.primerPairs.size();
	}
	
	public String toString(){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<this.primerPairs.size(); i++){
			buffer.append(((PrimerPair)this.primerPairs.getQuick(i)).toString() + "\n\n");
		}
		buffer.append("avgdOptPrimerPAir: " + this.getAvgDistOptPrimerPair() + "\n");
		buffer.append("Homogenity score: " + this.getHomogenityScore() + "\n");
		buffer.append("worst primer-  primer PA/PEA score (incl probes iff present): " + this.getMaxPairAlignScore() + " " + this.getMaxPairAlignEndScore());
		return buffer.toString();
	}
	
	public String toFormattedString(){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<this.primerPairs.size(); i++){
			buffer.append(((PrimerPair)this.primerPairs.getQuick(i)).toFormattedString() + "\n\n");
		}
		buffer.append("avgdOptPrimerPAir: " + this.getAvgDistOptPrimerPair() + "\n");
		buffer.append("Homogenity score: " + this.getHomogenityScore() + "\n");
		buffer.append("worst primer-  primer PA/PEA score (incl probes iff present): " + this.getMaxPairAlignScore() + " " + this.getMaxPairAlignEndScore());
		return buffer.toString();
	}
	
	public boolean equals(PrimerPairSet other){
		if((this == null && other != null) || (this != null && other == null)) return false;
		ObjectArrayList a = new ObjectArrayList();
		a.addAllOfFromTo(this.primerPairs, 0, this.primerPairs.size() - 1);
		ObjectArrayList b = new ObjectArrayList();
		b.addAllOfFromTo(other.primerPairs, 0, other.primerPairs.size() -1);
		a.sort();
		b.sort();
		return a.equals(b);
	}
	
	/**
	 * Returns the sum of the hash codes of the sequences of all primers included in this primer pair set.
	 * 
	 * @return the sum of the hash codes of the sequences of all primers included in this primer pair set
	 */
	public int getHashCode(){
		int result = 0;
		
		for(int i=0; i<this.primerPairs.size(); i++){
			result += ((PrimerPair)this.primerPairs.getQuick(i)).getForwardPrimer().getSequence().hashCode();
			result += ((PrimerPair)this.primerPairs.getQuick(i)).getReversePrimer().getSequence().hashCode();
			if(this.searchParams.isPickTaqManProbe()) result += ((PrimerPair)this.primerPairs.getQuick(i)).getHybridizationProbe().getSequence().hashCode();
		}
		return result;
	}
}
