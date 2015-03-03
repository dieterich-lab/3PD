package primerDesign.dsc;

import java.text.NumberFormat;

/**
 * Wraps some primer pair picking statistics.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPairPickingStatistics {
		private long emptyBestPair = 0;
		private long noImprovementInWeightedScoreSum = 0;
		private long invalidPairSize = 0;
		private long noImprovementOfAlignmentScore = 0;
		private long acceptedPairSets = 0;
		private long pairSetsScannedInit = 0;
		private long pairSetsScannedCheck = 0;
		
		private long pairsScanned = 0;
		private long pairsAccepted = 0;
		private long pairsRejected = 0;
		private int intraPairAlignmentAccept = 0;
		private long intraPairAlignmentReject = 0;
		private long interPairAlignmentAccept = 0;
		private long interPairAlignmentReject = 0;
		
		private long primerPrimerTMaccept = 0;
		private long primerPrimerTMreject = 0;
		private long primerProbeTMaccept = 0;
		private long primerProbeTMreject = 0;
		
		private long pairPrimerPrimerTMaccept = 0;
		private long pairPrimerPrimerTMreject = 0;
		private long pairPrimerProbeTMaccept = 0;
		private long pairPrimerProbeTMreject = 0;
		
		private long PrimerPrimerAlignmentAccept = 0;
		private long primerPrimerAlignmentReject = 0;
		private long primerProbeAlignmentAccept = 0;
		private long primerProbeAlignmentReject = 0;
		
		private long validPrimerPairsEnum = 0;
		private long invalidPrimerPairsEnum = 0;
		
		/**
		 * Combines statistics objects.
		 * 
		 * @param other the other statistic to merge into this one
		 */
		public void combineStats(PrimerPairPickingStatistics other){
			this.emptyBestPair += other.emptyBestPair;
			this.noImprovementInWeightedScoreSum += other.noImprovementInWeightedScoreSum;
			this.invalidPairSize += other.invalidPairSize;
			this.noImprovementOfAlignmentScore += other.noImprovementOfAlignmentScore;
			this.acceptedPairSets += other.acceptedPairSets;
			this.pairSetsScannedInit += other.pairSetsScannedInit;
			this.pairSetsScannedCheck += other.pairSetsScannedCheck;
			
			this.pairsScanned += other.pairsScanned;
			this.pairsAccepted += other.pairsAccepted;
			this.pairsRejected += other.pairsRejected;
			this.intraPairAlignmentAccept += other.intraPairAlignmentAccept;
			this.intraPairAlignmentReject += other.intraPairAlignmentReject;
			this.interPairAlignmentAccept += other.interPairAlignmentAccept;
			this.interPairAlignmentReject += other.interPairAlignmentReject;
			
			this.primerPrimerTMaccept += other.primerPrimerTMaccept;
			this.primerPrimerTMreject += other.primerPrimerTMreject;
			this.primerProbeTMaccept += other.primerProbeTMaccept;
			this.primerProbeTMreject += other.primerProbeTMreject;
			
			this.pairPrimerPrimerTMaccept += other.pairPrimerPrimerTMaccept;
			this.pairPrimerPrimerTMreject += other.pairPrimerPrimerTMreject;
			this.pairPrimerProbeTMaccept += other.pairPrimerProbeTMaccept;
			this.pairPrimerProbeTMreject += other.pairPrimerProbeTMreject;
			
			this.validPrimerPairsEnum += other.validPrimerPairsEnum;
		}
		
		public long getEmptyBestPair() {
			return emptyBestPair;
		}

		public long getNoImprovementInWeightedScoreSum() {
			return noImprovementInWeightedScoreSum;
		}

		public long getInvalidPairSize() {
			return invalidPairSize;
		}

		public long getNoImprovementOfAlignmentScore() {
			return noImprovementOfAlignmentScore;
		}

		public long getAcceptedPairSets() {
			return acceptedPairSets;
		}

		public long getPairSetsScannedInit() {
			return pairSetsScannedInit;
		}
		
		public long getPairSetsScannedCheck(){
			return this.pairSetsScannedCheck;
		}

		public long getPairsScanned() {
			return pairsScanned;
		}

		public long getPairsAccepted() {
			return pairsAccepted;
		}

		public long getPairsRejected() {
			return pairsRejected;
		}

		public long getIntraPairAlignmentAccept() {
			return intraPairAlignmentAccept;
		}

		public long getIntraPairAlignmentReject() {
			return intraPairAlignmentReject;
		}

		public long getInterPairAlignmentAccept() {
			return interPairAlignmentAccept;
		}

		public long getIntPairAlignmentReject() {
			return interPairAlignmentReject;
		}

		public void incEmptyBestPair(){
			this.emptyBestPair++;
		}
		
		public void incNoImprovementInWeightedSum(){
			this.noImprovementInWeightedScoreSum++;
		}
		
		public void incInvalidPairSize(){
			this.invalidPairSize++;
		}
		
		public void incNoImprovementOfAlignmentScore(){
			this.noImprovementOfAlignmentScore++;
		}
		
		public void incAcceptedPairSets(){
			this.acceptedPairSets++;
		}
		
		public void incPairSetsScanned(){
			this.pairSetsScannedInit++;
		}
		
		public void incPairSetsScannedCheck(){
			this.pairSetsScannedCheck++;
		}
		
		public void incPairsScanned(){
			this.pairsScanned++;
		}
		
		public void incPairsAccepted(){
			this.pairsAccepted++;
		}
		
		public void incPairsRejected(){
			this.pairsRejected++;
		}
		
		public void incIntraPairAlignmentAccept(){
			this.intraPairAlignmentAccept++;
		}
		
		public void incIntraPairAlignmentReject(){
			this.intraPairAlignmentReject++;
		}
		
		public void incInterPairAlignmentAccept(){
			this.interPairAlignmentAccept++;
		}
		
		public void incInterPairAlignmentReject(){
			this.interPairAlignmentReject++;
		}
		
		public void incPrimerPrimerTMaccept(){
			this.primerPrimerTMaccept++;
		}
		
		public void incPrimerPrimerTMreject(){
			this.primerPrimerTMreject++;
		}
		
		public void incPrimerProbeTMaccept(){
			this.primerProbeTMaccept++;
		}
		
		public void incPrimerProbeTMreject(){
			this.primerProbeTMreject++;
		}
		
		public void incPrimerPrimerAlignmentAccept(){
			this.PrimerPrimerAlignmentAccept++;
		}
		
		public void incPrimerPrimerAlignmentReject(){
			this.primerPrimerAlignmentReject++;
		}
		
		public void incPrimerProbeAlignmentAccept(){
			this.primerProbeAlignmentAccept++;
		}
		
		public void incPrimerProbeAlignmentReject(){
			this.primerProbeAlignmentReject++;
		}
		
		public void incPairPrimerPrimerTMaccept(){
			this.pairPrimerPrimerTMaccept++;
		}
		
		public void incPairPrimerPrimerTMreject(){
			this.pairPrimerPrimerTMreject++;
		}
		
		public void incPairPrimerProbeTMaccept(){
			this.pairPrimerProbeTMaccept++;
		}
		
		public void incPairPrimerProbeTMreject(){
			this.pairPrimerProbeTMreject++;
		}
		
		public void incValidPrimerPairsEnum(){
			this.validPrimerPairsEnum++;
		}
		
		public void incInvalidPrimerPairsEnum(){
			this.invalidPrimerPairsEnum++;
		}
		
		public String printStat(){
			StringBuffer buffy = new StringBuffer();
			NumberFormat format = NumberFormat.getInstance();
			
			buffy.append("Primer pair picking statistics:\n");
			buffy.append("-------------------------------\n");
			buffy.append("Pair Primer-Primer max TM diff check: " + format.format(this.primerPrimerTMaccept + this.primerPrimerTMreject) + " (" + format.format(this.primerPrimerTMaccept) + " OK/ " + format.format(this.primerPrimerTMreject) + " NOT OK)\n");
			buffy.append("Pair Primer-Probe min TM diff check: " + format.format(this.primerProbeTMaccept + this.primerProbeTMreject) + " (" + format.format(this.primerProbeTMaccept) + " OK/ " + format.format(this.primerProbeTMreject) + " NOT OK)\n");
			buffy.append("Pair Primer-Primer alignment check: " + format.format(this.PrimerPrimerAlignmentAccept + this.primerPrimerAlignmentReject) + " (" + format.format(this.PrimerPrimerAlignmentAccept) + " OK/ " + format.format(this.primerPrimerAlignmentReject) + " NOT OK)\n");
			buffy.append("Pair Primer-Probe alignment check: " + format.format(this.primerProbeAlignmentAccept + this.primerProbeAlignmentReject) + " (" + format.format(this.primerProbeAlignmentAccept) + " OK/ " + format.format(this.primerProbeAlignmentReject) + " NOT OK)\n");
			buffy.append("Pair valid primer pairs enumerated: " + format.format(this.validPrimerPairsEnum) + "\n");
			buffy.append("Pair invalid primer pairs enumerated: " + format.format(this.invalidPrimerPairsEnum) + "\n");
			buffy.append("Empty best pair: " + format.format(this.emptyBestPair) + "\n");
			buffy.append("Invalid pair set size: " + format.format(this.invalidPairSize) + "\n");
			buffy.append("No improvement in weighted sum score: " + format.format(this.noImprovementInWeightedScoreSum) + "\n");
			buffy.append("No improvement in alignment scores: " + format.format(this.noImprovementOfAlignmentScore) + "\n");
			buffy.append("Accepted pair sets: " + format.format(this.acceptedPairSets) + "\n");
			buffy.append("Pair sets scanned (init): " + format.format(this.pairSetsScannedInit) + "\n");
			buffy.append("Pair sets scanned (check): " + format.format(this.pairSetsScannedCheck) + "\n");
			buffy.append("\tPairs scanned: " + format.format(this.pairsScanned) + "\n");
			buffy.append("\tPair-Pair Primers max TM diff check: " + format.format(this.pairPrimerPrimerTMaccept + this.pairPrimerPrimerTMreject) + " (" + format.format(this.pairPrimerPrimerTMaccept) + " OK/ " + format.format(this.pairPrimerPrimerTMreject) + " NOT OK)\n");
			buffy.append("\tPair-Pair Probes min TM diff check: " + format.format(this.pairPrimerProbeTMaccept + this.pairPrimerProbeTMreject) + " (" + format.format(this.pairPrimerProbeTMaccept) + " OK/ " + format.format(this.pairPrimerProbeTMreject) + " NOT OK)\n");
			buffy.append("\tintra pair alignment reject: " + format.format(this.intraPairAlignmentReject) + "\n");
			buffy.append("\tintra pair alignment accepted: " + format.format(this.intraPairAlignmentAccept) + "\n");
			buffy.append("\tinter pair alignment reject: " + format.format(this.interPairAlignmentReject) + "\n");
			buffy.append("\tinter pair alignment accept: " + format.format(this.interPairAlignmentAccept) + "\n");
			buffy.append("\tAcceptable pairs accepted (improvement): " + format.format(this.pairsAccepted) + "\n");
			buffy.append("\tAcceptable pairs rejected (no improvement): " + format.format(this.pairsRejected) + "\n");
			
			return buffy.toString();
		}
		
		/**
		 * Returns a statistics for use with the webserver to figure out which parameters to refine.
		 * 
		 * @return a string of a statistics for use with the webserver to figure out which parameters to refine
		 */
		public String printWebserviceStat(){
			StringBuffer buffy = new StringBuffer();
			NumberFormat format = NumberFormat.getInstance();
			
			buffy.append("Primer pair picking statistics:\n");
			buffy.append("-------------------------------\n");
			buffy.append("Pair Primer-Primer max TM difference check: " + format.format(this.primerPrimerTMaccept + this.primerPrimerTMreject) + " (" + format.format(this.primerPrimerTMaccept) + " OK/ " + format.format(this.primerPrimerTMreject) + " NOT OK)\n");
			buffy.append("Pair Primer-Probe min TM difference check: " + format.format(this.primerProbeTMaccept + this.primerProbeTMreject) + " (" + format.format(this.primerProbeTMaccept) + " OK/ " + format.format(this.primerProbeTMreject) + " NOT OK)\n");
			buffy.append("Pair Primer-Primer alignment check: " + format.format(this.PrimerPrimerAlignmentAccept + this.primerPrimerAlignmentReject) + " (" + format.format(this.PrimerPrimerAlignmentAccept) + " OK/ " + format.format(this.primerPrimerAlignmentReject) + " NOT OK)\n");
			buffy.append("Pair Primer-Probe alignment check: " + format.format(this.primerProbeAlignmentAccept + this.primerProbeAlignmentReject) + " (" + format.format(this.primerProbeAlignmentAccept) + " OK/ " + format.format(this.primerProbeAlignmentReject) + " NOT OK)\n");
			buffy.append("Pair valid primer pairs enumerated: " + format.format(this.validPrimerPairsEnum) + "\n");
			buffy.append("Pair invalid primer pairs enumerated: " + format.format(this.invalidPrimerPairsEnum) + "\n");
			buffy.append("Pairs scanned: " + format.format(this.pairsScanned) + "\n");
			buffy.append("\tPair-Pair Primers max TM diff check: " + format.format(this.pairPrimerPrimerTMaccept + this.pairPrimerPrimerTMreject) + " (" + format.format(this.pairPrimerPrimerTMaccept) + " OK/ " + format.format(this.pairPrimerPrimerTMreject) + " NOT OK)\n");
			buffy.append("\tPair-Pair Probes min TM diff check: " + format.format(this.pairPrimerProbeTMaccept + this.pairPrimerProbeTMreject) + " (" + format.format(this.pairPrimerProbeTMaccept) + " OK/ " + format.format(this.pairPrimerProbeTMreject) + " NOT OK)\n");
			
			return buffy.toString();
		}
		
		public String toString(){
			return this.printStat();
		}
}
