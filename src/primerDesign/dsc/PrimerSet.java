package primerDesign.dsc;

import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.util.PrimerSearchParameters;
import weka.core.FastVector;


/**
 * This class implements a set of PCR forwardPrimers.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerSet implements Comparable{
	private FastVector forwardPrimers;
	private double avgMeltingTemp;
	private double deltaMeltingTemp;
	private double avgGcContent;
	private double deltaGcContent;
	private int avgLength;
	private int deltaLength;
	private int avgDistanceToRSS;
	private int deltaDistanceToRSS;
	private double maxSelfAlignmentScore;
	private double maxSelfEndAlignmentScore;
	private int maxPairAlignmentScore;
	private int maxPairEndAlignmentScore;
	private int sumFalsePositives;
	private int lengthShortestFalsePositive;
	private double score;
	private PrimerSearchParameters searchParams;
	public static PrimerAlignmentCalculation alignment;
	
	/**
	 * Initializes a primer set.
	 * 
	 * @param searchParams the 3PD search parameters
	 */
	public PrimerSet(PrimerSearchParameters searchParams){
		this.forwardPrimers = new FastVector();
		this.avgMeltingTemp = 0;
		this.deltaMeltingTemp = 0;
		this.avgGcContent = 0;
		this.deltaGcContent = 0;
		this.avgLength = 0;
		this.deltaLength = 0;
		this.avgDistanceToRSS = 0;
		this.deltaDistanceToRSS = 0;
		this.maxSelfAlignmentScore = 0;
		this.maxSelfEndAlignmentScore = 0;
		this.maxPairAlignmentScore = 0;
		this.maxPairEndAlignmentScore = 0;
		this.sumFalsePositives = 0;
		this.lengthShortestFalsePositive = Integer.MAX_VALUE;
		this.score = 0;
		this.searchParams = searchParams;
		alignment = searchParams.getPRIMER_ALIGNMENT_METHOD();
	}
	
	/**
	 * Returns the difference in distance to the restriction site.
	 * 
	 * @return the difference in distance to the restriction site
	 */
	public int getDeltaDistanceToRSS() {
		return deltaDistanceToRSS;
	}
	
	/**
	 * Returns the difference in GC content.
	 * 
	 * @return the difference in GC content
	 */
	public double getDeltaGcContent() {
		return deltaGcContent;
	}
	
	/**
	 * @return the deltaLength
	 */
	public int getDeltaLength() {
		return deltaLength;
	}
	
	/**
	 * @return the deltaMeltingTemp
	 */
	public double getDeltaMeltingTemp() {
		return deltaMeltingTemp;
	}

	/**
	 * @return the maxSelfAlignmentScore
	 */
	public double getMaxSelfAlignmentScore() {
		return maxSelfAlignmentScore;
	}

	/**
	 * @return the maxSelfEndAlignmentScore
	 */
	public double getMaxSelfEndAlignmentScore() {
		return maxSelfEndAlignmentScore;
	}

	/**
	 * @return the lengthShortestFalsePositive
	 */
	public int getLengthShortestFalsePositive() {
		return lengthShortestFalsePositive;
	}

	/**
	 * @return the forwardPrimers
	 */
	public FastVector getForwardPrimers() {
		return forwardPrimers;
	}
	/**
	 * @param forwardPrimers the forwardPrimers to set
	 */

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @return the sumFalsePositives
	 */
	public int getSumFalsePositives() {
		return sumFalsePositives;
	}
	
	/**
	 * @return the maxPairAlignmentScore
	 */
	public double getDeltaPairAlignmentScore() {
		return maxPairAlignmentScore;
	}
	
	/**
	 * Computes the properties of a DNA primer set.
	 * 
	 * The properties of a DNA primer set determine the homogeneousness of this primerset.
	 *
	 */
	private void computeProperties(){
		
		computePairAlignments();
		
		this.avgMeltingTemp = 0;
		this.avgGcContent = 0;
		this.avgLength = 0;
		this.avgDistanceToRSS = 0;
		double falsePositivesMean = 0;
		
		this.deltaMeltingTemp = 0;
		this.deltaGcContent = 0;
		this.deltaLength = 0;
		this.deltaDistanceToRSS = 0;
		this.maxSelfAlignmentScore = 0;
		this.maxSelfEndAlignmentScore = 0;
		this.sumFalsePositives = 0;
		this.lengthShortestFalsePositive = Integer.MAX_VALUE;
		this.score = 0;
		
		Primer primer;
		int nb_elements = forwardPrimers.size();
		// compute means
		for(int i=0; i<nb_elements; i++){
			primer = (Primer) forwardPrimers.elementAt(i);
			this.avgMeltingTemp += primer.getMeltingTemp();
			this.avgGcContent += primer.getGcContent();
			this.avgLength+= primer.getLength();
			this.avgDistanceToRSS += primer.getDistanceToRSS();
			falsePositivesMean += primer.getFalsePositiveMatches();
			// lengthShortestFalsePositiveMean = 0;
			if(primer.getSelfAlignmentScore() > this.maxSelfAlignmentScore) this.maxSelfAlignmentScore = primer.getSelfAlignmentScore();
			if(primer.getSelfEndAlignmentScore() > this.maxSelfEndAlignmentScore) this.maxSelfEndAlignmentScore = primer.getSelfEndAlignmentScore();
		}
		
		// normalize means
		this.avgMeltingTemp /= nb_elements;
		this.avgGcContent /= nb_elements;
		this.avgLength /= nb_elements;
		this.avgDistanceToRSS /= nb_elements;
		falsePositivesMean /= nb_elements;
		// lengthShortestFalsePositiveMean /= nb_elements;
		
		// compute deltas
		for(int i=0; i<nb_elements; i++){
			primer = (Primer) forwardPrimers.elementAt(i);
			this.deltaMeltingTemp += Math.abs(primer.getMeltingTemp() - this.avgMeltingTemp);
			this.deltaGcContent += Math.abs(primer.getGcContent() - this.avgGcContent);
			this.deltaLength += Math.abs(primer.getLength() - this.avgLength);
			this.deltaDistanceToRSS += Math.abs(primer.getDistanceToRSS() - this.avgDistanceToRSS);
			this.sumFalsePositives += primer.getFalsePositiveMatches();
			// lengthShortestFalsePositive = 0;
		}
		
		// normalize deltas
		this.deltaMeltingTemp /= nb_elements;
		this.deltaGcContent /= nb_elements;
		this.deltaLength /= nb_elements;
		this.deltaDistanceToRSS /= nb_elements;
		// this.lengthShortestFalsePositive /= nb_elements;
		
		// compute pair alignments
		PrimerAlignmentScores scores = null;
		for(int i=0; i<forwardPrimers.size(); i++){
			for(int j=i+1; j<forwardPrimers.size(); j++){
				 scores = alignment.computePairAlignment(((Primer) forwardPrimers.elementAt(i)).getSequence(), ((Primer) forwardPrimers.elementAt(j)).getSequence());
				 if(scores.getPairScore() > this.maxPairAlignmentScore) this.maxPairAlignmentScore = scores.getPairScore();
				 if(scores.getPairEndScore() > this.maxPairEndAlignmentScore) this.maxPairEndAlignmentScore = scores.getPairEndScore();
				 
			}
		}
		
		// compute score of this primer set
		this.score += this.searchParams.getPRIMER_DELTA_TM_WEIGHT() * this.deltaMeltingTemp;
		this.score += this.searchParams.getPRIMER_DELTA_GC_WEIGHT() * this.deltaGcContent;
		this.score += this.searchParams.getPRIMER_DELTA_LENGTH_WEIGHT() * this.deltaLength;
		this.score += this.searchParams.getPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT() * this.deltaDistanceToRSS;
		this.score += this.searchParams.getSELF_ALIGNMENT_WEIGHT() * this.maxSelfAlignmentScore;
		this.score += this.searchParams.getSELF_END_ALIGNMENT_WEIGHT() * this.maxSelfEndAlignmentScore;
		this.score += this.searchParams.getPAIR_ALIGNMENT_WEIGHT() * this.maxPairAlignmentScore;
		this.score += this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT() * this.maxPairEndAlignmentScore;
		this.score += this.searchParams.getPRIMER_FALSE_POSITIVES_WEIGHT() * this.sumFalsePositives;
		
		// if any of the forwardPrimers in this primer set has a false positive match on the genome AND the amplicon resulting from this FP is < 'Constants.SAFE_FALSE_POSITIVE_AMPLICON_LENGTH' the score is 'Double.MAX_VALUE'
		// and this primer set is to be rejected!
		if(this.sumFalsePositives > 0 && this.lengthShortestFalsePositive < this.searchParams.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH()) this.score = Double.MAX_VALUE;
	}	
	
	/**
	 * Adds the primer that is closest to this primer set to this primer set.
	 * 
	 * @param primers a list of primers to pick the closest primer from
	 * @param primerType the type of the primer to be added
	 */
	public void addBestPrimer(Primer[] primers, Enum primerType){
		double minScore = Double.MAX_VALUE;
		Primer bestPrimer = null;
		for(int i=0; i<primers.length; i++){
			Primer currentPrimer = (Primer) primers[i];
			double currentScore = this.scorePrimerSetTo(currentPrimer);
			if(currentScore < minScore){
				bestPrimer = currentPrimer;
				minScore = currentScore;
			}
		}
		if(primerType.equals(PrimerTypes.forwardPrimer)) this.forwardPrimers.addElement(bestPrimer);
		//else if(primerType.equals(PrimerTypes.reversePrimer))
		else throw new IllegalArgumentException("Unsupported primer type: " + primerType);
		this.computeProperties();
	}
	
	/**
	 * Adds a forward primer to this primer set.
	 * 
	 * @param primer a forward primer to be added to this primer set
	 */
	public void addForwardPrimer(Primer primer){
		this.forwardPrimers.addElement(primer);
		this.computeProperties();
	}
	
	/**
	 * Returns the valid forward primer at position 'index'.
	 * 
	 * @param index the index of the forward primer to be returned
	 * 
	 * @return the valid primer at position 'index'
	 */
	public Primer getForwardPrimer(int index){
		if(index < this.forwardPrimers.size())	return (Primer) this.forwardPrimers.elementAt(index);
		else throw new IllegalArgumentException("Invalid primer index!");
	}
	
	/**
	 * Returns the number of forward primer in this primer set.
	 * 
	 * @return the number of forward primer in this primer set
	 */
	public int getNumberOfForwardPrimers(){
		return this.forwardPrimers.size();
	}
	
	/**
	 * Returns the score of a primer set compared to another primer.
	 * 
	 * @param otherPrimer the other primer
	 * 
	 * @return the score of a primer set compared to another primer
	 */
	public double scorePrimerSetTo(Primer otherPrimer){
		double score = 0;
		
		int max_pa = this.maxPairAlignmentScore;
		int max_pea = this.maxPairEndAlignmentScore;
		for (int i = 0; i < forwardPrimers.size(); i++) {
			PrimerAlignmentScores scores = PrimerSet.alignment.computePairAlignment(((Primer)forwardPrimers.elementAt(i)).getSequence(), otherPrimer.getSequence());
			if(scores.getPairScore() > max_pa) max_pa = scores.getPairScore();
			if(scores.getPairEndScore() > max_pea) max_pea = scores.getPairEndScore();
		}
		
		score += this.searchParams.getPRIMER_DELTA_TM_WEIGHT() * Math.abs(this.avgMeltingTemp - otherPrimer.getMeltingTemp());
		score += this.searchParams.getPRIMER_DELTA_GC_WEIGHT() * Math.abs(this.avgGcContent - otherPrimer.getGcContent());
		score += this.searchParams.getPRIMER_DELTA_LENGTH_WEIGHT() * Math.abs(this.avgLength - otherPrimer.getLength());
		score += this.searchParams.getPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT() * Math.abs(this.avgDistanceToRSS - otherPrimer.getDistanceToRSS());
		score += this.searchParams.getSELF_ALIGNMENT_WEIGHT() * Math.max(this.maxSelfAlignmentScore, otherPrimer.getSelfAlignmentScore());
		score += this.searchParams.getSELF_END_ALIGNMENT_WEIGHT() * Math.max(this.maxSelfEndAlignmentScore, otherPrimer.getSelfEndAlignmentScore());
		score += this.searchParams.getPAIR_ALIGNMENT_WEIGHT() * max_pa;
		score += this.searchParams.getPAIR_END_ALIGNMENT_WEIGHT() * max_pea;
		score += this.searchParams.getPRIMER_FALSE_POSITIVES_WEIGHT() * (this.sumFalsePositives + otherPrimer.getFalsePositiveMatches());
		
		return score;
	}
	
	/**
	 * Computes primer pair alignments and primer pair end alignments.
	 *
	 */
	private void computePairAlignments(){
		if(forwardPrimers.size() > 1){
			for(int i=0; i < forwardPrimers.size()-1; i++){
				String primer1 = ((Primer) forwardPrimers.elementAt(i)).getSequence();
				for (int j = i+1; j < forwardPrimers.size(); j++) {
					PrimerSet.alignment.computePairAlignment(primer1, ((Primer) forwardPrimers.elementAt(j)).getSequence());
					int pa = PrimerSet.alignment.getPAScore();
					int pea = PrimerSet.alignment.getPEAScore();
					if(pa > this.maxPairAlignmentScore) this.maxPairAlignmentScore = pa;
					if(pea > this.maxPairEndAlignmentScore) this.maxPairEndAlignmentScore = pea;
				}
			}
		}
	}

	public int compareTo(Object o) {
		PrimerSet other = (PrimerSet) o;
		if(this.score < other.score) return -1;
		else if(this.score > other.score) return 1;
		else return 0;
	}

	/**
	 * @return the maxPairAlignmentScore
	 */
	public int getMaxPairAlignmentScore() {
		return maxPairAlignmentScore;
	}

	/**
	 * @return the maxPairEndAlignmentScore
	 */
	public int getMaxPairEndAlignmentScore() {
		return maxPairEndAlignmentScore;
	}
}
