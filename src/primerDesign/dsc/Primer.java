package primerDesign.dsc;

import java.util.HashMap;

import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.algo.PrimerMeltingTempCalculation;
import primerDesign.algo.RSSDPAligner;
import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;



/**
 * This class implements a PCR Primer.
 * 
 * @author Sebastian Fršhler
 *
 */
public class Primer implements Comparable{

	private String sequence; // the sequence must only be empty if this primer represents the virtual optimal primer!
	private double meltingTemp;
	private double gcContent;
	private int length;
	private int distanceToRSS;
	private int relativePosition;
	private double selfAlignmentScore;
	private double selfEndAlignmentScore;
	private int falsePositiveMatches;
	private double distanceToOptimalPrimer;
	private boolean isVirtualOptimalPrimer;
	private Enum<PrimerTypes> primerType;
//	private RSSDPAligner scanRegionAligner;
	private RestrictionSite restrictionSite;
	private int positionInScanSequence; // the relative! position within the sequence this primer was found in
	private PrimerAcceptanceLevel acceptanceLevel;
	
	private static PrimerAlignmentCalculation ALIGNMENT;
	private static PrimerMeltingTempCalculation TMCALCMETHOD = Constants.PRIMER_TM_CALC_METHOD;
	private static HashMap<String, Double> gc_content_hash = new HashMap<String, Double>();
	private static HashMap<String, Double> melting_temp_hash = new HashMap<String, Double>();
	
	// assign maximum values for score normalization such that each component of the score
	// scores equally when equal weights are assigned
	private static double DELTA_TM_MAX;
	private static double DELTA_TAQMAN_TM_MAX;
	private static double DELTA_GC_MAX;
	private static double DELTA_TAQMAN_GC_MAX;
	private static int DELTA_LENGTH_MAX;
	private static int DELTA_DIST_RSS_MAX;
	private static int DELTA_PRIMER_SA_MAX;
	private static int DELTA_PRIMER_SEA_MAX;
	private static int DELTA_PRIMER_PA_MAX;
	private static int DELTA_PRIMER_PEA_MAX;
	private static int DELTA_TAQMAN_SA_MAX;
	private static int DELTA_TAQMAN_SEA_MAX;
	private static int DELTA_TAQMAN_PA_MAX;
	private static int DELTA_TAQMAN_PEA_MAX;
	private static int MAX_FP;
	private static boolean isInitialized = false;
	
	/**
	 * Creates a new primer.
	 * 
	 * Primer sequence is assumed to be in 5'->3' direction!
	 * Primer alignments are computed as specified in 'Constants'
	 * 
	 * @param sequence
	 */
	public Primer(String sequence, Enum<PrimerTypes> primerType, PrimerSearchParameters searchParams){
		if(!isInitialized){
			initDeltas(searchParams);
			initFunctions(searchParams);
		}
		this.primerType = primerType;
		this.acceptanceLevel = PrimerAcceptanceLevel.NOT_TESTED;
		this.setSequence(sequence, true, searchParams);
		// computeProperties called by setSequence
		this.isVirtualOptimalPrimer = false;
	}
	
	/**
	 * Initializes a primer.
	 * 
	 * @param sequence the primer sequence
	 * @param primerType the primer type
	 * @param relPos the relative position of the primer within the scan region
	 * @param searchParams the 3PD primer search parameters
	 */
	public Primer(String sequence, Enum<PrimerTypes> primerType, int relPos, PrimerSearchParameters searchParams){
		if(!isInitialized){
			initDeltas(searchParams);
			initFunctions(searchParams);
		}
		this.primerType = primerType;
		this.acceptanceLevel = PrimerAcceptanceLevel.NOT_TESTED;
		this.setSequence(sequence, -1, relPos, 0, true, searchParams);
		// computeProperties called by setSequence
		this.isVirtualOptimalPrimer = false;
	}
	
	/**
	 * Initializes a primer.
	 * 
	 * @param sequence the primer sequence
	 * @param primerType the primer type
	 * @param relPos the relative position of the primer within the scan region
	 * @param sa the self alignment value of this primer
	 * @param sea the self-end alignment value of this primer
	 * @param searchParams the 3PD primer search parameters
	 */
	public Primer(String sequence, Enum<PrimerTypes> primerType, int relPos, int sa, int sea, PrimerSearchParameters searchParams){
		if(!isInitialized){
			initDeltas(searchParams);
			//initFunctions(searchParams);
		}
		
		this.primerType = primerType;
		this.acceptanceLevel = PrimerAcceptanceLevel.NOT_TESTED;
		this.selfAlignmentScore = sa;
		this.selfEndAlignmentScore = sea;
		this.setSequence(sequence, -1, relPos, 0, false, searchParams);
		this.isVirtualOptimalPrimer = false;
	}
	
	/**
	 * Initializes a primer.
	 * 
	 * @param sequence the primer sequence
	 * @param primerType the primer type
	 * @param scanRegionAligner the alignment method to be used for scan region alignments
	 * @param relPos the relative position of the primer within the scan region
	 * @param searchParams the 3PD primer search parameters
	 */
	public Primer(String sequence, Enum<PrimerTypes> primerType, RSSDPAligner scanRegionAligner, int relPos, PrimerSearchParameters searchParams){
		if(!isInitialized){
			initDeltas(searchParams);
			initFunctions(searchParams);
		}
		this.primerType = primerType;
		this.acceptanceLevel = PrimerAcceptanceLevel.NOT_TESTED;
//		this.scanRegionAligner = scanRegionAligner;
		this.setSequence(sequence, -1, relPos, 0, true, searchParams);
		// computeProperties called by setSequence
		this.isVirtualOptimalPrimer = false;
	}
	
	private Primer(PrimerSearchParameters searchParams){
		if(!isInitialized){
			initDeltas(searchParams);
			initFunctions(searchParams);
		}
		this.acceptanceLevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	private void initFunctions(PrimerSearchParameters searchParams){
		ALIGNMENT = searchParams.getPRIMER_ALIGNMENT_METHOD();
	}
	
	private void initDeltas(PrimerSearchParameters searchParams){
		// assign maximum values for score normalization such that each component of the score
		// scores equally when equal weights are assigned
		// except TM and GC values, all other values are integers, therefore division is bottomed at 1!
		// for TM and GC, division is bottomed by a virtualMinValue
		double virtualMinValue = 1E-09;
		DELTA_TM_MAX = Math.max(searchParams.getMAX_TM() - searchParams.getMIN_TM(), virtualMinValue);
		DELTA_TAQMAN_TM_MAX = Math.max(searchParams.getTAQMAN_MAX_TM() - searchParams.getTAQMAN_MIN_TM(), virtualMinValue);
		DELTA_GC_MAX = Math.max(searchParams.getMAX_GC() - searchParams.getMIN_GC(), virtualMinValue);
		DELTA_TAQMAN_GC_MAX = Math.max(searchParams.getTAQMAN_MAX_GC() - searchParams.getTAQMAN_MIN_GC(), virtualMinValue);
		DELTA_LENGTH_MAX = Math.max(searchParams.getMAX_PRIMER_LENGTH() - searchParams.getMIN_PRIMER_LENGTH(),1);
		DELTA_DIST_RSS_MAX = Math.max(searchParams.getMAX_AMPLICON_LENGTH()/2 - searchParams.getMIN_AMPLICON_LENGTH()/2,1);
		DELTA_PRIMER_SA_MAX = Math.max(searchParams.getMAX_PRIMER_SELF_ALIGNMENT_SCORE(),1);
		DELTA_PRIMER_SEA_MAX = Math.max(searchParams.getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE(),1);
		DELTA_PRIMER_PA_MAX = Math.max(searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE(),1);
		DELTA_PRIMER_PEA_MAX = Math.max(searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE(),1);
		DELTA_TAQMAN_SA_MAX = Math.max(searchParams.getMAX_TAQMAN_SELF_ALIGNMENT_SCORE(),1);
		DELTA_TAQMAN_SEA_MAX = Math.max(searchParams.getMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE(),1);
		DELTA_TAQMAN_PA_MAX = Math.max(searchParams.getMAX_TAQMAN_PAIR_ALIGNMENT_SCORE(),1);
		DELTA_TAQMAN_PEA_MAX = Math.max(searchParams.getMAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE(),1);
		MAX_FP = Math.max(searchParams.getMAX_PRIMER_MISPRIMING_CUTOFF(),1);
	}
	
	/**
	 * Computes/ sets all properties of a (new/modified) primer.
	 *
	 */
	private void computeProperties(boolean computeAlignments, PrimerSearchParameters searchParams){
		this.computeProperties(-1, -1, 0, computeAlignments, searchParams);
	}
	
	private void computeProperties(int dRSS, int relPos, int fpMatches, boolean computeAlignments, PrimerSearchParameters searchParams){
		this.length = this.sequence.length();
		this.computeMeltingTemp(searchParams);
		this.computeGCContent(searchParams);
		this.distanceToRSS = dRSS;
		this.relativePosition = relPos;
		if(computeAlignments) this.computeAlignments();
		this.falsePositiveMatches = fpMatches;
		this.computeDistanceToOptimalPrimer(searchParams);
	}
	
	/**
	 * Computes the GC content of the primer.
	 *
	 */
	private void computeGCContent(PrimerSearchParameters searchParams){
		if(searchParams.isUseIndices() && gc_content_hash.containsKey(this.sequence)){
			this.gcContent =  gc_content_hash.get(this.sequence);
		}
		else{
			double gcContent = 0;
			int maskedCharacters = 0; 
			char maskingChar = Constants.REPETITIVE_ELEMENT_CHARACTER.toCharArray()[0];
			for(int i=0; i < this.getLength(); i++){
				char base = sequence.charAt(i);
				if (base == 'G' || base == 'C'){
					gcContent++;
				}
				else if(base == maskingChar){
					maskedCharacters++;
				}
			}
			this.gcContent = gcContent/ (this.getLength() - maskedCharacters);
			if(searchParams.isUseIndices()) gc_content_hash.put(this.sequence, this.gcContent);
		}
	}
	
	/**
	 * Computes the melting temperature of the primer.
	 *
	 */
	private void computeMeltingTemp(PrimerSearchParameters searchParams){
		if(searchParams.isUseIndices() && Primer.melting_temp_hash.containsKey(this.sequence)){
			this.meltingTemp = Primer.melting_temp_hash.get(this.sequence);
		}else{
			// if short oligo: use 'cheap' Wallace method
			// Reference: Suggs et.al. 1981
			// (2¡C per {A,T}, 4¡C per {G,C})
			// http://en.wikipedia.org/wiki/DNA_melting
			// This formula is valid for oligos <14 bases and assumes that the reaction is carried out in the presence of 50mM monovalent cations. 
			// Source: http://www.promega.com/biomath/calc11.htm#disc
			if(this.getLength() <= 16){
				int a=0,t=0,g=0,c=0;
				for(int i=0; i< this.getLength(); i++){
					char base = this.getSequence().charAt(i);
					switch(base){
						case 'A': a++; break;
						case 'T': t++; break;
						case 'G': g++; break;
						case 'C': c++; break;					
					}
				}
				this.meltingTemp = 2*(a+t) + 4*(g+c);
			}
			// else use more accurate nearest-neighbor method, e.g.
			// SantaLucia JR. (1998). A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci., 95, 1460-65."
			else{
				//this.meltingTemp = Primer.lucia.computeTM(this.getSequence());
				if(this.primerType.equals(PrimerTypes.hybridizationProbe)) this.meltingTemp = Primer.TMCALCMETHOD.computeMonoAndDivalentCationCorrectedTM(this.getSequence(), searchParams.getTAQMAN_PROBE_CONCENTRATION(), searchParams.getMONOVALENT_CATION_CONCENTRATION(), searchParams.getDIVALENT_CATION_CONCENTRATION(), searchParams.getDNTP_CONCENTRATION());
				else this.meltingTemp = Primer.TMCALCMETHOD.computeMonoAndDivalentCationCorrectedTM(this.getSequence(), searchParams.getPRIMER_CONCENTRATION(), searchParams.getMONOVALENT_CATION_CONCENTRATION(), searchParams.getDIVALENT_CATION_CONCENTRATION(), searchParams.getDNTP_CONCENTRATION());
			}
			if(searchParams.isUseIndices()) Primer.melting_temp_hash.put(this.sequence, this.meltingTemp);
		}
	}
	
	/**
	 * Self alignment and self end alignment of a primer are computed.
	 *
	 */
	private void computeAlignments(){
		PrimerAlignmentScores scores = Primer.ALIGNMENT.computeSelfAlignment(this.getSequence()); // this.scanRegionAligner.computeSelfAlignment(this);
		this.selfAlignmentScore = scores.getPairScore();
		this.selfEndAlignmentScore = scores.getPairEndScore();
	}
	
	/**
	 * Computes the distance to a virtual optimal primer with optimal values
	 * as specified in Constants.
	 *
	 */
	private void computeDistanceToOptimalPrimer(PrimerSearchParameters searchParams){
		if(this.primerType.equals(PrimerTypes.forwardPrimer) || this.primerType.equals(PrimerTypes.reversePrimer)){
			this.distanceToOptimalPrimer = this.scoreTo(Primer.getVirtualOptimalPrimer(searchParams), searchParams);
		}
		else if(this.primerType.equals(PrimerTypes.hybridizationProbe)){
			this.distanceToOptimalPrimer = this.scoreTo(Primer.getVirtualOptimalProbe(searchParams), searchParams);
		}
		else throw new IllegalStateException("Unhandled case!");
	}
	
	/**
	 * Returns the distance of the Primer to the restriction site.
	 * 
	 * @return the distance of the Primer to the restriction site
	 */
	public int getDistanceToRSS() {
		return distanceToRSS;
	}
	/**
	 * Sets the distance of the 5' end of the primer to the restriction site.
	 * 
	 * @param distanceToRSS the distance of the 5' end of the Primer to the restriction site
	 */
	public void setDistanceToRSS(int distanceToRSS) {
		if(distanceToRSS >= 0){
			this.distanceToRSS = distanceToRSS;
		}
		else throw new IllegalArgumentException("Distance to restriction site must be >= 0!");
	}
	/**
	 * Gets the relative position of the primer in the target genomic region.
	 * 
	 * @return the relative position of the primer in the target genomic region
	 */
	public int getRelativePosition() {
		return this.relativePosition;
	}

	/**
	 * Sets the relative position of the primer in the target genomic region.
	 * 
	 * @param relativePosition the relative position of the primer in the target genomic region
	 */
	public void setRelativePosition(int relativePosition) {
		if(relativePosition >= 0){
			this.relativePosition = relativePosition;
		}
		else throw new IllegalArgumentException("The relative position of the primer must be >=0");
	}

	/**
	 * Returns the number of false-positive matches of the Primer.
	 * 
	 * @return the number of false-positive matches of the Primer
	 */
	public int getFalsePositiveMatches() {
		return falsePositiveMatches;
	}
	/**
	 * Sets the number of false-positive matches of the Primer.
	 * 
	 * @param falsePositiveMatches the number of false-positive matches of the Primer
	 */
	public void setFalsePositiveMatches(int falsePositiveMatches) {
		if(falsePositiveMatches >= 0){
			this.falsePositiveMatches = falsePositiveMatches;
		}
		else throw new IllegalArgumentException("falsePositiveMatches must be >= 0!");
	}
	/**
	 * Returns the GC content of the Primer.
	 * 
	 * @return the GC content of the Primer
	 */
	public double getGcContent() {
		return gcContent;
	}
	/**
	 * Returns the length of the Primer.
	 * 
	 * @return the length of the Primer
	 */
	public int getLength() {
		return length;
	}
	/**
	 * Returns the melting temperature of the Primer.
	 * 
	 * @return the melting temperature of the Primer
	 */
	public double getMeltingTemp() {
		return meltingTemp;
	}
	/**
	 * Returns the self-alignment score of the Primer.
	 * 
	 * @return the self-alignment score of the Primer
	 */
	public double getSelfAlignmentScore() {
		return selfAlignmentScore;
	}
	/**
	 * Returns the self-end-alignment score of the Primer.
	 * 
	 * @return the self-end-alignment score of the Primer
	 */
	public double getSelfEndAlignmentScore() {
		return selfEndAlignmentScore;
	}
	/**
	 * Returns the sequence of the Primer.
	 * 
	 * @return the sequence of the Primer
	 */
	public String getSequence() {
		if(sequence.equals("")) throw new IllegalStateException("No sequence can be queried for the virtual optimal primer!");
		return sequence;
	}
	/**
	 * Sets the sequence the sequence of the Primer.
	 * 
	 * @param sequence the sequence of the Primer
	 * 
	 * @throws IllegalArgumentException if sequence length is < 2
	 */
	public void setSequence(String sequence, boolean computeAlignments, PrimerSearchParameters searchParams) {
		this.setSequence(sequence, -1, -1, 0, computeAlignments, searchParams);
	}
	
	/**
	 * Returns a subsequence of the primer (suffix) starting at position 'i'.
	 * 
	 * @param start the position
	 * 
	 * @return a subsequence of the primer (suffix) starting at position 'i'
	 */
	public String getSubsequence(int start){
		return this.sequence.substring(start);
	}
	
	/**
	 * Returns the primer's sequence length.
	 * 
	 * @return the primer's sequence length
	 */
	public int getSequenceLength(){
		return this.sequence.length();
	}
	
	/**
	 * Sets the primer sequence.
	 * 
	 * @param sequence the primer sequence
	 * @param dRSS the distance to the restriction site
	 * @param relPos the relative position of the primer within the scan region
	 * @param fpMatches the number of false positive matches of this primer
	 * @param computeAlignments whether to compute alignments for this primer
	 * @param searchParams the 3PD primer search params
	 */
	private void setSequence(String sequence, int dRSS, int relPos, int fpMatches, boolean computeAlignments, PrimerSearchParameters searchParams){
		if(sequence.length() >= 2){
			this.sequence = sequence.toUpperCase();
			this.computeProperties(dRSS, relPos, fpMatches, computeAlignments, searchParams);
		}else throw new IllegalArgumentException("Sequence length must be >= 2!");
	}

	/**
	 * Returns the distance to a (virtual) optimal primer as specified in 'Constants'.
	 * 
	 * @return the distance to a (virtual) optimal primer as specified in 'Constants'
	 */
	public double getDistanceToOptimalPrimer() {
		return distanceToOptimalPrimer;
	}
	
	/**
	 * Returns the type of this primer.
	 * 
	 * @return the type of this primer
	 */
	public Enum<PrimerTypes> getPrimerType(){
		return this.primerType;
	}

	public int compareTo(Object o) {
		if(this.distanceToOptimalPrimer < ((Primer) o).distanceToOptimalPrimer) return -1;
		else if(this.distanceToOptimalPrimer > ((Primer) o).distanceToOptimalPrimer) return 1;
		else return 0;
	}
	
	public boolean equals(Primer otherPrimer){
		if(this.sequence.equals(otherPrimer.sequence) && this.meltingTemp == otherPrimer.meltingTemp
				&& this.gcContent == otherPrimer.gcContent && this.length == otherPrimer.length
				&& this.distanceToRSS == otherPrimer.distanceToRSS && this.distanceToOptimalPrimer == otherPrimer.distanceToOptimalPrimer
				&& this.selfAlignmentScore == otherPrimer.selfAlignmentScore && this.selfEndAlignmentScore == otherPrimer.selfEndAlignmentScore
				&& this.falsePositiveMatches == otherPrimer.falsePositiveMatches) return true;
		else return false;
	}
	
	/**
	 * Scores one primer to another of equal primer type!
	 * 
	 * @param other the other primer of equal primer type to score this primer to
	 * 
	 * @return the score of this primer to the other primer of equal type
	 */
	public double scoreTo(Primer other, PrimerSearchParameters searchParams){
		PrimerAlignmentScores scores;
		int items = 0;
		
		if(other.isVirtualOptimalPrimer || this.isVirtualOptimalPrimer) scores = PrimerAlignmentScores.getOptimalAlignmentScores();
		else scores = ALIGNMENT.computePairAlignment(this.sequence, other.getSequence());
		
		double score = 0;
		
		if(PrimerTypes.isCompatible(this.primerType, other.getPrimerType())){
			if(this.primerType.equals(PrimerTypes.hybridizationProbe) && other.getPrimerType().equals(PrimerTypes.hybridizationProbe)){
				score += searchParams.getPRIMER_DELTA_TM_WEIGHT() * Math.abs(this.meltingTemp - other.getMeltingTemp()) / Primer.DELTA_TAQMAN_TM_MAX
					+ searchParams.getPRIMER_DELTA_GC_WEIGHT() * Math.abs(this.gcContent - other.getGcContent()) / Primer.DELTA_TAQMAN_GC_MAX
					+ searchParams.getSELF_ALIGNMENT_WEIGHT() * Math.max(this.selfAlignmentScore, other.getSelfAlignmentScore()) / Primer.DELTA_TAQMAN_SA_MAX
					+ searchParams.getSELF_END_ALIGNMENT_WEIGHT() * Math.max(this.selfEndAlignmentScore, other.getSelfEndAlignmentScore()) / Primer.DELTA_TAQMAN_SEA_MAX
					+ searchParams.getPAIR_ALIGNMENT_WEIGHT() * scores.getPairScore() / Primer.DELTA_TAQMAN_PA_MAX
					+ searchParams.getPAIR_END_ALIGNMENT_WEIGHT() * scores.getPairEndScore() / Primer.DELTA_TAQMAN_PEA_MAX;
				items += 6;
			}else{
				score += searchParams.getPRIMER_DELTA_TM_WEIGHT() * Math.abs(this.meltingTemp - other.getMeltingTemp()) / Primer.DELTA_TM_MAX
					+ searchParams.getPRIMER_DELTA_GC_WEIGHT() * Math.abs(this.gcContent - other.getGcContent()) / Primer.DELTA_GC_MAX
					+ searchParams.getPRIMER_FALSE_POSITIVES_WEIGHT() * (this.falsePositiveMatches + other.getFalsePositiveMatches()) / Primer.MAX_FP
					+ searchParams.getPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT() * Math.abs(this.distanceToRSS - other.getDistanceToRSS()) / Primer.DELTA_DIST_RSS_MAX
					+ searchParams.getPRIMER_DELTA_LENGTH_WEIGHT() * Math.abs(this.length - other.getLength()) / Primer.DELTA_LENGTH_MAX
					+ searchParams.getSELF_ALIGNMENT_WEIGHT() * Math.max(this.selfAlignmentScore, other.getSelfAlignmentScore()) / Primer.DELTA_PRIMER_SA_MAX
					+ searchParams.getSELF_END_ALIGNMENT_WEIGHT() * Math.max(this.selfEndAlignmentScore, other.getSelfEndAlignmentScore()) / Primer.DELTA_PRIMER_SEA_MAX
					+ searchParams.getPAIR_ALIGNMENT_WEIGHT() * scores.getPairScore() / Primer.DELTA_PRIMER_PA_MAX
					+ searchParams.getPAIR_END_ALIGNMENT_WEIGHT() * scores.getPairEndScore() / Primer.DELTA_PRIMER_PEA_MAX;
				items += 9;
			}
		}
		else if(other.primerType.equals(PrimerTypes.hybridizationProbe)){
			score += Primer.scoreToHybProbe(this, other, scores, searchParams);
			items += 6;
		}
		else if(this.primerType.equals(PrimerTypes.hybridizationProbe)){
			score += Primer.scoreToHybProbe(other, this, scores, searchParams);
			items += 6;
		}
		else throw new IllegalArgumentException("Primer types " + this.primerType + " - " + other.primerType + " cannot be scored - not comparable!");
		
//		double score = Constants.PRIMER_DELTA_TM_WEIGHT * Math.abs(this.meltingTemp - other.getMeltingTemp()) / Primer.DELTA_TM_MAX
//						+ Constants.PRIMER_DELTA_GC_WEIGHT * Math.abs(this.gcContent - other.getGcContent()) / Primer.DELTA_GC_MAX
//						+ Constants.PRIMER_DELTA_LENGTH_WEIGHT * Math.abs(this.length - other.getLength()) / Primer.DELTA_LENGTH_MAX
//						+ Constants.PRIMER_FALSE_POSITIVES_WEIGHT * (this.falsePositiveMatches + other.getFalsePositiveMatches()) / Primer.MAX_FP;;
//		
//		if(this.primerType.equals(PrimerTypes.hybridizationProbe) || other.getPrimerType().equals(PrimerTypes.hybridizationProbe)){
//			score += Constants.PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT * Math.abs(this.distanceToRSS - other.getDistanceToRSS()) / Primer.DELTA_DIST_RSS_MAX;
//		}
//						
//		if(this.primerType.equals(PrimerTypes.hybridizationProbe) && other.getPrimerType().equals(PrimerTypes.hybridizationProbe)){
//			score += Constants.SELF_ALIGNMENT_WEIGHT * Math.max(this.selfAlignmentScore, other.getSelfAlignmentScore()) / Primer.DELTA_TAQMAN_SA_MAX
//			+ Constants.SELF_END_ALIGNMENT_WEIGHT * Math.max(this.selfEndAlignmentScore, other.getSelfEndAlignmentScore()) / Primer.DELTA_TAQMAN_SEA_MAX
//			+ Constants.PAIR_ALIGNMENT_WEIGHT * scores.getPairScore() / Primer.DELTA_TAQMAN_PA_MAX
//			+ Constants.PAIR_END_ALIGNMENT_WEIGHT * scores.getPairEndScore() / Primer.DELTA_TAQMAN_PEA_MAX;
//		}else{
//			score += Constants.SELF_ALIGNMENT_WEIGHT * Math.max(this.selfAlignmentScore, other.getSelfAlignmentScore()) / Primer.DELTA_PRIMER_SA_MAX
//				+ Constants.SELF_END_ALIGNMENT_WEIGHT * Math.max(this.selfEndAlignmentScore, other.getSelfEndAlignmentScore()) / Primer.DELTA_PRIMER_SEA_MAX
//				+ Constants.PAIR_ALIGNMENT_WEIGHT * scores.getPairScore() / Primer.DELTA_PRIMER_PA_MAX
//				+ Constants.PAIR_END_ALIGNMENT_WEIGHT * scores.getPairEndScore() / Primer.DELTA_PRIMER_PEA_MAX;
//		}
		
		return score/ items;
	}
	
	/**
	 * Scores a forward or reverse primer 'first' to a hybridization probe 'second'.
	 * 
	 * W.r.t alignments, primer-probe pairs are scored by the maximum primer-primer score thresholds!
	 * 
	 * @param primer the forward or reverse primer
	 * @param probe the hybridization probe
	 * @param scores the primer alignment scores between 'primer' and 'probe'
	 * 
	 * @return the score between the forward|reverse primer and the hybridization probe
	 */
	private static double scoreToHybProbe(Primer primer, Primer probe, PrimerAlignmentScores scores, PrimerSearchParameters searchParams){
		assert(probe.getPrimerType().equals(PrimerTypes.hybridizationProbe));
		double score = 0;
		if(probe.getMeltingTemp() < primer.getMeltingTemp() + searchParams.getMIN_TAQMAN_TM_DIFFERENCE()){
			score += searchParams.getPRIMER_DELTA_TM_WEIGHT() * 1;  // TM differences below MIN_TAQMAN_TM_DIFFERENCE always score maximum!
		}
		score += searchParams.getPRIMER_DELTA_GC_WEIGHT() * Math.abs(primer.gcContent - probe.getGcContent()) / Primer.DELTA_GC_MAX
				+ searchParams.getSELF_ALIGNMENT_WEIGHT() * Math.max(primer.selfAlignmentScore, probe.getSelfAlignmentScore()) / Primer.DELTA_PRIMER_SA_MAX
				+ searchParams.getSELF_END_ALIGNMENT_WEIGHT() * Math.max(probe.selfEndAlignmentScore, probe.getSelfEndAlignmentScore()) / Primer.DELTA_PRIMER_SEA_MAX
				+ searchParams.getPAIR_ALIGNMENT_WEIGHT() * scores.getPairScore() / Primer.DELTA_PRIMER_PA_MAX
				+ searchParams.getPAIR_END_ALIGNMENT_WEIGHT() * scores.getPairEndScore() / Primer.DELTA_PRIMER_PEA_MAX;
		
		return score;
	}
	
	/**
	 * Returns a textual representation of this primer.
	 * 
	 * @return a textual representation of this primer
	 */
	public String toString(){
		return this.relativePosition + "\t" + this.sequence + "\t" + this.length + "\t" + this.distanceToRSS + "\t" + this.gcContent + "\t" + this.meltingTemp + "\t" + this.distanceToOptimalPrimer + "\t" + this.selfAlignmentScore + "\t" + this.selfEndAlignmentScore;
	}
	
	public String toFormattedString(){
		return String.format("%-10s\t%-40s\t%s\t%s\t%.3f\t%.2f\t%.6f\t%s\t%s", new Object[]{this.relativePosition, this.sequence, this.length, this.distanceToRSS, this.gcContent, this.meltingTemp, this.distanceToOptimalPrimer, this.selfAlignmentScore, this.selfEndAlignmentScore});
	}
	
	/**
	 * Returns the format of the primer informations return by toString.
	 * 
	 * This format is supposed to be in sync with the 'toString' method!
	 * 
	 * @return the format of the primer informations return by toString
	 */
	public static String toStringDescription(){
		return "Position\tSequence\tlength\tdistRSS\t%GC\tTM\tdOptPrimer\tscoreSA\tscoreSEA";
	}
	
	public static String toFormattedStringDescription(){
		return String.format("%-10s\t%-40s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", new Object[]{"Position", "Sequence", "length", "distRSS", "%GC", "TM", "dOptPair", "scoreSA", "scoreSEA"});
	}
	
	/**
	 * Returns the virtual optimal primer for primer scoring.
	 * 
	 * @return the virtual optimal primer for primer scoring
	 */
	public static Primer getVirtualOptimalPrimer(PrimerSearchParameters searchParams){
		Primer primer = new Primer(searchParams);
		
		primer.primerType = PrimerTypes.forwardPrimer;
		primer.distanceToOptimalPrimer = 0;
		primer.distanceToRSS = searchParams.getOPT_DISTANCE_TO_RSS();
		primer.falsePositiveMatches = 0;
		primer.gcContent = searchParams.getOPT_GC();
		primer.length = searchParams.getOPT_PRIMER_LENGTH();
		primer.meltingTemp = searchParams.getOPT_TM();
		primer.relativePosition = 0;
		primer.selfAlignmentScore = 0;
		primer.selfEndAlignmentScore = 0;
		primer.sequence = "";
		primer.isVirtualOptimalPrimer = true;
		
		return primer;
	}
	
	/**
	 * Returns the virtual optimal probe for probe scoring.
	 * 
	 * @return the virtual optimal probe for probe scoring
	 */
	public static Primer getVirtualOptimalProbe(PrimerSearchParameters searchParams){
		Primer primer = Primer.getVirtualOptimalPrimer(searchParams);
		
		primer.primerType = PrimerTypes.hybridizationProbe;
		primer.gcContent = searchParams.getTAQMAN_OPT_GC();
		primer.length = searchParams.getTAQMAN_OPT_PRIMER_LENGTH();
		primer.meltingTemp = searchParams.getTAQMAN_OPT_TM();
		primer.isVirtualOptimalPrimer = true;
		
		return primer;
	}
	
	/**
	 * Returns true iff this primer is the virtual optimal primer.
	 * 
	 * @return true iff this primer is the virtual optimal primer
	 */
	public boolean isVirtualOptimalPrimer(){
		return this.isVirtualOptimalPrimer;
	}
	
	public static void main(String[] args){
		String[] primers = new String[]{"TCGACATCTCATTTGAGGC", "GGAAAGAAACAGTTCGAATGG", "AACGTTCCCGTCAGCAAC", "TGTAATCACATGAAGAAGTACTTGG", "AATGGTTTCAGAACCGCA", "GGTGAGAACTTGGGGACTTTA", "TGCTTGAGTGAATCGTTATAGTTAC", "CGAAAAGGCGTCGAATGTA", "AATACCGCGGCACAAATT"};
		PrimerSearchParameters searchParams = new PrimerSearchParameters();
		
		for(String sequence : primers){
			Primer primer = new Primer(sequence, PrimerTypes.forwardPrimer, searchParams);
			System.out.println(primer.toString());
		}
	}

	/**
	 * Returns the restrictionSite.
	 * 
	 * @return the restrictionSite
	 */
	public RestrictionSite getRestrictionSite() {
		return restrictionSite;
	}

	/**
	 * Sets the restrictionSite.
	 * 
	 * @param restrictionSite the restrictionSite to set
	 */
	public void setRestrictionSite(RestrictionSite restrictionSite) {
		this.restrictionSite = restrictionSite;
	}

	/**
	 * @return the positionInScanSequence
	 */
	public int getPositionInScanSequence() {
		return positionInScanSequence;
	}

	/**
	 * @param positionInScanSequence the positionInScanSequence to set
	 */
	public void setPositionInScanSequence(int positionInScanSequence) {
		if(positionInScanSequence < 0) throw new IllegalArgumentException("Position must be >= 0!");
		this.positionInScanSequence = positionInScanSequence;
	}
	
	public synchronized PrimerAcceptanceLevel getAcceptanceLevel(){
		return this.acceptanceLevel;
	}
	
	public synchronized void setAcceptanceLevel(PrimerAcceptanceLevel level){
		this.acceptanceLevel = level;
	}
}
