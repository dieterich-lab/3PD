package primerDesign.dsc;

import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.util.PrimerSearchParameters;

/**
 * Encapsulaes a set of: forward, reverse and hybridization probe primers.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPair implements Comparable{
	private Primer forwardPrimer;
	private Primer reversePrimer;
	private Primer hybridizationProbe;
	private PrimerAlignmentScores maxAlignmentScore;
	private double averageTMprimers;
	private double averageTMprobe;
	private double distanceToOptimalPrimerPair;
	private PrimerSearchParameters searchParams;
	private boolean wasModified;
	public static final double epsilon = 1E-12;
	private PrimerAlignmentCalculation alignment;
	private PrimerAcceptanceLevel acceptancelevel;

	/**
	 * Initializes a primer pair.
	 * 
	 * @param searchParams the 3PD primer search parameters
	 */
	public PrimerPair(PrimerSearchParameters searchParams){
		this.distanceToOptimalPrimerPair = -1;
		alignment = searchParams.getPRIMER_ALIGNMENT_METHOD();
		this.searchParams = searchParams;
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	/**
	 * Initializes a new primer set with the respective primers (s.b.).
	 * 
	 * @param forward the forward primer of this set
	 * @param reverse the reverse primer of this set
	 * @param probe the hybridization probe of this set
	 */
	public PrimerPair(Primer forward, Primer reverse, Primer probe, PrimerSearchParameters searchParams){
		this.forwardPrimer = forward;
		this.reversePrimer = reverse;
		this.hybridizationProbe = probe;
		this.distanceToOptimalPrimerPair = -1;
		this.averageTMprimers = forward.getMeltingTemp() + reverse.getMeltingTemp();
		this.averageTMprobe = probe.getMeltingTemp();
		this.searchParams = searchParams;
		alignment = searchParams.getPRIMER_ALIGNMENT_METHOD();
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;	
	}
	
	/**
	 * Initializes a new primer set with the respective primers (s.b.).
	 * 
	 * @param forward the forward primer of this set
	 * @param reverse the reverse primer of this set
	 * @param probe the hybridization probe of this set
	 * @param fwRev forward-reverse alignment scores
	 * @param fwProbe forward-probe alignment scores
	 * @param revProbe reverse-probe alignment scores
	 */
	public PrimerPair(Primer forward, Primer reverse, Primer probe, PrimerAlignmentScores fwRev, PrimerAlignmentScores fwProbe, PrimerAlignmentScores revProbe, PrimerSearchParameters searchParams){
		this.forwardPrimer = forward;
		this.reversePrimer = reverse;
		this.hybridizationProbe = probe;
		this.distanceToOptimalPrimerPair = -1;
		this.averageTMprimers = forward.getMeltingTemp() + reverse.getMeltingTemp();
		if(probe != null) this.averageTMprobe = probe.getMeltingTemp();
		else this.averageTMprobe = searchParams.getTAQMAN_OPT_TM();
		this.searchParams = searchParams;
		this.alignment = null;
		this.maxAlignmentScore = PrimerAlignmentScores.getMaxScores(new PrimerAlignmentScores[]{fwRev, fwProbe, revProbe});
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	/**
	 * Initializes a new primer set with the respective primers (s.b.).
	 * 
	 * @param forward the forward primer of this set
	 * @param reverse the reverse primer of this set
	 * @param probe the hybridization probe of this set
	 * @param fwRev forward-reverse alignment scores
	 * @param fwProbe forward-probe alignment scores
	 * @param revProbe reverse-probe alignment scores
	 */
	public PrimerPair(Primer forward, Primer reverse, PrimerAlignmentScores fwRev, PrimerSearchParameters searchParams){
		this.forwardPrimer = forward;
		this.reversePrimer = reverse;
		this.hybridizationProbe = null;
		this.distanceToOptimalPrimerPair = -1;
		this.averageTMprimers = forward.getMeltingTemp() + reverse.getMeltingTemp();
		if(this.hybridizationProbe != null) this.averageTMprobe = this.hybridizationProbe.getMeltingTemp();
		else this.averageTMprobe = searchParams.getTAQMAN_OPT_TM();
		this.searchParams = searchParams;
		this.alignment = null;
		this.maxAlignmentScore = fwRev; 
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	/**
	 * Copy constructor creating a shallow copy of an object.
	 * 
	 * This constructor returns a new PrimerPair instance with symbolic references to the single primers of the original pair!
	 * 
	 * @param other the other PrimerPair to generate the copy from
	 */
	public PrimerPair(PrimerPair other, PrimerSearchParameters searchParams){
		this.forwardPrimer = other.getForwardPrimer();
		this.reversePrimer = other.getReversePrimer();
		this.hybridizationProbe = other.getHybridizationProbe();
		this.averageTMprimers = other.averageTMprimers;
		this.averageTMprobe = other.averageTMprobe;
		this.searchParams = searchParams;
		alignment = other.alignment;
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	/**
	 * Returns the forwardPrimer.
	 * 
	 * @return the forwardPrimer
	 */
	public Primer getForwardPrimer() {
		return forwardPrimer;
	}

	/**
	 * Sets the forwardPrimer.
	 * 
	 * @param forwardPrimer the forwardPrimer to set
	 */
	public void setForwardPrimer(Primer forwardPrimer) {
		if(this.forwardPrimer != null) this.averageTMprimers -= this.forwardPrimer.getMeltingTemp();
		this.averageTMprimers += forwardPrimer.getMeltingTemp();
		this.forwardPrimer = forwardPrimer;
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}

	/**
	 * Returns the hybridizationProbe.
	 * 
	 * @return the hybridizationProbe
	 */
	public Primer getHybridizationProbe() {
		return hybridizationProbe;
	}

	/**
	 * Sets the hybridizationProbe.
	 * 
	 * @param hybridizationProbe the hybridizationProbe to set
	 */
	public void setHybridizationProbe(Primer hybridizationProbe) {
		if(hybridizationProbe != null) this.averageTMprobe = hybridizationProbe.getMeltingTemp();
		else this.averageTMprobe = searchParams.getTAQMAN_OPT_TM();
		this.hybridizationProbe = hybridizationProbe;
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}

	/**
	 * Returns the reversePrimer.
	 * 
	 * @return the reversePrimer
	 */
	public Primer getReversePrimer() {
		return reversePrimer;
	}

	/**
	 * Returns the reversePrimer.
	 * 
	 * @param reversePrimer the reversePrimer to set
	 */
	public void setReversePrimer(Primer reversePrimer) {
		if(this.reversePrimer != null) this.averageTMprimers -= this.reversePrimer.getMeltingTemp();
		this.averageTMprimers += reversePrimer.getMeltingTemp();
		this.reversePrimer = reversePrimer;
		this.wasModified = true;
		this.acceptancelevel = PrimerAcceptanceLevel.NOT_TESTED;
	}
	
	/**
	 * Computes the distance to the optimal primer pair.
	 */
	public void computeDistanceToOptimalPrimerPair(){
		this.distanceToOptimalPrimerPair = this.scoreTo(PrimerPair.getOptimalPrimerPair(this.searchParams));
	}
	
	/**
	 * Returns the distance to the optimal primer pair.
	 * 
	 * @return the distance to the optimal primer pair
	 */
	public double getDistanceToOptimalPrimerPair(){
		if(this.wasModified || this.distanceToOptimalPrimerPair == -1){
			this.computeDistanceToOptimalPrimerPair();
			if(this.alignment != null) this.computeMaxAlignmentScore();
			this.wasModified = false;
		}
		return this.distanceToOptimalPrimerPair;
	}
	
	/**
	 * Scores this primer pair against another one.
	 * 
	 * @param other the other primer pair to score this primer pair to
	 * 
	 * @return the score of this primer pair when compared to another primer pair
	 */
	public double scoreTo(PrimerPair other){
		double score = 0;
		
		score += this.getForwardPrimer().scoreTo(other.getForwardPrimer(), this.searchParams);
		score += this.getForwardPrimer().scoreTo(other.getReversePrimer(), this.searchParams);
		score += this.getForwardPrimer().scoreTo(other.getHybridizationProbe(), this.searchParams);
		
		score += this.getReversePrimer().scoreTo(other.getForwardPrimer(),this.searchParams);
		score += this.getReversePrimer().scoreTo(other.getReversePrimer(), this.searchParams);
		score += this.getReversePrimer().scoreTo(other.getHybridizationProbe(), this.searchParams);
		
		if(searchParams.isPickTaqManProbe()){
			score += this.getHybridizationProbe().scoreTo(other.getForwardPrimer(), this.searchParams);
			score += this.getHybridizationProbe().scoreTo(other.getReversePrimer(), this.searchParams);
			score += this.getHybridizationProbe().scoreTo(other.getHybridizationProbe(), this.searchParams);
			
			return score/9;
		}
		else return score/6;
	}
	
	/**
	 * Returns the alignment values for the primer pair.
	 * 
	 * @param pair the primer pair to return alignemnt values for
	 * 
	 * @return the alignment values for the primer pair
	 */
	public PrimerAlignmentScores getAlignmentValues(PrimerPair pair) {
		int sa_max = -1;
		int sea_max = -1;
		PrimerAlignmentScores scores;
		
		// check this pairs's score
		scores = this.getMaxPairElementsAlignmentScore();
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore();
		
		// check other pairs's score
		scores = pair.getMaxPairElementsAlignmentScore();
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore();
		
		// check combination of primers's score
		// TaqMan probe-TaqMan probe pairs need not to be aligned since in each reaction well exactly one of the TaqMan probes is contained!
		scores = alignment.computePairAlignment(this.forwardPrimer.getSequence(),pair.getForwardPrimer().getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
//		scores = alignment.computePairAlignment(this.forwardPrimer.getSequence(),pair.getReversePrimer().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
		scores = alignment.computePairAlignment(this.forwardPrimer.getSequence(),pair.getHybridizationProbe().getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
		
//		scores = alignment.computePairAlignment(this.reversePrimer.getSequence(),pair.getForwardPrimer().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
//		scores = alignment.computePairAlignment(this.reversePrimer.getSequence(),pair.getReversePrimer().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
//		scores = alignment.computePairAlignment(this.reversePrimer.getSequence(),pair.getHybridizationProbe().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
		
		scores = alignment.computePairAlignment(this.hybridizationProbe.getSequence(),pair.getForwardPrimer().getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
//		scores = alignment.computePairAlignment(this.hybridizationProbe.getSequence(),pair.getReversePrimer().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max > scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
//		scores = alignment.computePairAlignment(this.hybridizationProbe.getSequence(),pair.getHybridizationProbe().getSequence());
//		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
//		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore(); 
		
		return new PrimerAlignmentScores(sa_max,sea_max);
	}
	
	/**
	 * Computes alignment scores.
	 */
	public void computeMaxAlignmentScore(){
		PrimerAlignmentScores scores;
		int sa_max = -1;
		int sea_max = -1;
		
		scores = alignment.computePairAlignment(this.forwardPrimer.getSequence(),this.reversePrimer.getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore();
		
		scores = alignment.computePairAlignment(this.forwardPrimer.getSequence(),this.hybridizationProbe.getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore();
		
		scores = alignment.computePairAlignment(this.reversePrimer.getSequence(),this.hybridizationProbe.getSequence());
		if(sa_max < scores.getPairScore()) sa_max = scores.getPairScore();
		if(sea_max < scores.getPairEndScore()) sea_max = scores.getPairEndScore();
		
		this.maxAlignmentScore = new PrimerAlignmentScores(sa_max,sea_max);
	}
	
	/**
	 * Returns the maximum pair alignment score.
	 * 
	 * @return the maximum pair alignment score
	 */
	public PrimerAlignmentScores getMaxPairElementsAlignmentScore(){
		if(this.wasModified){
			if(this.alignment != null) this.computeMaxAlignmentScore();
			this.computeDistanceToOptimalPrimerPair();
			this.wasModified = false;
		}
		
		return this.maxAlignmentScore;
	}
	
	/**
	 * Returns the virtual optimal primer pair.
	 * 
	 * @param searchParams the 3PD search parameters
	 * @return the virtual optimal primer pair
	 */
	public static PrimerPair getOptimalPrimerPair(PrimerSearchParameters searchParams){
		return new PrimerPair(Primer.getVirtualOptimalPrimer(searchParams), Primer.getVirtualOptimalPrimer(searchParams), Primer.getVirtualOptimalProbe(searchParams), searchParams);
	}
	
	/**
	 * Returns the average primer melting temperature.
	 * 
	 * @return the average primer melting temperature
	 */
	public double getAverageTMPrimers(){
		return this.averageTMprimers / 2;
	}
	
	/**
	 * Returns the average probe melting temperature.
	 * 
	 * @return the average probe melting temperature
	 */
	public double getAverageTMProbe(){
		return this.averageTMprobe;
	}
	
	/**
	 * Returns true iff this primer is the virtual optimal primer.
	 * 
	 * @return true iff this primer is the virtual optimal primer
	 */
	public boolean isVirtualOptimalPrimerPair(){
		return this.forwardPrimer.isVirtualOptimalPrimer() && this.reversePrimer.isVirtualOptimalPrimer() && this.hybridizationProbe.isVirtualOptimalPrimer();
	}
	
	/**
	 * Returns true iff this primer pair contains unacceptable primers.
	 * 
	 * @return true iff this primer pair contains unacceptable primers
	 */
	public synchronized boolean containsUnacceptablePrimers(){
		if(this.forwardPrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.NOT_TESTED)){		
			if(this.searchParams.getPrimerMisprimingCheck().hasMisprimings(this.forwardPrimer)){
				this.forwardPrimer.setAcceptanceLevel(PrimerAcceptanceLevel.UNACCEPTABLE);
			}
			else{
				this.forwardPrimer.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
			}
		}
		if(this.reversePrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.NOT_TESTED)){
			if(this.searchParams.getPrimerMisprimingCheck().hasMisprimings(this.reversePrimer)){
				this.reversePrimer.setAcceptanceLevel(PrimerAcceptanceLevel.UNACCEPTABLE);
			}
			else{
				this.reversePrimer.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
			}
		}
		if(this.acceptancelevel.equals(PrimerAcceptanceLevel.NOT_TESTED)){
			if(this.forwardPrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE) || this.reversePrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE)){
				this.acceptancelevel = PrimerAcceptanceLevel.UNACCEPTABLE;
			}
			else if(this.forwardPrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.ACCEPTABLE) && this.reversePrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.ACCEPTABLE)){
				this.acceptancelevel = PrimerAcceptanceLevel.ACCEPTABLE;
			}
			else{
				throw new IllegalStateException("Illegal state!");
			}
		}
		return this.forwardPrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE) ||
		this.reversePrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE);
	}
	
	/**
	 * Returns true iff this primer set was specified being NOT acceptable.
	 * Returns false otherwise (not tested or acceptable!).
	 * 
	 * @return true iff this primer set was specified being NOT acceptable
	 */
	public synchronized boolean isNotAcceptablePair(){
		return this.forwardPrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE) || this.reversePrimer.getAcceptanceLevel().equals(PrimerAcceptanceLevel.UNACCEPTABLE);
	}
	
	/**
	 * Returns the primer acceptance level.
	 * 
	 * @return the primer acceptance level
	 */
	public synchronized PrimerAcceptanceLevel getAcceptanceLevel(){
		return this.acceptancelevel;
	}
	
	/**
	 * Sets the primer acceptance level.
	 * 
	 * @param level the primer acceptance level to set
	 */
	public synchronized void setAcceptanceLevel(PrimerAcceptanceLevel level){
		this.acceptancelevel = level;
	}

	/** 
	 * Returns a textual representation of this primer pair.
	 * 
	 * @return a textual representation of this primer pair
	 */
	public String toString(){
		if(searchParams.isPickTaqManProbe()){
			return "Forward:\t" + forwardPrimer.toString() + "\nProbe:\t" + hybridizationProbe.toString() + "\nReverse:\t" + reversePrimer.toString() + "\ndOptPair: " + this.getDistanceToOptimalPrimerPair() + "\nPAmax: " + maxAlignmentScore.getPairScore() + " PEA: " + maxAlignmentScore.getPairEndScore();
		}
		else{
			return "Forward:\t" + forwardPrimer.toString() + "\nReverse:\t" + reversePrimer.toString() + "\ndOptPair: " + this.getDistanceToOptimalPrimerPair() + "\nPAmax: " + maxAlignmentScore.getPairScore() + " PEA: " + maxAlignmentScore.getPairEndScore() ;
		}
	}
	
	/**
	 * Returns a string of this primer pair formatted in FASTA format.
	 * ,
	 * @return a string of this primer pair formatted in FASTA format
	 */
	public String toFormattedString(){
		if(searchParams.isPickTaqManProbe()){
			return String.format("%10s\t%s%s", "Forward:", forwardPrimer.toFormattedString(), "\n") + String.format("%10s\t%s%s", "Probe:", hybridizationProbe.toFormattedString(), "\n") + String.format("%10s\t%s%s", "Reverse:", reversePrimer.toFormattedString(), "\n") + "dOptPair: " + this.getDistanceToOptimalPrimerPair() + "\nPAmax: " + maxAlignmentScore.getPairScore() + " PEA: " + maxAlignmentScore.getPairEndScore() +"\nHashCode: " + this.getHashCode();
		}
		else{
			return String.format("%10s\t%s%s", "Forward:", forwardPrimer.toFormattedString(), "\n") + String.format("%10s\t%s%s", "Reverse:", reversePrimer.toFormattedString(), "\n") + "dOptPair: " + this.getDistanceToOptimalPrimerPair() + "\nPAmax: " + maxAlignmentScore.getPairScore() + " PEA: " + maxAlignmentScore.getPairEndScore();
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {
		PrimerPair other = (PrimerPair) o;
		if(this.getDistanceToOptimalPrimerPair() < other.getDistanceToOptimalPrimerPair()) return -1;
		else if(this.getDistanceToOptimalPrimerPair() > other.getDistanceToOptimalPrimerPair()) return 1;
		return 0;
	}
	
	public boolean equals(PrimerPair other){
		return this.forwardPrimer.getSequence().equals(other.getForwardPrimer().getSequence()) &&
				this.reversePrimer.getSequence().equals(other.getReversePrimer().getSequence()) &&
				this.hybridizationProbe.getSequence().equals(other.getHybridizationProbe().getSequence());
	}
	
	public int getHashCode(){
		if(this.searchParams.isPickTaqManProbe()) return this.forwardPrimer.getSequence().hashCode() + this.reversePrimer.getSequence().hashCode() + this.hybridizationProbe.getSequence().hashCode();
		else return this.forwardPrimer.getSequence().hashCode() + this.reversePrimer.getSequence().hashCode();
	}
}
