package primerDesign.dsc;

import org.biojava.bio.molbio.RestrictionEnzyme;

import primerDesign.algo.SequenceRegionAligner;
import primerDesign.util.Constants;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.PrimerSearchParameters;
import cern.colt.list.ObjectArrayList;

/**
 * This class implements a restriction site on a piece of dna.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RestrictionSite implements Comparable{

	private RestrictionEnzyme enzyme;
	private int position;
	private Primer[] validUpstreamPrimers;
	private Primer[] validDownstreamPrimers;
	private Primer[] validTaqManProbes;
	private ObjectArrayList primerPairs;
	private SequenceRegion sequenceRegion;
	private int distanceToIntervalMean;
	private PrimerSearchParameters searchParams;
	private boolean wasScannedForPrimers = false;
	private boolean hasEnumeratedPrimerPairs = false;
	private char[] forwardScanSequence; // the forward primer scan sequence in 5'->3' direction
	private char[] reverseScanSequence; // the reverse primer scan sequence in 5'->3' direction
	private char[] probeScanSequence; // the probe primer scan sequence in 5'->3' direction
	
	/** 
	 * Empty default constructor.
	 * 
	 * @param searchParams the 3PD search parameters
	 */
	public RestrictionSite(PrimerSearchParameters searchParams){
		this.searchParams = searchParams;
	}
	
	/**
	 * Parameterized constructur.
	 * 
	 * @param position the position of the restriction site
	 * @param enzyme the enzyme generating this restriction site
	 * @param searchParams the 3PD search parameters
	 */
	public RestrictionSite(int position, RestrictionEnzyme enzyme, PrimerSearchParameters searchParams){
		
		if(position < 0) throw new IllegalArgumentException("Position must be >= 0!");
		if(enzyme.getRecognitionSite().length() == 0) throw new IllegalArgumentException("Restriction recognition sequence of restriction enzyme is not set!");
		
		this.setPosition(position);
		this.enzyme = enzyme;
		this.validUpstreamPrimers = new Primer[0];
		this.validDownstreamPrimers = new Primer[0];
		this.validTaqManProbes = new Primer[0];
		this.primerPairs = new ObjectArrayList();
		this.sequenceRegion = new SequenceRegion(null, 0,0);
		this.distanceToIntervalMean = Integer.MIN_VALUE;
		this.searchParams = searchParams;
	}
	
	/**
	 * Parameterized constructur.
	 * 
	 * @param position the position of the restriction site
	 * @param enzyme the enzyme generating this restriction site
	 * @param searchParams the 3PD search parameters
	 */
	public RestrictionSite(SequenceRegion region, int position, RestrictionEnzyme enzyme, PrimerSearchParameters searchParams){
		
		if(position < 0) throw new IllegalArgumentException("Position must be >= 0!");
		if(enzyme.getRecognitionSite().length() == 0) throw new IllegalArgumentException("Restriction recognition sequence of restriction enzyme is not set!");
		
		this.setPosition(position);
		this.enzyme = enzyme;
		this.validUpstreamPrimers = new Primer[0];
		this.validDownstreamPrimers = new Primer[0];
		this.validTaqManProbes = new Primer[0];
		this.primerPairs = new ObjectArrayList();
		this.sequenceRegion = region;
		this.distanceToIntervalMean = Integer.MIN_VALUE;
		this.searchParams = searchParams;
	}
	
	public boolean wasScannedForPrimers(){
		return this.wasScannedForPrimers;
	}
	
	public boolean hasEnumeratedPrimerPairs(){
		return this.hasEnumeratedPrimerPairs;
	}
	
	/**
	 * Returns all valid upstream primers of a restriction site.
	 * 
	 * @return the validUpstreamPrimers of a restriction site
	 */
	public Primer[] getValidUpstreamPrimers() {
		if(validUpstreamPrimers.length == 0) throw new EmptyResultSetException("No upstream-primers were set for this restriction site.");
		return validUpstreamPrimers;
	}
	
	/**
	 * Returns all valid downstream primers of a restriction site.
	 * 
	 * @return the validDownstreamPrimers of a restriction site
	 */
	public Primer[] getValidDownstreamPrimers(){
		if(validDownstreamPrimers.length == 0) throw new EmptyResultSetException("No downstream-primers were set for this restriction site.");
		return validDownstreamPrimers;
	}
	
	/**
	 * Returns all valid TaqMan probes of a restriction site.
	 * 
	 * @return the validTaqManProbes of a restriction site
	 */
	public Primer[] getValidTaqManProbes(){
		if(validTaqManProbes.length == 0) throw new EmptyResultSetException("No TaqMan probes were set for this restriction site.");
		return validTaqManProbes;
	}
	
	/**
	 * Returns the number of valid upstream primers at this restriction site.
	 * 
	 * @return the number of valid upstream primers at this restriction site
	 */
	public int getNumberOfValidUpstreamPrimers(){
		return this.validUpstreamPrimers.length;
	}
	
	/**
	 * Returns the number of valid downstream primers at this restriction site.
	 * 
	 * @return the number of valid downstream primers at this restriction site
	 */
	public int getNumberOfValidDownstreamPrimers(){
		return this.validDownstreamPrimers.length;
	}
	
	/**
	 * Returns the number of valid TaqMan probes at this restriction site.
	 * 
	 * @return the number of valid TaqMan probes at this restriction site
	 */
	public int getNumberOfValidTaqManProbes(){
		return this.validTaqManProbes.length;
	}
	
	/**
	 * Returns the upstream primer at position 'index'.
	 * 
	 * @param index index the index of the upstream primer to be returned
	 * 
	 * @return the upstream primer at position 'index'
	 */
	public Primer getUpstreamPrimer(int index){
		if(index < this.validUpstreamPrimers.length) return (Primer) this.validUpstreamPrimers[index];
		else throw new IllegalArgumentException("Invalid primer index!");
	}
	
	/**
	 * Returns the downstream primer at position 'index'.
	 * 
	 * @param index index the index of the downstream primer to be returned
	 * 
	 * @return the downstream primer at position 'index'
	 */
	public Primer getDownstreamPrimer(int index){
		if(index < this.validDownstreamPrimers.length) return (Primer) this.validDownstreamPrimers[index];
		else throw new IllegalArgumentException("Invalid primer index!");
	}
	
	/**
	 * Returns the TaqMan probe at position 'index'.
	 * 
	 * @param index index the index of the TaqMan probe to be returned
	 * 
	 * @return the TaqMan probe at position 'index'
	 */
	public Primer getTaqManProbe(int index){
		if(index < this.validTaqManProbes.length) return (Primer) this.validTaqManProbes[index];
		else throw new IllegalArgumentException("Invalid primer index!");
	}
	
	/**
	 * Sets the valid upstream primers of a restriction site.
	 * 
	 * @param primers a vector of valid upstream primers
	 */
	public void setValidUpstreamPrimers(Primer[] primers){
		if(primers.length == 0) throw new IllegalArgumentException("The upstream-primer vector 'primers' is empty! - no primers can be set!");
		this.validUpstreamPrimers = primers;
//		if(this.validUpstreamPrimers.length > 0 && this.validDownstreamPrimers.length > 0 && this.validTaqManProbes.length > 0){
//			this.wasScannedForPrimers = true;
//		}
		this.wasScannedForPrimers = true;
	}
	
	/**
	 * Sets the valid downstream primers of a restriction site.
	 * 
	 * @param primers a vector of valid downstream primers
	 */
	public void setValidDownstreamPrimers(Primer[] primers){
		if(primers.length == 0) throw new IllegalArgumentException("The downstream-primer vector 'primers' is empty! - no primers can be set!");
		this.validDownstreamPrimers = primers;
//		if(this.validUpstreamPrimers.length > 0 && this.validDownstreamPrimers.length > 0 && this.validTaqManProbes.length > 0){
//			this.wasScannedForPrimers = true;
//		}
		this.wasScannedForPrimers = true;
	}
	
	/**
	 * Sets the valid TaqMan probes of a restriction site.
	 * 
	 * @param primers a vector of valid TaqMan probes
	 */
	public void setValidTaqManProbes(Primer[] primers){
		if(primers.length == 0) throw new IllegalArgumentException("The TaqMan probe vector 'primers' is empty! - no primers can be set!");
		this.validTaqManProbes = primers;
//		if(this.validUpstreamPrimers.length > 0 && this.validDownstreamPrimers.length > 0 && this.validTaqManProbes.length > 0){
//			this.wasScannedForPrimers = true;
//		}
		this.wasScannedForPrimers = true;
	}
	
	public void wasScannedForPrimers(boolean value){
		if(value == true && this.validUpstreamPrimers.length > 0 && this.validDownstreamPrimers.length > 0 && (!searchParams.isPickTaqManProbe() || this.validTaqManProbes.length > 0)) this.wasScannedForPrimers = value;
		else if(value == false) this.wasScannedForPrimers = value;
		else throw new IllegalStateException("Not every primer list contains at least one entry!");
	}

	/**
	 * Returns the restriction enzyme which generated this restriction site.
	 * 
	 * @return the restriction enzyme which generated this restriction site
	 */
	public RestrictionEnzyme getEnzyme() {
		return enzyme;
	}

	/**
	 * Sets the enzyme which generated thius restriction site.
	 * 
	 * @param enzyme the enzyme which generated thius restriction site
	 */
	public void setEnzyme(RestrictionEnzyme enzyme) {
		if(enzyme.getRecognitionSite().length() == 0) throw new IllegalArgumentException("Restriction recognition sequence of restriction enzyme is not set!");
		this.enzyme = enzyme;
	}

	/**
	 * Returns the position of the restriction site.
	 * 
	 * @return the position of the restriction site
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Sets the position of the retriction site.
	 * 
	 * @param position the position of the retriction site
	 */
	public void setPosition(int position) {
		if(position < 0) throw new IllegalArgumentException("Position must be >= 0");
		this.position = position;
	}

	/**
	 * Returns the distance to the mean of the interval, this restriction site is located on.
	 * 
	 * @return the distance to the mean of the interval, this restriction site is located on
	 */
	public int getDistanceToIntervalMean() {
		return distanceToIntervalMean;
	}

	/**
	 * Sets the distance to the mean of the interval, this restriction site is located on
	 * 
	 * @param distanceToIntervalMean the distance to the mean of the interval, this restriction site is located on
	 */
	public void setDistanceToIntervalMean(int distanceToIntervalMean) {
		if(distanceToIntervalMean < 0) throw new IllegalArgumentException("Distance to interval mean must be >= 0!");
		this.distanceToIntervalMean = distanceToIntervalMean;
	}

	/**
	 * Returns the sequence region, this restriction site is located on.
	 * 
	 * @return the sequence region, this restriction site is located on.
	 */
	public SequenceRegion getSequenceRegion() {
		return sequenceRegion;
	}

	/**
	 * Sets the sequence region, this restriction site is located on.
	 * 
	 * @param sequenceRegion the sequence region, this restriction site is located on.
	 */
	public void setSequenceRegion(SequenceRegion sequenceRegion) {
		if(!(sequenceRegion.getSeqRegionStart() >= 0 && sequenceRegion.getSeqRegionEnd() >= 0 && sequenceRegion.getSeqRegionEnd() >= sequenceRegion.getSeqRegionStart())) throw new IllegalArgumentException("sequenceRegion was not properly initialized!");
		this.sequenceRegion = sequenceRegion;
	}
	
	/**
	 * This method enumerates all valid primer pairs each consisting of one upstream, one downstream and one hybridization probe primer.
	 * 
	 * After enumeration, the list of valid primer pairs is orderd by increasing score w.r.t the virtual optimal primer.
	 *
	 */
	public void enumeratePrimerPairs(){
		if(this.validUpstreamPrimers.length == 0 || this.validDownstreamPrimers.length == 0 || (searchParams.isPickTaqManProbe() && this.validTaqManProbes.length == 0)) throw new IllegalStateException("At least one list of valid primers is empty!");
		
		PrimerPair currentPair;
		SequenceRegionAlignment fwRev = SequenceRegionAligner.alignSequenceRegions(this.getForwardScanSequence(), this.getReverseScanSequence(), this.searchParams.getA_t_basepair_score(), this.searchParams.getG_c_basepair_score());
		SequenceRegionAlignment fwProbe = null;
		SequenceRegionAlignment revProbe = null;
		if(searchParams.isPickTaqManProbe()){
			fwProbe = SequenceRegionAligner.alignSequenceRegions(this.getForwardScanSequence(), this.getProbeScanSequence(), this.searchParams.getA_t_basepair_score(), this.searchParams.getG_c_basepair_score());
			revProbe = SequenceRegionAligner.alignSequenceRegions(this.getReverseScanSequence(), this.getProbeScanSequence(), this.searchParams.getA_t_basepair_score(), this.searchParams.getG_c_basepair_score());
		}
		for(int i=0; i<this.validUpstreamPrimers.length; i++){
			for(int j=0; j<this.validDownstreamPrimers.length; j++){
				if(searchParams.isPickTaqManProbe()){
					for(int k=0; k<this.validTaqManProbes.length; k++){
						currentPair = new PrimerPair(this.validUpstreamPrimers[i], this.validDownstreamPrimers[j], this.validTaqManProbes[k], fwRev.getGlobalAlignmentValues(this.validUpstreamPrimers[i].getPositionInScanSequence(), this.validUpstreamPrimers[i].getLength(), this.validDownstreamPrimers[j].getPositionInScanSequence(), this.validDownstreamPrimers[j].getLength()), fwProbe.getGlobalAlignmentValues(this.validUpstreamPrimers[i].getPositionInScanSequence(), this.validUpstreamPrimers[i].getLength(), this.validTaqManProbes[k].getPositionInScanSequence(), this.validTaqManProbes[k].getLength()), revProbe.getGlobalAlignmentValues(this.validDownstreamPrimers[j].getPositionInScanSequence(), this.validDownstreamPrimers[j].getLength(), this.validTaqManProbes[k].getPositionInScanSequence(), this.validTaqManProbes[k].getLength()), this.searchParams);
						currentPair.computeDistanceToOptimalPrimerPair();
						if(isValidPrimerPair(currentPair, fwRev, fwProbe, revProbe)){
							if(Constants.doEarlyMMScan) currentPair.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
							this.primerPairs.add(currentPair);
							this.searchParams.getPickingStat().incValidPrimerPairsEnum();
						}
						else this.searchParams.getPickingStat().incInvalidPrimerPairsEnum();
					}
				}
				else{
					currentPair = new PrimerPair(this.validUpstreamPrimers[i], this.validDownstreamPrimers[j], fwRev.getGlobalAlignmentValues(this.validUpstreamPrimers[i].getPositionInScanSequence(), this.validUpstreamPrimers[i].getLength(), this.validDownstreamPrimers[j].getPositionInScanSequence(), this.validDownstreamPrimers[j].getLength()), this.searchParams);
					currentPair.computeDistanceToOptimalPrimerPair();
					if(isValidPrimerPair(currentPair, fwRev, fwProbe, revProbe)){
						if(Constants.doEarlyMMScan) currentPair.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
						this.primerPairs.add(currentPair);
						this.searchParams.getPickingStat().incValidPrimerPairsEnum();
					}
					else this.searchParams.getPickingStat().incInvalidPrimerPairsEnum();
				}
			}
		}
		if(this.primerPairs.size() > 0){
			this.hasEnumeratedPrimerPairs = true;
		}
		//else throw new EmptyResultSetException("No primer pair appears to be valid according to the pair constraints specified!");
		//this.primerPairs.quickSort();
	}
	
	public int getNumberOfValidPrimerPairs(){
		return this.primerPairs.size();
	}
	
	public PrimerPair getPrimerPair(int i){
		if(!this.hasEnumeratedPrimerPairs) this.enumeratePrimerPairs();
		return (PrimerPair)this.primerPairs.get(i);
	}
	
	public ObjectArrayList getPrimerPairs(){
		if(!this.hasEnumeratedPrimerPairs) this.enumeratePrimerPairs();
		return this.primerPairs;
	}
	
	@Deprecated
	private boolean isValidPrimerPair(PrimerPair pair){
		// melting temperatures of forward and reverse primers should be 'similar'
		// AND TAQMAN probe's TM should be 'sufficiently higher' than the TMs of the forward and reverse primer
		pair.computeMaxAlignmentScore();
		return Math.abs(pair.getForwardPrimer().getMeltingTemp() - pair.getReversePrimer().getMeltingTemp()) <= this.searchParams.getMAX_PRIMER_TM_DIFFERENCE()
			&& pair.getHybridizationProbe().getMeltingTemp() >= pair.getForwardPrimer().getMeltingTemp() + this.searchParams.getMIN_TAQMAN_TM_DIFFERENCE()
			&& pair.getHybridizationProbe().getMeltingTemp() >= pair.getReversePrimer().getMeltingTemp() + this.searchParams.getMIN_TAQMAN_TM_DIFFERENCE()
			&& pair.getMaxPairElementsAlignmentScore().getPairScore() <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE()
			&& pair.getMaxPairElementsAlignmentScore().getPairEndScore() <= this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE();
	}
	
	private boolean isValidPrimerPair(PrimerPair pair, SequenceRegionAlignment fwRev, SequenceRegionAlignment fwProbe, SequenceRegionAlignment revProbe){
		// Test - scan for homogenous forward and reverse primer TMs
		if(Math.abs(pair.getForwardPrimer().getMeltingTemp() - pair.getReversePrimer().getMeltingTemp()) > searchParams.getMAX_PRIMER_TM_DIFFERENCE()){
			this.searchParams.getPickingStat().incPrimerPrimerTMreject();
			return false;
		}
		else{
			this.searchParams.getPickingStat().incPrimerPrimerTMaccept();
		}
		if(searchParams.isPickTaqManProbe() && (Math.abs(pair.getForwardPrimer().getMeltingTemp() - pair.getHybridizationProbe().getMeltingTemp()) < searchParams.getMIN_TAQMAN_TM_DIFFERENCE()
				|| Math.abs(pair.getReversePrimer().getMeltingTemp() - pair.getHybridizationProbe().getMeltingTemp()) < searchParams.getMIN_TAQMAN_TM_DIFFERENCE())){
			this.searchParams.getPickingStat().incPrimerProbeTMreject();
			return false;
		}
		else{
			this.searchParams.getPickingStat().incPrimerProbeTMaccept();
		}
			
		PrimerAlignmentScores fwRevScores = fwRev.getGlobalAlignmentValues(pair.getForwardPrimer().getPositionInScanSequence(), pair.getForwardPrimer().getLength(),pair.getReversePrimer().getPositionInScanSequence(), pair.getReversePrimer().getLength());
		int fwRevPA = fwRevScores.getPairScore();
		int fwRevPEA = fwRevScores.getPairEndScore();
		int fwProbePA = (fwProbe == null)? 0 : fwProbe.getGlobalAlignmentValues(pair.getForwardPrimer().getPositionInScanSequence(), pair.getForwardPrimer().getLength(), pair.getHybridizationProbe().getPositionInScanSequence(), pair.getHybridizationProbe().getLength()).getPairScore();
		int revProbePA = (revProbe == null)? 0 : revProbe.getGlobalAlignmentValues(pair.getReversePrimer().getPositionInScanSequence(), pair.getReversePrimer().getLength(), pair.getHybridizationProbe().getPositionInScanSequence(), pair.getHybridizationProbe().getLength()).getPairScore();
		
//		return Math.abs(pair.getForwardPrimer().getMeltingTemp() - pair.getReversePrimer().getMeltingTemp()) <= this.searchParams.getMAX_PRIMER_TM_DIFFERENCE()
//		&& pair.getHybridizationProbe().getMeltingTemp() >= pair.getForwardPrimer().getMeltingTemp() + this.searchParams.getMIN_TAQMAN_TM_DIFFERENCE()
//		&& pair.getHybridizationProbe().getMeltingTemp() >= pair.getReversePrimer().getMeltingTemp() + this.searchParams.getMIN_TAQMAN_TM_DIFFERENCE()
//		&& fwRevPA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && fwRevPEA <= this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()
//		&& fwProbePA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && revProbePA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE();

		// do stat
		if(fwRevPA > this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || fwRevPEA > this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()) this.searchParams.getPickingStat().incPrimerPrimerAlignmentReject();
		else this.searchParams.getPickingStat().incPrimerPrimerAlignmentAccept();
		
		if(fwProbePA > this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || revProbePA > this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE()) this.searchParams.getPickingStat().incPrimerProbeAlignmentReject();
		else this.searchParams.getPickingStat().incPrimerProbeAlignmentAccept();
		
		// compatibility of melting temperatures has already been checked above!
		return fwRevPA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && fwRevPEA <= this.searchParams.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()
			&& fwProbePA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() && revProbePA <= this.searchParams.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE();
	}
	
	/**
	 * @return the forwardScanSequence
	 */
	public char[] getForwardScanSequence() {
		return forwardScanSequence;
	}

	/**
	 * @param forwardScanSequence the forwardScanSequence to set
	 */
	public void setForwardScanSequence(char[] forwardScanSequence) {
		this.forwardScanSequence = forwardScanSequence;
	}

	/**
	 * @return the reverseScanSequence
	 */
	public char[] getReverseScanSequence() {
		return reverseScanSequence;
	}

	/**
	 * @param reverseScanSequence the reverseScanSequence to set
	 */
	public void setReverseScanSequence(char[] reverseScanSequence) {
		this.reverseScanSequence = reverseScanSequence;
	}

	/**
	 * @return the probeScanSequence
	 */
	public char[] getProbeScanSequence() {
		return probeScanSequence;
	}

	/**
	 * @param probeScanSequence the probeScanSequence to set
	 */
	public void setProbeScanSequence(char[] probeScanSequence) {
		this.probeScanSequence = probeScanSequence;
	}

	/**
	 * Sorts the primer pairs into ascending distance to the optimal primer pair (with distance 0 to itself :-).
	 */
	public void sortPrimerPairs(){
		this.primerPairs.sort();
	}

	public int compareTo(Object o) {
		RestrictionSite otherSite = (RestrictionSite) o;
		if(this.distanceToIntervalMean < otherSite.getDistanceToIntervalMean()) return -1;
		else if(this.distanceToIntervalMean > otherSite.getDistanceToIntervalMean()) return 1;
		else return 0;
	}	
	
	public boolean equals(Object other){
		RestrictionSite otherSite = (RestrictionSite) other;
		if(this.enzyme.equals(otherSite.enzyme) && this.position == otherSite.position
				&& this.distanceToIntervalMean == otherSite.distanceToIntervalMean
				&& this.validUpstreamPrimers.equals(otherSite.validUpstreamPrimers) && this.sequenceRegion.equals(otherSite.sequenceRegion)) return true;
		else return false;
	}
	
}
