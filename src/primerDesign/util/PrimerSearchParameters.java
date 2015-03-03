/**
 * 
 */
package primerDesign.util;

import java.util.Properties;
import java.util.regex.Pattern;

import org.biojava.bio.molbio.RestrictionEnzyme;

import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.algo.PrimerSearch;
import primerDesign.algo.SimpleAlignment;
import primerDesign.dsc.PrimerPairPickingStatistics;
import primerDesign.dsc.PrimerSearchStatistics;
import primerDesign.dsc.indexStructures.TargetOrganisms;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck;

/**
 * Class for storing 3PD parameters.
 * 
 * @author Sebastian Fr�hler
 */
public class PrimerSearchParameters {
	private boolean useIndices = true;  // re-use pre-computed values
	private boolean PRINT_DEBUG_LOG = true; // print some debugging output about rejected primers
	private boolean computeScanningStatistics = true;
	private boolean computePickingStatistics = true;
	
	private PrimerSearchStatistics searchStat;
	private PrimerPairPickingStatistics pickingStat;
	
	// constants for melting p1. calculation - primer concentrations in Mol!!!
	private double PRIMER_CONCENTRATION = 500E-09; // Dekker 2007 protocol (http://www.nature.com/nprot/journal/v2/n7/full/nprot.2007.243.html)
	private double TAQMAN_PROBE_CONCENTRATION = 150E-09; // Dekker 2007 protocol 
	private double MONOVALENT_CATION_CONCENTRATION = 0; //50E-03;
	private double DIVALENT_CATION_CONCENTRATION = 2.5E-03; // Roche probes master: Roche support: ~2,5mM Mg2+
	private double DNTP_CONCENTRATION = 0.4E-03; //?? as in std. PCR ??
	
	// parameters for primer alignment computation
	private int a_t_basepair_score = 2;
	private int g_c_basepair_score = 4;
	private PrimerAlignmentCalculation PRIMER_ALIGNMENT_METHOD = new SimpleAlignment(a_t_basepair_score, g_c_basepair_score); //new KaempkePrimerAlignment();
	private boolean checkInterPairAlignments = true;
	private boolean checkIntraPairAlignments = true;
	
	private PrimerMisprimingCheck primerMisprimingCheck;
	private Properties configFile;
	
	// parameters used for primer design
	// constants used for primer design ('optimal primer properties')
	// 2Do: setup parameters like in: http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html ??
	
	// primer3 primer parameters are for breslauer TM calc method!!! -> Santa Lucia TM ~ Breslauer TM + 3-5�C 
	// source: http://frodo.wi.mit.edu/
	
	private int MAX_PRIMER_LENGTH = 30; // primer3: 27
	private int MIN_PRIMER_LENGTH = 16; // primer3: 18
	private int OPT_PRIMER_LENGTH = 20; // primer3: 20
	
	private int TAQMAN_MAX_PRIMER_LENGTH = 35; // 30 // primer3: 27
	private int TAQMAN_MIN_PRIMER_LENGTH = 18; // primer3: 18
	private int TAQMAN_OPT_PRIMER_LENGTH = 25; // ?? // primer3: 20
	
	private int MAX_AMPLICON_LENGTH = 200; //200 //250;  //130 // 120?
	private int MIN_AMPLICON_LENGTH = 140; //140;   // 70 // 80?
	private int OPT_AMPLICON_LENGTH = 150; //150 //200;  // 100?
	
	// primer mispriming parameters
	private int SAFE_FALSE_POSITIVE_AMPLICON_LENGTH = 500; //1000;
	private int MIN_MISPRIMING_TM_DIFFERENCE = 10;
	// there are two possible mispriming scans:
	// - exact mismatch of primer end of given length
	private int PRIMER_END_MISMATCH_SCAN_LENGTH = 11;  // as in blast for DNA
	// - approximate mismatch search analogous to blast (with DP verify - slow!)
	private int MAX_PRIMER_MISPRIMING_CUTOFF = 0;
	private int SEED_WORD_SIZE = 11;
	private int SEED_STEP_SIZE = 1;
	private int THREE_PRIME_DP_LENGTH = 5;
	private int THREE_PRIME_DP_THRESHOLD = 1;
	private int WHOLE_DP_THRESHOLD = (MAX_PRIMER_LENGTH + MIN_PRIMER_LENGTH) / 4;
	
	private int OPT_DISTANCE_TO_RSS = OPT_AMPLICON_LENGTH/2;
	private int MIN_RESTRICTION_FRAGMENT_LENGTH = SAFE_FALSE_POSITIVE_AMPLICON_LENGTH; // default: 1000
	private int MAX_RESTRICTION_FRAGMENT_LENGTH = Integer.MAX_VALUE;
	private int MAX_POLY_X_LENGTH = 5; // default: 5, maximum length+1 of single nocleotide repeats // primer3: 5
	private Pattern SINGLE_NUCLEOTIDE_REPEAT_PATTERN = Pattern.compile(".*([A]{" + MAX_POLY_X_LENGTH +  ",}|[T]{" + MAX_POLY_X_LENGTH +  ",}|[G]{" + MAX_POLY_X_LENGTH +  ",}|[C]{" + MAX_POLY_X_LENGTH +  ",}).*");
	private int GC_CLAMP_LENGTH = 5; // default: 5
	private int MIN_GC_CLAMP = 0; // default: 1, the minimum number of G|C within the last GC_CLAMP_LENGTH basepairs of a primer (3' end)
	private int MAX_GC_CLAMP = 5; // default: 3
	
	private double MAX_TM = 63; // 63 ?? // primer3: 63
	private double MIN_TM = 55; // 55 ?? // primer3: 57
	private double OPT_TM = 60; // 59 ?? // primer3: 60
	
	private double TAQMAN_MAX_TM = 75; // 70 // primer3: 63
	private double TAQMAN_MIN_TM = 65; //65; //63?? 68 // primer3: 57
	private double TAQMAN_OPT_TM = 70; // 69 // primer3: 60
	
	private double MAX_PRIMER_TM_DIFFERENCE = 1;  // primer pair (FW and Rev) TMs should not differ by more than x �C
	private double MIN_TAQMAN_TM_DIFFERENCE = 5;  // TAQMAN TM in a primer set has to be at least x �C higher than TMs of primer pairs
	
	private double MAX_GC = 0.7; // ?? // primer3: 80
	private double MIN_GC = 0.3; // ?? // primer3: 20
	private double OPT_GC = 0.5; // ?? // primer3: NOT specified
	
	private double TAQMAN_MAX_GC = 0.7; // 0.80 // primer3: 80
	private double TAQMAN_MIN_GC = 0.3; // 0.30 // primer3: 20
	private double TAQMAN_OPT_GC = 0.5; // 0.55 // primer3: NOT specified
	
	private int MAX_PRIMER_SELF_ALIGNMENT_SCORE = 30; //30; // 24 //??
	private int OPT_PRIMER_SELF_ALIGNMENT_SCORE = 0;
	private int MAX_PRIMER_SELF_END_ALIGNMENT_SCORE = 10; //10; // 8//??
	private int OPT_PRIMER_SELF_END_ALIGNMENT_SCORE = 0;
	
	private int MAX_PRIMER_PAIR_ALIGNMENT_SCORE = 30; //30; // 24 //??
	private int OPT_PRIMER_PAIR_ALIGNMENT_SCORE = 0;
	private int MAX_PRIMER_PAIR_END_ALIGNMENT_SCORE = 10; //10; // 8 //??
	private int OPT_PRIMER_PAIR_END_ALIGNMENT_SCORE = 0;
	
	private int MAX_TAQMAN_SELF_ALIGNMENT_SCORE = 36; //??
	private int OPT_TAQMAN_SELF_ALIGNMENT_SCORE = 0;
	private int MAX_TAQMAN_SELF_END_ALIGNMENT_SCORE = 36; //24??
	private int OPT_TAQMAN_SELF_END_ALIGNMENT_SCORE = 0;
	
	private int MAX_TAQMAN_PAIR_ALIGNMENT_SCORE = 36; //??
	private int OPT_TAQMAN_PAIR_ALIGNMENT_SCORE = 0;
	private int MAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE = 36; //24??
	private int OPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE = 0;
	
	// weights for score & distance computation ('optimal primer properties weights') - equal scores lead to equal weighting of the respective terms
	private double PRIMER_DELTA_TM_WEIGHT = 1;  // 1
	private double PRIMER_DELTA_GC_WEIGHT = 0.5; // 0.5
	private double PRIMER_DELTA_LENGTH_WEIGHT = 0.5; // 0.5
	private double PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT = 1; // 1
	private double SELF_ALIGNMENT_WEIGHT = 0.5; // 0.1
	private double SELF_END_ALIGNMENT_WEIGHT = 1; // 0.2 // w(SEA) sould be 2 w(SA)! (K�mpke et.al.)
	private double PAIR_ALIGNMENT_WEIGHT = 0.5; // 0.1
	private double PAIR_END_ALIGNMENT_WEIGHT = 1; // 0.2 // w(PEA) sould be 2 w(PA)! (K�mpke et.al.) 
	private double PRIMER_FALSE_POSITIVES_WEIGHT = 1; // 1
	
	private double PRIMER_PAIR_HOMOGENITY_WEIGHT = 1;
	private double PRIMER_PAIR_DOPT_WEIGHT = 1;
	
	private Enum<TargetOrganisms> targetOrganism;
	private SimpleContig[] contigs;
	private int numPrimers;
	private RestrictionEnzyme enzyme;
	
	private PrimerSearch primerSearch;
	
	private boolean pickTaqManProbe = false;
	
	public PrimerSearchParameters(){
		this.searchStat = new PrimerSearchStatistics();
		this.pickingStat = new PrimerPairPickingStatistics();
	}
	
	public PrimerSearchParameters(PrimerSearch search){
		this.primerSearch = search;
		this.searchStat = new PrimerSearchStatistics();
		this.pickingStat = new PrimerPairPickingStatistics();
	}
	
	public PrimerSearch getPrimerSearch() {
		return primerSearch;
	}
	
	public void setPrimerSearch(PrimerSearch search){
		this.primerSearch = search;
	}
	
	public int getNumPrimers() {
		return numPrimers;
	}
	public void setNumPrimers(int numPrimers) {
		if(numPrimers < 1) throw new IllegalArgumentException("The number of primers must be >= 1!");
		this.numPrimers = numPrimers;
	}
	public RestrictionEnzyme getEnzyme() {
		return enzyme;
	}
	public void setEnzyme(RestrictionEnzyme enzyme) {
		if(enzyme == null) throw new NullPointerException();
		this.enzyme = enzyme;
	}
	
//	public String getSequence() {
//		return sequence;
//	}
//	4
//	public void setSequence(String sequence) {
//		this.sequence = sequence;
//	}
	
	public void setContigs(SimpleContig[] contigs){
		this.contigs = contigs;
	}
	public SimpleContig[] getContigs(){
		return this.contigs;
	}
	public int getNumContigs(){
		return this.contigs.length;
	}
	
	/**
	 * @return the useIndices
	 */
	public boolean isUseIndices() {
		return useIndices;
	}
	/**
	 * @param useIndices the useIndices to set
	 */
	public void setUseIndices(boolean useIndices) {
		this.useIndices = useIndices;
	}
	/**
	 * @return the pRINT_DEBUG_LOG
	 */
	public boolean isPRINT_DEBUG_LOG() {
		return PRINT_DEBUG_LOG;
	}
	/**
	 * @param print_debug_log the pRINT_DEBUG_LOG to set
	 */
	public void setPRINT_DEBUG_LOG(boolean print_debug_log) {
		PRINT_DEBUG_LOG = print_debug_log;
	}
	/**
	 * @return the pRIMER_CONCENTRATION
	 */
	public double getPRIMER_CONCENTRATION() {
		return PRIMER_CONCENTRATION;
	}
	/**
	 * @param primer_concentration the pRIMER_CONCENTRATION to set
	 */
	public void setPRIMER_CONCENTRATION(double primer_concentration) {
		PRIMER_CONCENTRATION = primer_concentration;
	}
	/**
	 * @return the tAQMAN_PROBE_CONCENTRATION
	 */
	public double getTAQMAN_PROBE_CONCENTRATION() {
		return TAQMAN_PROBE_CONCENTRATION;
	}
	/**
	 * @param taqman_probe_concentration the tAQMAN_PROBE_CONCENTRATION to set
	 */
	public void setTAQMAN_PROBE_CONCENTRATION(double taqman_probe_concentration) {
		TAQMAN_PROBE_CONCENTRATION = taqman_probe_concentration;
	}
	/**
	 * @return the mONOVALENT_CATION_CONCENTRATION
	 */
	public double getMONOVALENT_CATION_CONCENTRATION() {
		return MONOVALENT_CATION_CONCENTRATION;
	}
	/**
	 * @param monovalent_cation_concentration the mONOVALENT_CATION_CONCENTRATION to set
	 */
	public void setMONOVALENT_CATION_CONCENTRATION(
			double monovalent_cation_concentration) {
		MONOVALENT_CATION_CONCENTRATION = monovalent_cation_concentration;
	}
	/**
	 * @return the dIVALENT_CATION_CONCENTRATION
	 */
	public double getDIVALENT_CATION_CONCENTRATION() {
		return DIVALENT_CATION_CONCENTRATION;
	}
	/**
	 * @param divalent_cation_concentration the dIVALENT_CATION_CONCENTRATION to set
	 */
	public void setDIVALENT_CATION_CONCENTRATION(
			double divalent_cation_concentration) {
		DIVALENT_CATION_CONCENTRATION = divalent_cation_concentration;
	}
	/**
	 * @return the dNTP_CONCENTRATION
	 */
	public double getDNTP_CONCENTRATION() {
		return DNTP_CONCENTRATION;
	}
	/**
	 * @param dntp_concentration the dNTP_CONCENTRATION to set
	 */
	public void setDNTP_CONCENTRATION(double dntp_concentration) {
		DNTP_CONCENTRATION = dntp_concentration;
	}
	/**
	 * @return the pRIMER_ALIGNMENT_METHOD
	 */
	public PrimerAlignmentCalculation getPRIMER_ALIGNMENT_METHOD() {
		return PRIMER_ALIGNMENT_METHOD;
	}
	/**
	 * @param primer_alignment_method the pRIMER_ALIGNMENT_METHOD to set
	 */
	public void setPRIMER_ALIGNMENT_METHOD(
			PrimerAlignmentCalculation primer_alignment_method) {
		PRIMER_ALIGNMENT_METHOD = primer_alignment_method;
	}
	/**
	 * @return the a_t_basepair_score
	 */
	public int getA_t_basepair_score() {
		return a_t_basepair_score;
	}
	/**
	 * @param a_t_basepair_score the a_t_basepair_score to set
	 */
	public void setA_t_basepair_score(int a_t_basepair_score) {
		this.a_t_basepair_score = a_t_basepair_score;
	}
	/**
	 * @return the g_c_basepair_score
	 */
	public int getG_c_basepair_score() {
		return g_c_basepair_score;
	}
	/**
	 * @param g_c_basepair_score the g_c_basepair_score to set
	 */
	public void setG_c_basepair_score(int g_c_basepair_score) {
		this.g_c_basepair_score = g_c_basepair_score;
	}
	/**
	 * @return the mAX_PRIMER_LENGTH
	 */
	public int getMAX_PRIMER_LENGTH() {
		return MAX_PRIMER_LENGTH;
	}
	/**
	 * @param max_primer_length the mAX_PRIMER_LENGTH to set
	 */
	public void setMAX_PRIMER_LENGTH(int max_primer_length) {
		MAX_PRIMER_LENGTH = max_primer_length;
	}
	/**
	 * @return the mIN_PRIMER_LENGTH
	 */
	public int getMIN_PRIMER_LENGTH() {
		return MIN_PRIMER_LENGTH;
	}
	/**
	 * @param min_primer_length the mIN_PRIMER_LENGTH to set
	 */
	public void setMIN_PRIMER_LENGTH(int min_primer_length) {
		MIN_PRIMER_LENGTH = min_primer_length;
	}
	/**
	 * @return the oPT_PRIMER_LENGTH
	 */
	public int getOPT_PRIMER_LENGTH() {
		return OPT_PRIMER_LENGTH;
	}
	/**
	 * @param opt_primer_length the oPT_PRIMER_LENGTH to set
	 */
	public void setOPT_PRIMER_LENGTH(int opt_primer_length) {
		OPT_PRIMER_LENGTH = opt_primer_length;
	}
	/**
	 * @return the tAQMAN_MAX_PRIMER_LENGTH
	 */
	public int getTAQMAN_MAX_PRIMER_LENGTH() {
		return TAQMAN_MAX_PRIMER_LENGTH;
	}
	/**
	 * @param taqman_max_primer_length the tAQMAN_MAX_PRIMER_LENGTH to set
	 */
	public void setTAQMAN_MAX_PRIMER_LENGTH(int taqman_max_primer_length) {
		TAQMAN_MAX_PRIMER_LENGTH = taqman_max_primer_length;
	}
	/**
	 * @return the tAQMAN_MIN_PRIMER_LENGTH
	 */
	public int getTAQMAN_MIN_PRIMER_LENGTH() {
		return TAQMAN_MIN_PRIMER_LENGTH;
	}
	/**
	 * @param taqman_min_primer_length the tAQMAN_MIN_PRIMER_LENGTH to set
	 */
	public void setTAQMAN_MIN_PRIMER_LENGTH(int taqman_min_primer_length) {
		TAQMAN_MIN_PRIMER_LENGTH = taqman_min_primer_length;
	}
	/**
	 * @return the tAQMAN_OPT_PRIMER_LENGTH
	 */
	public int getTAQMAN_OPT_PRIMER_LENGTH() {
		return TAQMAN_OPT_PRIMER_LENGTH;
	}
	/**
	 * @param taqman_opt_primer_length the tAQMAN_OPT_PRIMER_LENGTH to set
	 */
	public void setTAQMAN_OPT_PRIMER_LENGTH(int taqman_opt_primer_length) {
		TAQMAN_OPT_PRIMER_LENGTH = taqman_opt_primer_length;
	}
	/**
	 * @return the mAX_AMPLICON_LENGTH
	 */
	public int getMAX_AMPLICON_LENGTH() {
		return MAX_AMPLICON_LENGTH;
	}
	/**
	 * @param max_amplicon_length the mAX_AMPLICON_LENGTH to set
	 */
	public void setMAX_AMPLICON_LENGTH(int max_amplicon_length) {
		MAX_AMPLICON_LENGTH = max_amplicon_length;
	}
	/**
	 * @return the mIN_AMPLICON_LENGTH
	 */
	public int getMIN_AMPLICON_LENGTH() {
		return MIN_AMPLICON_LENGTH;
	}
	/**
	 * @param min_amplicon_length the mIN_AMPLICON_LENGTH to set
	 */
	public void setMIN_AMPLICON_LENGTH(int min_amplicon_length) {
		MIN_AMPLICON_LENGTH = min_amplicon_length;
	}
	/**
	 * @return the oPT_AMPLICON_LENGTH
	 */
	public int getOPT_AMPLICON_LENGTH() {
		return OPT_AMPLICON_LENGTH;
	}
	/**
	 * @param opt_amplicon_length the oPT_AMPLICON_LENGTH to set
	 */
	public void setOPT_AMPLICON_LENGTH(int opt_amplicon_length) {
		OPT_AMPLICON_LENGTH = opt_amplicon_length;
	}
	/**
	 * @return the sAFE_FALSE_POSITIVE_AMPLICON_LENGTH
	 */
	public int getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH() {
		return SAFE_FALSE_POSITIVE_AMPLICON_LENGTH;
	}
	/**
	 * @param safe_false_positive_amplicon_length the sAFE_FALSE_POSITIVE_AMPLICON_LENGTH to set
	 */
	public void setSAFE_FALSE_POSITIVE_AMPLICON_LENGTH(
			int safe_false_positive_amplicon_length) {
		SAFE_FALSE_POSITIVE_AMPLICON_LENGTH = safe_false_positive_amplicon_length;
	}
	/**
	 * @return the pRIMER_END_MISMATCH_SCAN_LENGTH
	 */
	public int getPRIMER_END_MISMATCH_SCAN_LENGTH() {
		return PRIMER_END_MISMATCH_SCAN_LENGTH;
	}
	/**
	 * @param primer_end_mismatch_scan_length the pRIMER_END_MISMATCH_SCAN_LENGTH to set
	 */
	public void setPRIMER_END_MISMATCH_SCAN_LENGTH(
			int primer_end_mismatch_scan_length) {
		PRIMER_END_MISMATCH_SCAN_LENGTH = primer_end_mismatch_scan_length;
	}
	/**
	 * @return the mAX_PRIMER_MISPRIMING_CUTOFF
	 */
	public int getMAX_PRIMER_MISPRIMING_CUTOFF() {
		return MAX_PRIMER_MISPRIMING_CUTOFF;
	}
	/**
	 * @param max_primer_mispriming_cutoff the mAX_PRIMER_MISPRIMING_CUTOFF to set
	 */
	public void setMAX_PRIMER_MISPRIMING_CUTOFF(int max_primer_mispriming_cutoff) {
		MAX_PRIMER_MISPRIMING_CUTOFF = max_primer_mispriming_cutoff;
	}
	/**
	 * @return the sEED_WORD_SIZE
	 */
	public int getSEED_WORD_SIZE() {
		return SEED_WORD_SIZE;
	}
	/**
	 * @param seed_word_size the sEED_WORD_SIZE to set
	 */
	public void setSEED_WORD_SIZE(int seed_word_size) {
		SEED_WORD_SIZE = seed_word_size;
	}
	/**
	 * @return the sEED_STEP_SIZE
	 */
	public int getSEED_STEP_SIZE() {
		return SEED_STEP_SIZE;
	}
	/**
	 * @param seed_step_size the sEED_STEP_SIZE to set
	 */
	public void setSEED_STEP_SIZE(int seed_step_size) {
		SEED_STEP_SIZE = seed_step_size;
	}
	/**
	 * @return the tHREE_PRIME_DP_LENGTH
	 */
	public int getTHREE_PRIME_DP_LENGTH() {
		return THREE_PRIME_DP_LENGTH;
	}
	/**
	 * @param three_prime_dp_length the tHREE_PRIME_DP_LENGTH to set
	 */
	public void setTHREE_PRIME_DP_LENGTH(int three_prime_dp_length) {
		THREE_PRIME_DP_LENGTH = three_prime_dp_length;
	}
	/**
	 * @return the tHREE_PRIME_DP_THRESHOLD
	 */
	public int getTHREE_PRIME_DP_THRESHOLD() {
		return THREE_PRIME_DP_THRESHOLD;
	}
	/**
	 * @param three_prime_dp_threshold the tHREE_PRIME_DP_THRESHOLD to set
	 */
	public void setTHREE_PRIME_DP_THRESHOLD(int three_prime_dp_threshold) {
		THREE_PRIME_DP_THRESHOLD = three_prime_dp_threshold;
	}
	/**
	 * @return the wHOLE_DP_THRESHOLD
	 */
	public int getWHOLE_DP_THRESHOLD() {
		return WHOLE_DP_THRESHOLD;
	}
	/**
	 * @param whole_dp_threshold the wHOLE_DP_THRESHOLD to set
	 */
	public void setWHOLE_DP_THRESHOLD(int whole_dp_threshold) {
		WHOLE_DP_THRESHOLD = whole_dp_threshold;
	}
	/**
	 * @return the oPT_DISTANCE_TO_RSS
	 */
	public int getOPT_DISTANCE_TO_RSS() {
		return OPT_DISTANCE_TO_RSS;
	}
	/**
	 * @param opt_distance_to_rss the oPT_DISTANCE_TO_RSS to set
	 */
	public void setOPT_DISTANCE_TO_RSS(int opt_distance_to_rss) {
		OPT_DISTANCE_TO_RSS = opt_distance_to_rss;
	}
	/**
	 * @return the mIN_RESTRICTION_FRAGMENT_LENGTH
	 */
	public int getMIN_RESTRICTION_FRAGMENT_LENGTH() {
		return MIN_RESTRICTION_FRAGMENT_LENGTH;
	}
	/**
	 * @param min_restriction_fragment_length the mIN_RESTRICTION_FRAGMENT_LENGTH to set
	 */
	public void setMIN_RESTRICTION_FRAGMENT_LENGTH(int min_restriction_fragment_length) {
		MIN_RESTRICTION_FRAGMENT_LENGTH = min_restriction_fragment_length;
	}
	/**
	 * @return returns the the maximum size of acceptable restriction fragments
	 */
	public int getMAX_RESTRICTION_FRAGMENT_LENGTH(){
		return this.MAX_RESTRICTION_FRAGMENT_LENGTH;
	}
	/**
	 * 
	 * @param max_restriction_fragment_length the maximum size of acceptable restriction fragments
	 */
	public void setMAX_RESTRICTION_FRAGMENT_LENGTH(int max_restriction_fragment_length){
		this.MAX_RESTRICTION_FRAGMENT_LENGTH = max_restriction_fragment_length;
	}
	
	/**
	 * @return the mAX_POLY_X_LENGTH
	 */
	public int getMAX_POLY_X_LENGTH() {
		return MAX_POLY_X_LENGTH;
	}
	/**
	 * @param max_poly_x_length the mAX_POLY_X_LENGTH to set
	 */
	public void setMAX_POLY_X_LENGTH(int max_poly_x_length) {
		MAX_POLY_X_LENGTH = max_poly_x_length;
	}
	/**
	 * @return the gC_CLAMP_LENGTH
	 */
	public int getGC_CLAMP_LENGTH() {
		return GC_CLAMP_LENGTH;
	}
	/**
	 * @param gc_clamp_length the gC_CLAMP_LENGTH to set
	 */
	public void setGC_CLAMP_LENGTH(int gc_clamp_length) {
		GC_CLAMP_LENGTH = gc_clamp_length;
	}
	/**
	 * @return the mIN_GC_CLAMP
	 */
	public int getMIN_GC_CLAMP() {
		return MIN_GC_CLAMP;
	}
	/**
	 * @param min_gc_clamp the mIN_GC_CLAMP to set
	 */
	public void setMIN_GC_CLAMP(int min_gc_clamp) {
		MIN_GC_CLAMP = min_gc_clamp;
	}
	/**
	 * @return the mAX_GC_CLAMP
	 */
	public int getMAX_GC_CLAMP() {
		return MAX_GC_CLAMP;
	}
	/**
	 * @param max_gc_clamp the mAX_GC_CLAMP to set
	 */
	public void setMAX_GC_CLAMP(int max_gc_clamp) {
		MAX_GC_CLAMP = max_gc_clamp;
	}
	/**
	 * @return the mAX_TM
	 */
	public double getMAX_TM() {
		return MAX_TM;
	}
	/**
	 * @param max_tm the mAX_TM to set
	 */
	public void setMAX_TM(double max_tm) {
		MAX_TM = max_tm;
	}
	/**
	 * @return the mIN_TM
	 */
	public double getMIN_TM() {
		return MIN_TM;
	}
	/**
	 * @param min_tm the mIN_TM to set
	 */
	public void setMIN_TM(double min_tm) {
		MIN_TM = min_tm;
	}
	/**
	 * @return the oPT_TM
	 */
	public double getOPT_TM() {
		return OPT_TM;
	}
	/**
	 * @param opt_tm the oPT_TM to set
	 */
	public void setOPT_TM(double opt_tm) {
		OPT_TM = opt_tm;
	}
	/**
	 * @return the tAQMAN_MAX_TM
	 */
	public double getTAQMAN_MAX_TM() {
		return TAQMAN_MAX_TM;
	}
	/**
	 * @param taqman_max_tm the tAQMAN_MAX_TM to set
	 */
	public void setTAQMAN_MAX_TM(double taqman_max_tm) {
		TAQMAN_MAX_TM = taqman_max_tm;
	}
	/**
	 * @return the tAQMAN_MIN_TM
	 */
	public double getTAQMAN_MIN_TM() {
		return TAQMAN_MIN_TM;
	}
	/**
	 * @param taqman_min_tm the tAQMAN_MIN_TM to set
	 */
	public void setTAQMAN_MIN_TM(double taqman_min_tm) {
		TAQMAN_MIN_TM = taqman_min_tm;
	}
	/**
	 * @return the tAQMAN_OPT_TM
	 */
	public double getTAQMAN_OPT_TM() {
		return TAQMAN_OPT_TM;
	}
	/**
	 * @param taqman_opt_tm the tAQMAN_OPT_TM to set
	 */
	public void setTAQMAN_OPT_TM(double taqman_opt_tm) {
		TAQMAN_OPT_TM = taqman_opt_tm;
	}
	/**
	 * @return the mAX_PRIMER_TM_DIFFERENCE
	 */
	public double getMAX_PRIMER_TM_DIFFERENCE() {
		return MAX_PRIMER_TM_DIFFERENCE;
	}
	/**
	 * @param max_primer_tm_difference the mAX_PRIMER_TM_DIFFERENCE to set
	 */
	public void setMAX_PRIMER_TM_DIFFERENCE(double max_primer_tm_difference) {
		MAX_PRIMER_TM_DIFFERENCE = max_primer_tm_difference;
	}
	/**
	 * @return the mIN_TAQMAN_TM_DIFFERENCE
	 */
	public double getMIN_TAQMAN_TM_DIFFERENCE() {
		return MIN_TAQMAN_TM_DIFFERENCE;
	}
	/**
	 * @param min_taqman_tm_difference the mIN_TAQMAN_TM_DIFFERENCE to set
	 */
	public void setMIN_TAQMAN_TM_DIFFERENCE(double min_taqman_tm_difference) {
		MIN_TAQMAN_TM_DIFFERENCE = min_taqman_tm_difference;
	}
	/**
	 * @return the mAX_GC
	 */
	public double getMAX_GC() {
		return MAX_GC;
	}
	/**
	 * @param max_gc the mAX_GC to set
	 */
	public void setMAX_GC(double max_gc) {
		MAX_GC = max_gc;
	}
	/**
	 * @return the mIN_GC
	 */
	public double getMIN_GC() {
		return MIN_GC;
	}
	/**
	 * @param min_gc the mIN_GC to set
	 */
	public void setMIN_GC(double min_gc) {
		MIN_GC = min_gc;
	}
	/**
	 * @return the oPT_GC
	 */
	public double getOPT_GC() {
		return OPT_GC;
	}
	/**
	 * @param opt_gc the oPT_GC to set
	 */
	public void setOPT_GC(double opt_gc) {
		OPT_GC = opt_gc;
	}
	/**
	 * @return the tAQMAN_MAX_GC
	 */
	public double getTAQMAN_MAX_GC() {
		return TAQMAN_MAX_GC;
	}
	/**
	 * @param taqman_max_gc the tAQMAN_MAX_GC to set
	 */
	public void setTAQMAN_MAX_GC(double taqman_max_gc) {
		TAQMAN_MAX_GC = taqman_max_gc;
	}
	/**
	 * @return the tAQMAN_MIN_GC
	 */
	public double getTAQMAN_MIN_GC() {
		return TAQMAN_MIN_GC;
	}
	/**
	 * @param taqman_min_gc the tAQMAN_MIN_GC to set
	 */
	public void setTAQMAN_MIN_GC(double taqman_min_gc) {
		TAQMAN_MIN_GC = taqman_min_gc;
	}
	/**
	 * @return the tAQMAN_OPT_GC
	 */
	public double getTAQMAN_OPT_GC() {
		return TAQMAN_OPT_GC;
	}
	/**
	 * @param taqman_opt_gc the tAQMAN_OPT_GC to set
	 */
	public void setTAQMAN_OPT_GC(double taqman_opt_gc) {
		TAQMAN_OPT_GC = taqman_opt_gc;
	}
	/**
	 * @return the mAX_PRIMER_SELF_ALIGNMENT_SCORE
	 */
	public int getMAX_PRIMER_SELF_ALIGNMENT_SCORE() {
		return MAX_PRIMER_SELF_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_primer_self_alignment_score the mAX_PRIMER_SELF_ALIGNMENT_SCORE to set
	 */
	public void setMAX_PRIMER_SELF_ALIGNMENT_SCORE(
			int max_primer_self_alignment_score) {
		MAX_PRIMER_SELF_ALIGNMENT_SCORE = max_primer_self_alignment_score;
	}
	/**
	 * @return the oPT_PRIMER_SELF_ALIGNMENT_SCORE
	 */
	public int getOPT_PRIMER_SELF_ALIGNMENT_SCORE() {
		return OPT_PRIMER_SELF_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_primer_self_alignment_score the oPT_PRIMER_SELF_ALIGNMENT_SCORE to set
	 */
	public void setOPT_PRIMER_SELF_ALIGNMENT_SCORE(
			int opt_primer_self_alignment_score) {
		OPT_PRIMER_SELF_ALIGNMENT_SCORE = opt_primer_self_alignment_score;
	}
	/**
	 * @return the mAX_PRIMER_SELF_END_ALIGNMENT_SCORE
	 */
	public int getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE() {
		return MAX_PRIMER_SELF_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_primer_self_end_alignment_score the mAX_PRIMER_SELF_END_ALIGNMENT_SCORE to set
	 */
	public void setMAX_PRIMER_SELF_END_ALIGNMENT_SCORE(
			int max_primer_self_end_alignment_score) {
		MAX_PRIMER_SELF_END_ALIGNMENT_SCORE = max_primer_self_end_alignment_score;
	}
	/**
	 * @return the oPT_PRIMER_SELF_END_ALIGNMENT_SCORE
	 */
	public int getOPT_PRIMER_SELF_END_ALIGNMENT_SCORE() {
		return OPT_PRIMER_SELF_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_primer_self_end_alignment_score the oPT_PRIMER_SELF_END_ALIGNMENT_SCORE to set
	 */
	public void setOPT_PRIMER_SELF_END_ALIGNMENT_SCORE(
			int opt_primer_self_end_alignment_score) {
		OPT_PRIMER_SELF_END_ALIGNMENT_SCORE = opt_primer_self_end_alignment_score;
	}
	/**
	 * @return the mAX_PRIMER_PAIR_ALIGNMENT_SCORE
	 */
	public int getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() {
		return MAX_PRIMER_PAIR_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_primer_pair_alignment_score the mAX_PRIMER_PAIR_ALIGNMENT_SCORE to set
	 */
	public void setMAX_PRIMER_PAIR_ALIGNMENT_SCORE(
			int max_primer_pair_alignment_score) {
		MAX_PRIMER_PAIR_ALIGNMENT_SCORE = max_primer_pair_alignment_score;
	}
	/**
	 * @return the oPT_PRIMER_PAIR_ALIGNMENT_SCORE
	 */
	public int getOPT_PRIMER_PAIR_ALIGNMENT_SCORE() {
		return OPT_PRIMER_PAIR_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_primer_pair_alignment_score the oPT_PRIMER_PAIR_ALIGNMENT_SCORE to set
	 */
	public void setOPT_PRIMER_PAIR_ALIGNMENT_SCORE(
			int opt_primer_pair_alignment_score) {
		OPT_PRIMER_PAIR_ALIGNMENT_SCORE = opt_primer_pair_alignment_score;
	}
	/**
	 * @return the mAX_PRIMER_PAIR_END_ALIGNMENT_SCORE
	 */
	public int getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() {
		return MAX_PRIMER_PAIR_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_primer_pair_end_alignment_score the mAX_PRIMER_PAIR_END_ALIGNMENT_SCORE to set
	 */
	public void setMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE(
			int max_primer_pair_end_alignment_score) {
		MAX_PRIMER_PAIR_END_ALIGNMENT_SCORE = max_primer_pair_end_alignment_score;
	}
	/**
	 * @return the oPT_PRIMER_PAIR_END_ALIGNMENT_SCORE
	 */
	public int getOPT_PRIMER_PAIR_END_ALIGNMENT_SCORE() {
		return OPT_PRIMER_PAIR_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_primer_pair_end_alignment_score the oPT_PRIMER_PAIR_END_ALIGNMENT_SCORE to set
	 */
	public void setOPT_PRIMER_PAIR_END_ALIGNMENT_SCORE(
			int opt_primer_pair_end_alignment_score) {
		OPT_PRIMER_PAIR_END_ALIGNMENT_SCORE = opt_primer_pair_end_alignment_score;
	}
	/**
	 * @return the mAX_TAQMAN_SELF_ALIGNMENT_SCORE
	 */
	public int getMAX_TAQMAN_SELF_ALIGNMENT_SCORE() {
		return MAX_TAQMAN_SELF_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_taqman_self_alignment_score the mAX_TAQMAN_SELF_ALIGNMENT_SCORE to set
	 */
	public void setMAX_TAQMAN_SELF_ALIGNMENT_SCORE(
			int max_taqman_self_alignment_score) {
		MAX_TAQMAN_SELF_ALIGNMENT_SCORE = max_taqman_self_alignment_score;
	}
	/**
	 * @return the oPT_TAQMAN_SELF_ALIGNMENT_SCORE
	 */
	public int getOPT_TAQMAN_SELF_ALIGNMENT_SCORE() {
		return OPT_TAQMAN_SELF_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_taqman_self_alignment_score the oPT_TAQMAN_SELF_ALIGNMENT_SCORE to set
	 */
	public void setOPT_TAQMAN_SELF_ALIGNMENT_SCORE(
			int opt_taqman_self_alignment_score) {
		OPT_TAQMAN_SELF_ALIGNMENT_SCORE = opt_taqman_self_alignment_score;
	}
	/**
	 * @return the mAX_TAQMAN_SELF_END_ALIGNMENT_SCORE
	 */
	public int getMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE() {
		return MAX_TAQMAN_SELF_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_taqman_self_end_alignment_score the mAX_TAQMAN_SELF_END_ALIGNMENT_SCORE to set
	 */
	public void setMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE(
			int max_taqman_self_end_alignment_score) {
		MAX_TAQMAN_SELF_END_ALIGNMENT_SCORE = max_taqman_self_end_alignment_score;
	}
	/**
	 * @return the oPT_TAQMAN_SELF_END_ALIGNMENT_SCORE
	 */
	public int getOPT_TAQMAN_SELF_END_ALIGNMENT_SCORE() {
		return OPT_TAQMAN_SELF_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_taqman_self_end_alignment_score the oPT_TAQMAN_SELF_END_ALIGNMENT_SCORE to set
	 */
	public void setOPT_TAQMAN_SELF_END_ALIGNMENT_SCORE(
			int opt_taqman_self_end_alignment_score) {
		OPT_TAQMAN_SELF_END_ALIGNMENT_SCORE = opt_taqman_self_end_alignment_score;
	}
	/**
	 * @return the mAX_TAQMAN_PAIR_ALIGNMENT_SCORE
	 */
	public int getMAX_TAQMAN_PAIR_ALIGNMENT_SCORE() {
		return MAX_TAQMAN_PAIR_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_taqman_pair_alignment_score the mAX_TAQMAN_PAIR_ALIGNMENT_SCORE to set
	 */
	public void setMAX_TAQMAN_PAIR_ALIGNMENT_SCORE(
			int max_taqman_pair_alignment_score) {
		MAX_TAQMAN_PAIR_ALIGNMENT_SCORE = max_taqman_pair_alignment_score;
	}
	/**
	 * @return the oPT_TAQMAN_PAIR_ALIGNMENT_SCORE
	 */
	public int getOPT_TAQMAN_PAIR_ALIGNMENT_SCORE() {
		return OPT_TAQMAN_PAIR_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_taqman_pair_alignment_score the oPT_TAQMAN_PAIR_ALIGNMENT_SCORE to set
	 */
	public void setOPT_TAQMAN_PAIR_ALIGNMENT_SCORE(
			int opt_taqman_pair_alignment_score) {
		OPT_TAQMAN_PAIR_ALIGNMENT_SCORE = opt_taqman_pair_alignment_score;
	}
	/**
	 * @return the mAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE
	 */
	public int getMAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE() {
		return MAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param max_taqman_pair_end_alignment_score the mAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE to set
	 */
	public void setMAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE(
			int max_taqman_pair_end_alignment_score) {
		MAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE = max_taqman_pair_end_alignment_score;
	}
	/**
	 * @return the oPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE
	 */
	public int getOPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE() {
		return OPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE;
	}
	/**
	 * @param opt_taqman_pair_end_alignment_score the oPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE to set
	 */
	public void setOPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE(
			int opt_taqman_pair_end_alignment_score) {
		OPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE = opt_taqman_pair_end_alignment_score;
	}
	/**
	 * @return the pRIMER_DELTA_TM_WEIGHT
	 */
	public double getPRIMER_DELTA_TM_WEIGHT() {
		return PRIMER_DELTA_TM_WEIGHT;
	}
	/**
	 * @param primer_delta_tm_weight the pRIMER_DELTA_TM_WEIGHT to set
	 */
	public void setPRIMER_DELTA_TM_WEIGHT(double primer_delta_tm_weight) {
		PRIMER_DELTA_TM_WEIGHT = primer_delta_tm_weight;
	}
	/**
	 * @return the pRIMER_DELTA_GC_WEIGHT
	 */
	public double getPRIMER_DELTA_GC_WEIGHT() {
		return PRIMER_DELTA_GC_WEIGHT;
	}
	/**
	 * @param primer_delta_gc_weight the pRIMER_DELTA_GC_WEIGHT to set
	 */
	public void setPRIMER_DELTA_GC_WEIGHT(double primer_delta_gc_weight) {
		PRIMER_DELTA_GC_WEIGHT = primer_delta_gc_weight;
	}
	/**
	 * @return the pRIMER_DELTA_LENGTH_WEIGHT
	 */
	public double getPRIMER_DELTA_LENGTH_WEIGHT() {
		return PRIMER_DELTA_LENGTH_WEIGHT;
	}
	/**
	 * @param primer_delta_length_weight the pRIMER_DELTA_LENGTH_WEIGHT to set
	 */
	public void setPRIMER_DELTA_LENGTH_WEIGHT(double primer_delta_length_weight) {
		PRIMER_DELTA_LENGTH_WEIGHT = primer_delta_length_weight;
	}
	/**
	 * @return the pRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT
	 */
	public double getPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT() {
		return PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT;
	}
	/**
	 * @param primer_delta_distance_to_rss_weight the pRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT to set
	 */
	public void setPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT(
			double primer_delta_distance_to_rss_weight) {
		PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT = primer_delta_distance_to_rss_weight;
	}
	/**
	 * @return the sELF_ALIGNMENT_WEIGHT
	 */
	public double getSELF_ALIGNMENT_WEIGHT() {
		return SELF_ALIGNMENT_WEIGHT;
	}
	/**
	 * @param self_alignment_weight the sELF_ALIGNMENT_WEIGHT to set
	 */
	public void setSELF_ALIGNMENT_WEIGHT(double self_alignment_weight) {
		SELF_ALIGNMENT_WEIGHT = self_alignment_weight;
	}
	/**
	 * @return the sELF_END_ALIGNMENT_WEIGHT
	 */
	public double getSELF_END_ALIGNMENT_WEIGHT() {
		return SELF_END_ALIGNMENT_WEIGHT;
	}
	/**
	 * @param self_end_alignment_weight the sELF_END_ALIGNMENT_WEIGHT to set
	 */
	public void setSELF_END_ALIGNMENT_WEIGHT(double self_end_alignment_weight) {
		SELF_END_ALIGNMENT_WEIGHT = self_end_alignment_weight;
	}
	/**
	 * @return the pAIR_ALIGNMENT_WEIGHT
	 */
	public double getPAIR_ALIGNMENT_WEIGHT() {
		return PAIR_ALIGNMENT_WEIGHT;
	}
	/**
	 * @param pair_alignment_weight the pAIR_ALIGNMENT_WEIGHT to set
	 */
	public void setPAIR_ALIGNMENT_WEIGHT(double pair_alignment_weight) {
		PAIR_ALIGNMENT_WEIGHT = pair_alignment_weight;
	}
	/**
	 * @return the pAIR_END_ALIGNMENT_WEIGHT
	 */
	public double getPAIR_END_ALIGNMENT_WEIGHT() {
		return PAIR_END_ALIGNMENT_WEIGHT;
	}
	/**
	 * @param pair_end_alignment_weight the pAIR_END_ALIGNMENT_WEIGHT to set
	 */
	public void setPAIR_END_ALIGNMENT_WEIGHT(double pair_end_alignment_weight) {
		PAIR_END_ALIGNMENT_WEIGHT = pair_end_alignment_weight;
	}
	/**
	 * @return the pRIMER_FALSE_POSITIVES_WEIGHT
	 */
	public double getPRIMER_FALSE_POSITIVES_WEIGHT() {
		return PRIMER_FALSE_POSITIVES_WEIGHT;
	}
	/**
	 * @param primer_false_positives_weight the pRIMER_FALSE_POSITIVES_WEIGHT to set
	 */
	public void setPRIMER_FALSE_POSITIVES_WEIGHT(
			double primer_false_positives_weight) {
		PRIMER_FALSE_POSITIVES_WEIGHT = primer_false_positives_weight;
	}
	/**
	 * @return the pRIMER_PAIR_HOMOGENITY_WEIGHT
	 */
	public double getPRIMER_PAIR_HOMOGENITY_WEIGHT() {
		return PRIMER_PAIR_HOMOGENITY_WEIGHT;
	}
	/**
	 * @param primer_pair_homogenity_weight the pRIMER_PAIR_HOMOGENITY_WEIGHT to set
	 */
	public void setPRIMER_PAIR_HOMOGENITY_WEIGHT(
			double primer_pair_homogenity_weight) {
		PRIMER_PAIR_HOMOGENITY_WEIGHT = primer_pair_homogenity_weight;
	}
	/**
	 * @return the pRIMER_PAIR_DOPT_WEIGHT
	 */
	public double getPRIMER_PAIR_DOPT_WEIGHT() {
		return PRIMER_PAIR_DOPT_WEIGHT;
	}
	/**
	 * @param primer_pair_dopt_weight the pRIMER_PAIR_DOPT_WEIGHT to set
	 */
	public void setPRIMER_PAIR_DOPT_WEIGHT(double primer_pair_dopt_weight) {
		PRIMER_PAIR_DOPT_WEIGHT = primer_pair_dopt_weight;
	}
	
	public Pattern getSINGLE_NUCLEOTIDE_REPEAT_PATTERN() {
		return SINGLE_NUCLEOTIDE_REPEAT_PATTERN;
	}

	public String getEnzymeName() {
		return this.enzyme.getName();
	}

	public Pattern getEnzymeForwardPattern() {
		return Pattern.compile(this.enzyme.getForwardRegex());
	}

	public Pattern getEnzymeReversePattern() {
		return Pattern.compile(this.enzyme.getReverseRegex());
	}

	public int getMIN_MISPRIMING_TM_DIFFERENCE() {
		return MIN_MISPRIMING_TM_DIFFERENCE;
	}

	public void setMIN_MISPRIMING_TM_DIFFERENCE(int min_mispriming_tm_difference) {
		MIN_MISPRIMING_TM_DIFFERENCE = min_mispriming_tm_difference;
	}

	public boolean isComputeScanningStatistics() {
		return computeScanningStatistics;
	}

	public boolean isComputePickingStatistics() {
		return computePickingStatistics;
	}

	public boolean isCheckInterPairAlignments() {
		return checkInterPairAlignments;
	}

	public boolean isCheckIntraPairAlignments() {
		return checkIntraPairAlignments;
	}

	/**
	 * @return the targetOrganism
	 */
	public String getTargetOrganism() {
		return targetOrganism.toString();
	}

	/**
	 * @param targetOrganism the targetOrganism to set
	 */
	public void setTargetOrganism(Enum<TargetOrganisms> targetOrganism) {
		this.targetOrganism = targetOrganism;
	}

	/**
	 * @return the misprimingCheck
	 */
	public PrimerMisprimingCheck getPrimerMisprimingCheck() {
		return primerMisprimingCheck;
	}

	/**
	 * @param misprimingCheck the misprimingCheck to set
	 */
	public void setPrimerMisprimingCheck(PrimerMisprimingCheck misprimingCheck) {
		this.primerMisprimingCheck = misprimingCheck;
	}
	
	public void setConfigFile(Properties configFile){
		this.configFile = configFile;
	}
	
	public String getValueFromConfigFile(String key){
		return this.configFile.getProperty(key);
	}
	
	public PrimerSearchStatistics getSearchStat(){
		return this.searchStat;
	}
	
	public PrimerPairPickingStatistics getPickingStat(){
		return this.pickingStat;
	}

	/**
	 * @return the pickTaqManProbe
	 */
	public boolean isPickTaqManProbe() {
		return pickTaqManProbe;
	}

	/**
	 * @param pickTaqManProbe the pickTaqManProbe to set
	 */
	public void setPickTaqManProbe(boolean pickTaqManProbe) {
		this.pickTaqManProbe = pickTaqManProbe;
	}
	
	/**
	 * Returns all parameters for logging purposes.
	 * 
	 * @return a string containing all parameters for logging purposes
	 */
	public String getLogString(){
		return "PARAMS:\torganism: " + this.targetOrganism.toString() + 
		"\tpairs: " + this.numPrimers + 
		"\tEnzyme: " + this.enzyme.getName() +
		"\tinclude probes: " + this.pickTaqManProbe +
		"\tPRIMERS: " + 
		"\tlength: " + this.MIN_PRIMER_LENGTH + "/ " + this.OPT_PRIMER_LENGTH + "/ " + this.MAX_PRIMER_LENGTH + 
		"\tTM: " + this.MIN_TM + "/ " + this.OPT_TM + "/ " + this.MAX_TM + 
		"\tmax TM diff: " + this.MAX_PRIMER_TM_DIFFERENCE + 
		"\tSA/ SEA: " + this.MAX_PRIMER_SELF_ALIGNMENT_SCORE + "/ " + this.MAX_PRIMER_SELF_END_ALIGNMENT_SCORE + 
		"\tPA/ PEA: " + this.MAX_PRIMER_PAIR_ALIGNMENT_SCORE + "/ " + this.MAX_PRIMER_PAIR_END_ALIGNMENT_SCORE +
		"\tPROBES: " +
		"\tlength: " + this.TAQMAN_MIN_PRIMER_LENGTH + "/ " + this.TAQMAN_OPT_PRIMER_LENGTH + "/ " + this.TAQMAN_MAX_PRIMER_LENGTH + 
		"\tTM: " + this.TAQMAN_MIN_TM + "/ " + this.TAQMAN_OPT_TM + "/ " + this.TAQMAN_MAX_TM + 
		"\tmin TAQMAN TM diff: " + this.MIN_TAQMAN_TM_DIFFERENCE + 
		"\tSA/ SEA: " + this.MAX_TAQMAN_SELF_ALIGNMENT_SCORE + "/ " + this.MAX_TAQMAN_SELF_END_ALIGNMENT_SCORE + 
		"\tPA/ PEA: " + this.MAX_TAQMAN_PAIR_ALIGNMENT_SCORE + "/ " + this.MAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE +
		"\tamplicon length: " + this.MIN_AMPLICON_LENGTH + "/ " + this.OPT_AMPLICON_LENGTH + "/ " + this.MAX_AMPLICON_LENGTH +
		"\tsafe FP length: " + this.SAFE_FALSE_POSITIVE_AMPLICON_LENGTH + 
		"\tTM weight: " + this.PRIMER_DELTA_TM_WEIGHT + 
		"\tGC weight: " + this.PRIMER_DELTA_GC_WEIGHT + 
		"\tlength weight: " + this.PRIMER_DELTA_LENGTH_WEIGHT +
		"\tdRSS weight: " + this.PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT + 
		"\tSA weight: " + this.SELF_ALIGNMENT_WEIGHT + 
		"\tSEA weight: " + this.SELF_END_ALIGNMENT_WEIGHT + 
		"\tPA weight: " + this.PAIR_ALIGNMENT_WEIGHT + 
		"\tPEA weight: " + this.PAIR_END_ALIGNMENT_WEIGHT + 
		"\tFP weight: " + this.PRIMER_FALSE_POSITIVES_WEIGHT;
	}
}
