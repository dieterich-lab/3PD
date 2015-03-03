package primerDesign.algo;

import java.io.File;
import java.io.FileInputStream;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerAcceptanceLevel;
import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.PrimerSearchStatistics;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.SequenceRegionAlignment;
import primerDesign.dsc.indexStructures.TargetOrganisms;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheckDeserializer;
import primerDesign.util.Constants;
import primerDesign.util.DuplicateUseOfRestrictionSiteException;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import primerDesign.web.PrimerDesignWebProperites;
import weka.core.FastVector;
import cern.colt.list.ObjectArrayList;

/**
 * This class enumerates all valid primers in a given interval.
 * 
 * @author Sebastian Fr�hler
 *
 */
public class PrimerSearch {
	
	private static final Pattern repetitiveElementPattern = Pattern.compile(".*" + Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase() + ".*");
	private static final Pattern ambiguityCodePattern = Pattern.compile(".*[MRWSYKVHDB].*");
	private Pattern singleNucleotideRepeatPattern;
	private static final Pattern newGcClampPattern = Pattern.compile("[GC]");
	private Pattern restrictionEnzymeForwardPattern;
	private Pattern restrictionEnzymeReversePattern;
	private static RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
	private static Sequence rssScanRegion;
	private static int primerCount = 0;
	private static int primerIndexCount = 0;
	private PrimerSearchStatistics stat = new PrimerSearchStatistics();
	private boolean doStat = true;
	
	/**
	 * Enumerates all 'potential' primers and returns all valid primers.
	 * 
	 * This method computes 'one column in the primer matrix'. Forward and reverse primers are checked for misprimings and safe distance to restriction sites at mispriming positions.
	 * Hybridization probes are not tested for misprimings since those misprimings do not affect the experiment (w.r.t. false positive fluorescence emission)
	 * 
	 * @param sequence the sequence the primer can be placed in
	 * @param primerType the type of the primer
	 * @param regionStart the start of 'sequence' w.r.t. the genomic sequence
	 * @param returnOrdered whether valid primers should be returned ordered by their distance to the '(virtual) optimal' primer
	 * @param fivePrimeDistToRSS the five prime distance of the scan region to the restriction site
	 * @param enzyme the restriction enzyme which generated the restriction site adjacent to the scan region
	 * 
	 * @return the set of all valid primers for sequence 'sequence'
	 */
	public Primer[] enumeratePrimers(String sequence, PrimerTypes primerType, int minlength, int maxlength, int regionStart, boolean returnOrdered, int fivePrimeDistToRSS, RestrictionSite restrictionSite, RestrictionEnzyme enzyme, PrimerSearchParameters searchParams){
		// if valid primers for this sequence have been precomputed, return these
//		if(enumerated_primers_hash.containsKey(sequence)){
//			return enumerated_primers_hash.get(sequence);
//		}else{
		//System.out.println("Enumerate primers for sequence " + sequence);
		
		if(this.doStat) this.stat.incEnumeratePrimersCount();
		
		this.singleNucleotideRepeatPattern = searchParams.getSINGLE_NUCLEOTIDE_REPEAT_PATTERN();
		this.restrictionEnzymeForwardPattern = searchParams.getEnzymeForwardPattern();
		this.restrictionEnzymeReversePattern = searchParams.getEnzymeReversePattern();
		
			PrimerMisprimingCheck misprimingCheck = searchParams.getPrimerMisprimingCheck();
			Primer currentPrimer; 
			ObjectArrayList primers = new ObjectArrayList();
			//Vector<Primer> primers = new Vector<Primer>();
			
			//RSSDPAligner scanRegionAligner = new RSSDPAligner(sequence, regionStart);
			int relativePosition;
			if(primerType.equals(PrimerTypes.forwardPrimer)){
				restrictionSite.setForwardScanSequence(sequence.toCharArray());
			}
			else if(primerType.equals(PrimerTypes.reversePrimer)){
				restrictionSite.setReverseScanSequence(sequence.toCharArray());
			}
			else if(primerType.equals(PrimerTypes.hybridizationProbe)){
				restrictionSite.setProbeScanSequence(sequence.toCharArray());
			}
			else{
				throw new IllegalArgumentException("Unsupported primer type!");
			}
			
			if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Scan - enumerating primers in sequence: " + sequence);
			
			SequenceRegionAlignment alignment = SequenceRegionAligner.alignSequenceRegions(sequence.toCharArray(), sequence.toCharArray(), searchParams.getA_t_basepair_score(), searchParams.getG_c_basepair_score());
			PrimerAlignmentScores scores;
			String primerSequence;
			
			for(int i=minlength; i<=maxlength ; i++){
				//for(int j=0; j<sequence.length()-maxlength; j++){
				for(int j=0; j<=sequence.length()-minlength; j++){
					if(j+i > sequence.length()) continue;
					// if primer does not contain masked basepairs = repetitive dna -> primers containing repetitive dna are prone to mispriming!
					// AND if primer sequence does NOT contain a restriction site for the specified restriction enzyme (-> would prevent full-length priming iff cut by restriction enzyme)
					primerSequence = sequence.substring(j, j+i).toUpperCase();
					
					if(this.doStat) this.stat.incPrimerCount();
					PrimerSearch.primerCount++;
					
					// skip primers containing ambiguity characters as bases
					if(ambiguityCodePattern.matcher(primerSequence).matches()){
						this.stat.incPrimerAmbiguityCode();
						continue;
					}
					
					if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Scan - candidate primer: " + primerSequence + " of length " + primerSequence.length());
					
					// parameters: http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
					// Poly-X exclude: .*([A]{5,}|[T]{5,}|[G]{5,}|[C]{5,}).*
					if(singleNucleotideRepeatPattern.matcher(primerSequence).matches()){
						if(this.doStat) this.stat.incSingleNucleotideRepeatPattern();
						if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Reject - repetitive single nucleotides: " + primerSequence);
						continue;
					}
					if(!PrimerSearch.matchesGCClampPattern(primerSequence, searchParams)){
						if(this.doStat) this.stat.incGCClampCount();
						if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Reject - NO GC-Clamp: " + primerSequence);
						continue;
					}
					
					//System.out.println("Current primer: " + primerSequence);
					
					if(!(repetitiveElementPattern.matcher(primerSequence).matches()) && !(restrictionEnzymeForwardPattern.matcher(primerSequence).matches()) && !(restrictionEnzymeReversePattern.matcher(primerSequence).matches())){
						if(primerType.equals(PrimerTypes.forwardPrimer)) relativePosition = regionStart + j;
						// due to 3C experiment setup, hybridization probes are designed to located at the reverse strand 5' of the restriction site
						else if(primerType.equals(PrimerTypes.reversePrimer) || primerType.equals(PrimerTypes.hybridizationProbe)) relativePosition = regionStart - j;
						else throw new IllegalArgumentException("Unsupported primer type! " + primerType);
						
						scores = alignment.getGlobalAlignmentValues(j, primerSequence.length(), j, primerSequence.length());
						
						currentPrimer = new Primer(primerSequence, primerType, relativePosition, scores.getPairScore(), scores.getPairEndScore(), searchParams); //new Primer(primerSequence, primerType, scanRegionAligner, relativePosition); //new Primer(primerSequence, primerType);
						currentPrimer.setRestrictionSite(restrictionSite);
						currentPrimer.setPositionInScanSequence(j);
						int currentPrimerLength = currentPrimer.getSequence().length();
						
						// debug trace
						//if(primerType.equals(PrimerTypes.hybridizationProbe)) System.out.println("Seq " + currentPrimer.getSequence() + " Pos " + currentPrimer.getRelativePosition() + " TM: " + currentPrimer.getMeltingTemp() + " GC " + currentPrimer.getGcContent());
						
						if(primerType.equals(PrimerTypes.forwardPrimer) || primerType.equals(PrimerTypes.reversePrimer)) currentPrimer.setDistanceToRSS(fivePrimeDistToRSS - j);
						else if(primerType.equals(PrimerTypes.hybridizationProbe)) currentPrimer.setDistanceToRSS(fivePrimeDistToRSS + j);
						else throw new IllegalArgumentException("Unsupported primer type! " + primerType);
						
						if((primerType.equals(PrimerTypes.forwardPrimer) || primerType.equals(PrimerTypes.reversePrimer)) && isValidPrimer(currentPrimer, searchParams)){
							//int maxMispriming = 0;
							boolean primerMisprimingDoesNOTAffectExperiment = true;
							//int five_prime_start = currentPrimerLength - MyExtendedMath.floor(currentPrimerLength * Constants.SAFE_MISPRIMING_SEQUENCE_PRECENT_CUTOFF) -1;
							// scan primer mispriming property in increasing order, starting form three prime end - offset 'five_prime_start
							// primer mismatch counts decrease monotonically with increasing length of 'primerEnd'
//							for(int k= five_prime_start; k > 0; k--){
//								String primerEnd = currentPrimer.getSequence().substring(k, currentPrimerLength);
//								float hits = 0;
//								if(backgroundIndex.includesScanRegion()){
//									hits = backgroundIndex.searchNbMatchesInIndex(primerEnd);
//								}else{
//									hits = primerScanRegionIndex.searchNbMatchesInIndex(primerEnd) + backgroundIndex.searchNbMatchesInIndex(primerEnd);
//								}								
//								if(hits == 1){
//									currentPrimer.setFalsePositiveMatches(0);
//									break;
//								}else if(hits > 1){
//									hits--; // besides a TP hit, there are hits-1 FP hits in the sequence!
//									if(hits > currentPrimer.getFalsePositiveMatches()) currentPrimer.setFalsePositiveMatches((int) hits);
//									if(hits > maxMispriming) maxMispriming = (int) hits;
//									hits++; // re-increment hits for proper break condition
//								}
//								if(hits == 1) break; // if a primer-end can only be found once, each longer stretch of that same primer end can also only be found at max once!
//							}
//							// if this primer can produce spurious amplicons besides the desired amplicon, reject this primer
//							if(!PrimerMispriming.hasSafeDistanceToNextRSS(currentPrimer.getSequence(), primerType, enzyme, backgroundIndex)){
//								maxMispriming = Integer.MAX_VALUE;
//								primerMisprimingAffectsExperiment = true;
//								currentPrimer.setFalsePositiveMatches(maxMispriming);
//								//break;
//							}
							
//							ObjectArrayList hits;
//							if(currentPrimer.getPrimerType().equals(PrimerTypes.forwardPrimer)){
//								//hits = ApproxHitSearcher.findHits(currentPrimer.getSequence(), backgroundIndex, searchParams.getSEED_WORD_SIZE(), searchParams.getSEED_STEP_SIZE(), searchParams.getTHREE_PRIME_DP_LENGTH(), searchParams.getTHREE_PRIME_DP_THRESHOLD(), searchParams.getWHOLE_DP_THRESHOLD());
//								hits = backgroundIndex.findHitPositions(currentPrimer.getSequence().substring(currentPrimer.getLength() - searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
//								if(this.doStat) this.stat.incPrimerIndexCount();
//								//System.err.println("Found " + hits.size() + " hits for forward primer " + currentPrimer.getSequence());
//							}
//							else{
//								//hits = ApproxHitSearcher.findHits(SeqTools.complementDNA(currentPrimer.getSequence().toCharArray()), backgroundIndex, Constants.SEED_WORD_SIZE, Constants.SEED_STEP_SIZE, Constants.THREE_PRIME_DP_LENGTH, Constants.THREE_PRIME_DP_THRESHOLD, Constants.WHOLE_DP_THRESHOLD);
//								hits = backgroundIndex.findHitPositions(currentPrimer.getSequence().substring(currentPrimer.getLength()  - searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
//								if(this.doStat) this.stat.incPrimerIndexCount();
//								//System.err.println("Found " + hits.size() + " hits for reverse primer " + currentPrimer.getSequence());
//							}
//							//ObjectArrayList hits = new ObjectArrayList();
//							if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Found " + hits.size() + " hits for primer " + currentPrimer.getSequence() + " in index!");
//							primerMisprimingDoesNOTAffectExperiment = PrimerMispriming.hasSafeDistanceToNextRSS(currentPrimer, enzyme, hits, searchParams);
//							maxMispriming = hits.size() - 1;
							
							if(Constants.doEarlyMMScan){
								primerMisprimingDoesNOTAffectExperiment = !misprimingCheck.hasMisprimings(currentPrimer);
								
								if(primerMisprimingDoesNOTAffectExperiment) currentPrimer.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
								else currentPrimer.setAcceptanceLevel(PrimerAcceptanceLevel.UNACCEPTABLE);
								
								if(this.doStat) this.stat.incPrimerIndexCount();
								
								primerIndexCount++;
							}
							
							//if(maxMispriming <= Constants.MAX_PRIMER_MISPRIMING_CUTOFF || !primerMisprimingAffectsExperiment)	primers.add(currentPrimer);			
							if(primerMisprimingDoesNOTAffectExperiment){
								primers.add(currentPrimer);
								if(this.doStat) this.stat.incNoMispriming();
							}
							else{
								if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Reject - Primer mispriming does affect experiment: " + !primerMisprimingDoesNOTAffectExperiment + " found FP hits! for primer " + currentPrimer.toString());
								if(this.doStat) this.stat.incMispriming();
							}
						}
						else if(primerType.equals(PrimerTypes.hybridizationProbe) && isValidHybridizationProbe(currentPrimer, searchParams)){
							currentPrimer.setAcceptanceLevel(PrimerAcceptanceLevel.ACCEPTABLE);
							primers.add(currentPrimer);
						}
						
						else if(searchParams.isPRINT_DEBUG_LOG() && (primerType.equals(PrimerTypes.forwardPrimer) || primerType.equals(PrimerTypes.reversePrimer)) && !isValidPrimer(currentPrimer, searchParams)) System.err.println("Reject - invalid primer: " + currentPrimer.toString());
						else if(searchParams.isPRINT_DEBUG_LOG() && primerType.equals(PrimerTypes.hybridizationProbe) && !isValidHybridizationProbe(currentPrimer, searchParams)) System.err.println("Reject - invalid probe: " + currentPrimer.toString());
					}
					else{
						if(searchParams.isPRINT_DEBUG_LOG() && repetitiveElementPattern.matcher(primerSequence).matches()) System.err.println("Reject - masked basepairs (repeats): " + primerSequence);
						else if(searchParams.isPRINT_DEBUG_LOG() && restrictionEnzymeForwardPattern.matcher(primerSequence).matches()) System.err.println("Reject - restriction site in forward sequence found: " + primerSequence);
						else if(searchParams.isPRINT_DEBUG_LOG() && restrictionEnzymeReversePattern.matcher(primerSequence).matches()) System.err.println("Reject - restriction site in reverse sequence found: " + primerSequence);
						if(this.doStat){
							if(repetitiveElementPattern.matcher(primerSequence).matches()) this.stat.incRepetitiveElementPattern();
							if(restrictionEnzymeForwardPattern.matcher(primerSequence).matches() || restrictionEnzymeReversePattern.matcher(primerSequence).matches()) this.stat.incRestrictionEnzymePattern();
						}
					}
				}
			}
			// include primer mispriming verification as batch after all primers have been enumerated
			
			
			if(returnOrdered){
				// sort primers in ascending distance to '(virtual) optimal' primer
				//Collections.sort(primers);
				primers.quickSortFromTo(0, primers.size()-1);
			}
			Primer[] result = new Primer[primers.size()];
			if(primers.size() > 0){
				primers.toArray(result);
			}
			return result;
		//}
	}
	
	/**
	 * Matches the primer end to check for GC-clamp criterion.
	 * 
	 * The occurrence of a specific number of G and/or C bases at the primer end is supporsed to lead
	 * to a more specific binding of the primer end to the target dna.
	 * The length of the primer end is: Constants.GC_CLAMP_LENGTH.
	 * 
	 * @param primer the primer to scan for GC-clamp
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return true iff the primer end contains between Constants.MIN_GC_CLAMP (inclusive) and Constants.MAX_GC_CLAMP (inclusive)
	 */
	public static boolean matchesGCClampPattern(String primer, PrimerSearchParameters searchParams){
		if(primer.length() < searchParams.getGC_CLAMP_LENGTH()) throw new IllegalArgumentException("GC_CLAMP_LENGTH has to be <= than primer length!");
		Matcher matcher = newGcClampPattern.matcher(primer.substring(primer.length() - searchParams.getGC_CLAMP_LENGTH()));
		
		int pos;
		boolean found = matcher.find();
		if(!found && searchParams.getMIN_GC_CLAMP() == 0) return true;
		else if(found){
			pos = matcher.start();
			int matches = 1;
			while(matcher.find(++pos)){
				pos = matcher.start();
				matches++;
			}
			return matches >= searchParams.getMIN_GC_CLAMP() && matches <= searchParams.getMAX_GC_CLAMP();
		}
		else return false;
	}
	
	/**
	 * Determines whether a primer is acceptable according to the preconditions specified.
	 * 
	 * @param primer the primer
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return true: iff all preconditions are met, false: else
	 */
	private boolean isValidPrimer(Primer primer, PrimerSearchParameters searchParams){
		if(primer.getGcContent() >= searchParams.getMIN_GC() &&
				primer.getGcContent() <= searchParams.getMAX_GC() &&
				primer.getMeltingTemp() >= searchParams.getMIN_TM() &&
				primer.getMeltingTemp() <= searchParams.getMAX_TM() &&
				primer.getSelfAlignmentScore() <= searchParams.getMAX_PRIMER_SELF_ALIGNMENT_SCORE() &&
				primer.getSelfEndAlignmentScore() <= searchParams.getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE()
				){
			if(this.doStat) this.stat.incValidPrimer();
			return true;
		}else{
			if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Reject - Invalid primer: " + primer.toString());
			if(this.doStat){
				this.stat.incInvalidPrimer();
				this.stat.evaluatePrimer(primer, searchParams);
			}
			return false;
		}
	}
	
	/**
	 * Determines whether a hybridization probe is valid according to the preconditions specified.
	 * 
	 * @param primer the hybridization probe
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return true iff hybridization probe 'primer' is valid
	 */
	private boolean isValidHybridizationProbe(Primer primer, PrimerSearchParameters searchParams){
		if(primer.getGcContent() >= searchParams.getTAQMAN_MIN_GC() &&
				primer.getGcContent() <= searchParams.getTAQMAN_MAX_GC() &&
				primer.getMeltingTemp() >= searchParams.getTAQMAN_MIN_TM() &&
				primer.getMeltingTemp() <= searchParams.getTAQMAN_MAX_TM() &&
				primer.getSelfAlignmentScore() <= searchParams.getMAX_TAQMAN_SELF_ALIGNMENT_SCORE() &&
				primer.getSelfEndAlignmentScore() <= searchParams.getMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE() &&
				primer.getSequence().charAt(0) != 'G'
				){
			if(this.doStat) this.stat.incValidTaqManPrimer();
			return true;
		}else{
			if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Reject - Invalid TaqMan probe: " + primer.toString());
			if(this.doStat){
				this.stat.incInvalidTaqManPrimer();
				this.stat.evaluateProbe(primer, searchParams);
			}
			return false;
		}
	}
	
	/**
	 * Naively searches the search space for 'nbOfPrimers' best primers for sequence 'sequence' and restriction sites of enzyme 'enzyme'.
	 * 
	 * A restriction enzyme (endo-nuclease) cuts a sequence generating several restriction sites. A set of 'nbOfPrimers' most homogenuously distributed
	 * restriction sites of those are first searched. Thereafter, primers are to be searched in the upstream region
	 * of each restriction site and the most homogenuous set of these primers is returned.
	 * 
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return a set of sequence regions with valid primers (according to the constraints specified)
	 */
	public RestrictionSite[] naivePrimerSearch(PrimerSearchParameters searchParams){		
		this.stat = searchParams.getSearchStat();
		//this.singleNucleotideRepeatPattern = Pattern.compile(".*([A]{" + searchParams.getMAX_POLY_X_LENGTH() +  ",}|[T]{" + searchParams.getMAX_POLY_X_LENGTH() +  ",}|[G]{" + searchParams.getMAX_POLY_X_LENGTH() +  ",}|[C]{" + searchParams.getMAX_POLY_X_LENGTH() +  ",}).*");
		this.singleNucleotideRepeatPattern = searchParams.getSINGLE_NUCLEOTIDE_REPEAT_PATTERN();
		
		int numberOfPrimers = searchParams.getNumPrimers();
		RestrictionEnzyme enzyme = searchParams.getEnzyme();
//		DNASequenceIndex primerScanRegionIndex = searchParams.getScanRegionIndex();
//		primerDesign.dsc.indexStructures.DNASequenceIndex backgroundIndex = searchParams.getBackgroundIndex();
		PrimerMisprimingCheck misprimingCheck = searchParams.getPrimerMisprimingCheck();
		
		// scan for most homogenuous restriction sites
		boolean returnOrderedPrimerLists = false;
		RestrictionSiteSearch search = new RestrictionSiteSearch();
		FastVector optimalRSSs = new FastVector();
		optimalRSSs.appendElements(search.getMostHomogenuousRSSs(enzyme, numberOfPrimers, searchParams));
		
		RestrictionEnzymeManager.register(enzyme, new TreeSet());
		PrimerSearch.mapper.clearEnzymes();
		PrimerSearch.mapper.addEnzyme(enzyme);
		
		String sequence;
		
		// compute primers for each restriction site
		for(int i=0; i<optimalRSSs.size(); i++){
			System.out.println("Site: " + (i+1) + "/" + optimalRSSs.size());
			RestrictionSite site = (RestrictionSite) optimalRSSs.elementAt(i);
			int position = site.getPosition();
			// compute start and end for 5' scan region for primer search
			int forwardStart =  Math.max(0,position - searchParams.getMAX_AMPLICON_LENGTH()/2);
			int forwardEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
			
			// compute start and end for 3' scan region for primer search
			int reverseStart = Math.min(position + searchParams.getMAX_AMPLICON_LENGTH()/2, site.getSequenceRegion().getContigSequenceLength() - 1);
			int reverseEnd = position + (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
			
			// compute start and end for 5' scan region for hybridization probe search
			int hybProbeStart = position - 1;
			int hybProbeEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
			if(searchParams.isPickTaqManProbe() && searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH() < searchParams.getTAQMAN_MAX_PRIMER_LENGTH()) throw new IllegalArgumentException("MIN_AMPLICON_LENGTH/2 - MAX_PRIMER_LENGTH must be > TAQMAN_MAX_PRIMER_LENGTH! Scan regions for forward primer and HybProbe must NOT overlap!");
			
			assert(forwardEnd+1 > forwardStart);
			assert(reverseStart+1 > reverseEnd);
			if(searchParams.isPickTaqManProbe()) assert(hybProbeStart+1 > hybProbeEnd);
			
			sequence = site.getSequenceRegion().getContigSequence();
			
			if(forwardStart < 0 || forwardEnd+1<=forwardStart || reverseEnd >= sequence.length()) continue;
			
			String forwardScanRegion = sequence.substring(forwardStart, forwardEnd+1);
			String reverseScanRegion = SeqTools.revcompDNA(sequence.substring(reverseEnd, reverseStart+1).toCharArray());
			String hybProbeScanRegion = "";
			if(searchParams.isPickTaqManProbe()) hybProbeScanRegion = SeqTools.revcompDNA(sequence.substring(hybProbeEnd, hybProbeStart+1).toCharArray());
			
			Primer[] validUpstreamPrimers = new Primer[0];
			Primer[] validDownstreamPrimers = new Primer[0];
			Primer[] validHybridizationProbes = new Primer[0];
			
			String wholeAmplicon = sequence.substring(position - searchParams.getMAX_AMPLICON_LENGTH()/2, position + searchParams.getMAX_AMPLICON_LENGTH()/2);
			
			// if another restriction site for enzyme was found in this amplicon, reject this amplicon and 'search remaining restriction sites in this region'
			if(doesNotContainAnotherRSS(wholeAmplicon)){
				// compute and store all valid 5' primers in sequence 'scanRegion' of current restriction site
				validUpstreamPrimers = this.enumeratePrimers(forwardScanRegion, PrimerTypes.forwardPrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), forwardStart, returnOrderedPrimerLists, position - forwardStart, site, enzyme, searchParams);
				if(validUpstreamPrimers.length != 0) validDownstreamPrimers = this.enumeratePrimers(reverseScanRegion, PrimerTypes.reversePrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), reverseStart, returnOrderedPrimerLists, -position + reverseStart, site, enzyme, searchParams);
				else if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Empty upstream primer list");
				if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && searchParams.isPickTaqManProbe()) validHybridizationProbes = this.enumeratePrimers(hybProbeScanRegion, PrimerTypes.hybridizationProbe, searchParams.getTAQMAN_MIN_PRIMER_LENGTH(), searchParams.getTAQMAN_MAX_PRIMER_LENGTH(), hybProbeStart, returnOrderedPrimerLists, 1, site, enzyme, searchParams);
				else if(searchParams.isPRINT_DEBUG_LOG() && validUpstreamPrimers.length == 0) System.out.println("Empty upstream primer list!"); 
				else if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Empty downstream primer list");
				if(this.doStat){
					if(validUpstreamPrimers.length == 0) this.stat.incEmptyUpstreamPrimersList();
					else this.stat.incNonemptyUpstreamPrimersList();
					if(validDownstreamPrimers.length == 0) this.stat.incEmptyDownstreamPrimersList();
					else this.stat.incNonEmptyDownstreamPrimersList();
					if(searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0) this.stat.incEmptyTaqManPrimersList();
					else if(searchParams.isPickTaqManProbe()) this.stat.incNonEmptyTaqManPrimersList();
				}
			}
			else{
				if(this.doStat) this.stat.incAnotherRSSinAmpliconCount();
				if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Found another RSS in this amplicon!");
			}
			
			if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && (!searchParams.isPickTaqManProbe() || validHybridizationProbes.length != 0)){
				site.setValidUpstreamPrimers(validUpstreamPrimers); 
				site.setValidDownstreamPrimers(validDownstreamPrimers);
				if(searchParams.isPickTaqManProbe()) site.setValidTaqManProbes(validHybridizationProbes);
				site.wasScannedForPrimers(true);
				if(this.doStat) this.stat.incNonEmptyUpDpwnTaqManList();
			}else{
				if(this.doStat) this.stat.incEmptyUpDownTaqManList();
				// if no acceptable 5', 3' or TaqMan primers can be found in this 'scanRegion' (at current restriction site), scan next best restriction site's sequence in current sequence region
				Iterator allRSSIterator = site.getSequenceRegion().getRestrictionSitesIterator(true);
				RestrictionSite nextBestSite;
				while((validUpstreamPrimers.length == 0 || validDownstreamPrimers.length == 0 || (searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0)) && allRSSIterator.hasNext()){
					nextBestSite = (RestrictionSite) allRSSIterator.next();
					
					// if this restriction site is not already used as optimal restriction site
					if(!optimalRSSs.contains(nextBestSite)){
						position = nextBestSite.getPosition();
						
						if(position - searchParams.getMAX_AMPLICON_LENGTH()/2 < 0 || position + searchParams.getMAX_AMPLICON_LENGTH()/2 >= sequence.length()) continue; // terminate search with empty result (last element!)
						wholeAmplicon = sequence.substring(position - searchParams.getMAX_AMPLICON_LENGTH()/2, position + searchParams.getMAX_AMPLICON_LENGTH()/2);
						
						if(doesNotContainAnotherRSS(wholeAmplicon)){
							// compute start and end for scan region for 5' primer search
							forwardStart = Math.max(0, position - searchParams.getMAX_AMPLICON_LENGTH()/2);
							forwardEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
							
							reverseStart = Math.min(position + searchParams.getMAX_AMPLICON_LENGTH()/2, sequence.length() - 1);
							reverseEnd = position + (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
							
							// compute start and end for 5' scan region for hybridization probe search
							hybProbeStart = position - 1;
							hybProbeEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
							
							if(forwardStart < 0 || forwardEnd >= sequence.length() || reverseStart < 0 || reverseEnd >= sequence.length() || hybProbeStart < 0 || hybProbeStart >= sequence.length()) break; // terminate search with empty result (last element!)
							if(searchParams.isPickTaqManProbe() && searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH() < searchParams.getTAQMAN_MAX_PRIMER_LENGTH()) throw new IllegalArgumentException("MIN_AMPLICON_LENGTH/2 - MAX_PRIMER_LENGTH must be > TAQMAN_MAX_PRIMER_LENGTH! Scan regions for forward primer and HybProbe must NOT overlap!");
							
							assert(forwardEnd+1 > forwardStart);
							assert(reverseStart+1 > reverseEnd);
							if(searchParams.isPickTaqManProbe()) assert(hybProbeStart+1 > hybProbeEnd);
							
							if(forwardStart < 0 || forwardEnd+1<=forwardStart || reverseEnd >= sequence.length()) continue;
							
							forwardScanRegion = sequence.substring(forwardStart, forwardEnd+1);
							reverseScanRegion = SeqTools.revcompDNA(sequence.substring(reverseEnd, reverseStart+1).toCharArray());
							if(searchParams.isPickTaqManProbe()) hybProbeScanRegion = SeqTools.revcompDNA(sequence.substring(hybProbeEnd, hybProbeStart+1).toCharArray());
							
							validUpstreamPrimers = this.enumeratePrimers(forwardScanRegion, PrimerTypes.forwardPrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), forwardStart, returnOrderedPrimerLists, position - forwardStart, nextBestSite, enzyme, searchParams);
							if(validUpstreamPrimers.length != 0) validDownstreamPrimers = this.enumeratePrimers(reverseScanRegion, PrimerTypes.reversePrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), reverseStart, returnOrderedPrimerLists, -position + reverseStart, nextBestSite, enzyme, searchParams);
							else{
								validDownstreamPrimers = new Primer[0];
								if(this.doStat) this.stat.incEmptyUpstreamPrimersList();
								if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Empty upstream primer list");
							}
							if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && searchParams.isPickTaqManProbe()) validHybridizationProbes = this.enumeratePrimers(hybProbeScanRegion, PrimerTypes.hybridizationProbe, searchParams.getTAQMAN_MIN_PRIMER_LENGTH(), searchParams.getTAQMAN_MAX_PRIMER_LENGTH(), hybProbeStart, returnOrderedPrimerLists, 1, nextBestSite, enzyme, searchParams);
							else{
								validHybridizationProbes = new Primer[0];
								if(this.doStat) this.stat.incEmptyDownstreamPrimersList();
								if(searchParams.isPRINT_DEBUG_LOG() && validUpstreamPrimers.length == 0) System.out.println("Empty upstream primer list!"); 
								else if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Empty downstream primer list");
							}
							
							if(this.doStat && (searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0)) this.stat.incEmptyTaqManPrimersList();
							if(searchParams.isPRINT_DEBUG_LOG() && searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0) System.err.println("Empty TaqMan probe list!");
							
							// remember new, optimal restriction site 'nextBestSite' and all of its valid primers instead of 'old', unsuitable restriction site 'site'
							if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && (!searchParams.isPickTaqManProbe() || validHybridizationProbes.length != 0)){
								if(this.doStat) this.stat.incNonEmptyUpDpwnTaqManList();
								nextBestSite.setValidUpstreamPrimers(validUpstreamPrimers);
								nextBestSite.setValidDownstreamPrimers(validDownstreamPrimers);
								if(searchParams.isPickTaqManProbe()) nextBestSite.setValidTaqManProbes(validHybridizationProbes);
								nextBestSite.wasScannedForPrimers(true);
								// replace 'site' by 'nextBestSite' (in result vector) and continue with next 'optimal restriction site'
								optimalRSSs.setElementAt(nextBestSite, optimalRSSs.indexOf(site));
								break;
							}
						}else{
							if(this.doStat) this.stat.incAnotherRSSinAmpliconCount();
							if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Found another RSS in this amplicon!");
						}
					}
					else throw new DuplicateUseOfRestrictionSiteException("Restriction site at position " + nextBestSite.getPosition() + " is used multiple times!");
				}
			}
			// base case if no primers can be found at all in this sequence region
			if(validUpstreamPrimers.length == 0) throw new EmptyResultSetException("No valid upstream primers could be found in Sequence " + i + " in SequenceRegion: " + site.getSequenceRegion().getSeqRegionStart() + "-" + site.getSequenceRegion().getSeqRegionEnd() + " (at any restriction site) please check your selection parameters/ change enzyme!");
			else if(validDownstreamPrimers.length == 0) throw new EmptyResultSetException("No valid downstream primers could be found in Sequence " + i + " inSequenceRegion: " + site.getSequenceRegion().getSeqRegionStart() + "-" + site.getSequenceRegion().getSeqRegionEnd() + " (at any restriction site) please check your selection parameters/ change enzyme!");
			else if(searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0) throw new EmptyResultSetException("No valid hybridization probes could be found in Sequence " + i + " in SequenceRegion: " + site.getSequenceRegion().getSeqRegionStart() + "-" + site.getSequenceRegion().getSeqRegionEnd() + " (at any restriction site) please check your selection parameters/ change enzyme!");
		}
		RestrictionSite[] result = new RestrictionSite[optimalRSSs.size()];
		System.arraycopy(optimalRSSs.toArray(), 0, result, 0, optimalRSSs.size());
		System.out.println("return");
		return result;
	}
	
	/**
	 * Checks whether the query sequence contains more than one instance of a specific restriction site.
	 * 
	 * If yes: no valid qPCR amplicon can be formed! So discard this query sequence!
	 * If no: accept this sequence
	 * 
	 * @param sequence the query sequence to scan for RSSs
	 * 
	 * @return true if query sequence only contains one RSS for enzyme specified by user 
	 */
	private boolean doesNotContainAnotherRSS(String sequence){
		try{
			PrimerSearch.rssScanRegion = DNATools.createDNASequence(sequence, "");
		}catch(IllegalSymbolException e){
			e.printStackTrace();
		}
		PrimerSearch.rssScanRegion = PrimerSearch.mapper.annotate(PrimerSearch.rssScanRegion);
		
		return PrimerSearch.rssScanRegion.countFeatures() == 1;
	}
	
	/**
	 * Sets the restriction enzyme pattern.
	 * 
	 * @param restrictionEnzymeForwardPattern the forward pattern
	 * @param restrictionEnzymeReversePattern the reverse pattern
	 */
	public void setEnzymePatterns(Pattern restrictionEnzymeForwardPattern, Pattern restrictionEnzymeReversePattern){
		this.restrictionEnzymeForwardPattern = restrictionEnzymeForwardPattern;
		this.restrictionEnzymeReversePattern = restrictionEnzymeReversePattern;
	}
	
	/**
	 * Refines a sequence region, returns an updated list of restriction sites with RSS at position 'refinePosition' being refined.
	 * 
	 * @param sites a list of restriction sites
	 * @param refinePosition the position of the list of restriction sites to be refined
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return an updated list of restriction sites with RSS at position 'refinePosition' being refined
	 */
	public RestrictionSite[] refineSeqRegion(RestrictionSite[] sites, int refinePosition, PrimerSearchParameters searchParams){
		assert(refinePosition >= 0 && refinePosition < sites.length);
		
		String sequence = sites[refinePosition].getSequenceRegion().getContigSequence();
		RestrictionEnzyme enzyme = searchParams.getEnzyme();
		
		Iterator<RestrictionSite> allRSSIterator = sites[refinePosition].getSequenceRegion().getCurrentRestrictionSitesIterator();
		if(allRSSIterator == null) allRSSIterator = sites[refinePosition].getSequenceRegion().getRestrictionSitesIterator(true);
		RestrictionSite nextBestSite = null;
		Primer[] validUpstreamPrimers = new Primer[0];
		Primer[] validDownstreamPrimers = new Primer[0];
		Primer[] validHybridizationProbes = new Primer[0];
		
		boolean returnOrderedPrimerLists = false;
		
		while((validUpstreamPrimers.length == 0 || validDownstreamPrimers.length == 0 || (searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0)) && allRSSIterator.hasNext()){
			nextBestSite = allRSSIterator.next();
			
			// if this restriction site is not already used as optimal restriction site
			if(!sites[refinePosition].equals(nextBestSite)){
				int position = nextBestSite.getPosition();
				
				if(position - searchParams.getMAX_AMPLICON_LENGTH()/2 < 0 || position + searchParams.getMAX_AMPLICON_LENGTH()/2 >= sequence.length()) continue; // terminate search with empty result (last element!)
				String wholeAmplicon = sequence.substring(position - searchParams.getMAX_AMPLICON_LENGTH()/2, position + searchParams.getMAX_AMPLICON_LENGTH()/2);
				
				if(doesNotContainAnotherRSS(wholeAmplicon)){
					// compute start and end for scan region for 5' primer search
					int forwardStart = position - searchParams.getMAX_AMPLICON_LENGTH()/2;
					int forwardEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
					
					int reverseStart = position + searchParams.getMAX_AMPLICON_LENGTH()/2;
					int reverseEnd = position + (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
					
					// compute start and end for 5' scan region for hybridization probe search
					int hybProbeStart = position - 1;
					int hybProbeEnd = position - (searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH());
					
					if(forwardStart < 0 || forwardEnd >= sequence.length() || reverseStart < 0 || reverseEnd >= sequence.length() || hybProbeStart < 0 || hybProbeStart >= sequence.length()) continue; //break; // terminate search with empty result (last element!)
					if(searchParams.getMIN_AMPLICON_LENGTH()/2 - searchParams.getMAX_PRIMER_LENGTH() < searchParams.getTAQMAN_MAX_PRIMER_LENGTH()) throw new IllegalArgumentException("MIN_AMPLICON_LENGTH/2 - MAX_PRIMER_LENGTH must be > TAQMAN_MAX_PRIMER_LENGTH! Scan regions for forward primer and HybProbe must NOT overlap!");
					
					assert(forwardEnd+1 > forwardStart);
					assert(reverseStart+1 > reverseEnd);
					assert(hybProbeStart+1 > hybProbeEnd);
					
					if(forwardStart < 0 || forwardEnd+1<=forwardStart || reverseEnd >= sequence.length()) continue;
					
					String forwardScanRegion = sequence.substring(forwardStart, forwardEnd+1);
					String reverseScanRegion = SeqTools.revcompDNA(sequence.substring(reverseEnd, reverseStart+1).toCharArray());
					String hybProbeScanRegion = SeqTools.revcompDNA(sequence.substring(hybProbeEnd, hybProbeStart+1).toCharArray());
					
					validUpstreamPrimers = this.enumeratePrimers(forwardScanRegion, PrimerTypes.forwardPrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), forwardStart, returnOrderedPrimerLists, position - forwardStart, nextBestSite, enzyme, searchParams);
					if(validUpstreamPrimers.length != 0) validDownstreamPrimers = this.enumeratePrimers(reverseScanRegion, PrimerTypes.reversePrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), reverseStart, returnOrderedPrimerLists, -position + reverseStart, nextBestSite, enzyme, searchParams);
					else{
						validDownstreamPrimers = new Primer[0];
						if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Empty upstream primer list");
					}
					if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && searchParams.isPickTaqManProbe()) validHybridizationProbes = this.enumeratePrimers(hybProbeScanRegion, PrimerTypes.hybridizationProbe, searchParams.getTAQMAN_MIN_PRIMER_LENGTH(), searchParams.getTAQMAN_MAX_PRIMER_LENGTH(), hybProbeStart, returnOrderedPrimerLists, 1, nextBestSite, enzyme, searchParams);
					else{
						validHybridizationProbes = new Primer[0];
						if(searchParams.isPRINT_DEBUG_LOG() && validUpstreamPrimers.length == 0) System.out.println("Empty upstream primer list!"); 
						else if(searchParams.isPRINT_DEBUG_LOG() && validDownstreamPrimers.length == 0) System.err.println("Empty downstream primer list");
					}
					
					if(searchParams.isPRINT_DEBUG_LOG() && searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0) System.err.println("Empty TaqMan probe list!");
					
					// remember new, optimal restriction site 'nextBestSite' and all of its valid primers instead of 'old', unsuitable restriction site 'site'
					if(validUpstreamPrimers.length != 0 && validDownstreamPrimers.length != 0 && (!searchParams.isPickTaqManProbe() || validHybridizationProbes.length != 0)){
						nextBestSite.setValidUpstreamPrimers(validUpstreamPrimers);
						nextBestSite.setValidDownstreamPrimers(validDownstreamPrimers);
						if(searchParams.isPickTaqManProbe()) nextBestSite.setValidTaqManProbes(validHybridizationProbes);
						nextBestSite.wasScannedForPrimers(true);
						// replace 'site' by 'nextBestSite' (in result vector) and continue with next 'optimal restriction site'
						sites[refinePosition] = nextBestSite;
						break;
					}
				}else{
					if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Found another RSS in this amplicon!");
				}
			}
			//else throw new DuplicateUseOfRestrictionSiteException("Restriction site at position " + nextBestSite.getPosition() + " is used multiple times!");
		}
		if(nextBestSite != null && (validUpstreamPrimers.length == 0 || validDownstreamPrimers.length == 0 || (searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0)) && !allRSSIterator.hasNext()) throw new EmptyResultSetException("There are no more restriction sites to scan for in sequence region " + refinePosition + ": " + nextBestSite.getSequenceRegion().getSeqRegionStart() + "-" + nextBestSite.getSequenceRegion().getSeqRegionEnd() + " (currRSS: " + (nextBestSite.getSequenceRegion().getRestricionSiteIndex(nextBestSite)+1) + "/" + nextBestSite.getSequenceRegion().getNumRestrictionSites() + ")");
		else if(nextBestSite == null && (validUpstreamPrimers.length == 0 || validDownstreamPrimers.length == 0 || (searchParams.isPickTaqManProbe() && validHybridizationProbes.length == 0)) && !allRSSIterator.hasNext()) throw new EmptyResultSetException("There are no more restriction sites to scan for in sequence region " + refinePosition);
		
		return sites;
	}
	
	public int getPrimerCount(){
		return this.primerCount;
	}
	
	public int getPrimerIndexCount(){
		return this.primerIndexCount;
	}
	
	public void setDoStat(boolean stat){
		this.doStat = stat;
	}
	
	public static void main(String[] args){
		if(args.length < 7){
			System.out.println("Usage: PrimerSearch <Target Region Sequence> <Background Region Index> <# Primer Pairs> <Enzyme Name> <Enzyme Site> <Enzyme Forward Cut Position> <Enzyme Reverse Cut Position>");
			System.out.println();
			System.exit(1);
		}
		PrimerSearch search = new PrimerSearch();
		PrimerSearchParameters searchParameters = null;
		
		String scanRegionFilename = args[0];
		String backgroundFilename = args[1];
		int nbPrimers = Integer.parseInt(args[2]);
		String enzymeName = args[3];
		String enzymeSite = args[4];
		search.setEnzymePatterns(Pattern.compile(".*" + enzymeSite.toUpperCase() + ".*"), Pattern.compile(".*" + SeqTools.revcompDNA(enzymeSite.toUpperCase().toCharArray()) + ".*"));
		int enzymeFwCutPos = Integer.parseInt(args[5]);
		int enzymeRevCutPos = Integer.parseInt(args[6]);
		
		if(args.length > 10)  Constants.doSortedDOPTScreen = Boolean.parseBoolean(args[9]);
		boolean benchmark = (args.length > 10) ? Boolean.parseBoolean(args[10]) : false;
		
//		Constants.PRIMER_DELTA_TM_WEIGHT = Double.parseDouble(args[7]);
//		Constants.PRIMER_DELTA_GC_WEIGHT = Double.parseDouble(args[8]);
//		Constants.PRIMER_DELTA_LENGTH_WEIGHT = Double.parseDouble(args[9]);
//		Constants.PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT = Double.parseDouble(args[10]);
//		Constants.SELF_ALIGNMENT_WEIGHT = Double.parseDouble(args[11]);
//		Constants.SELF_END_ALIGNMENT_WEIGHT = Double.parseDouble(args[12]);
//		Constants.PAIR_ALIGNMENT_WEIGHT = Double.parseDouble(args[13]);
//		Constants.PAIR_END_ALIGNMENT_WEIGHT = Double.parseDouble(args[14]);
//		Constants.PRIMER_FALSE_POSITIVES_WEIGHT = Double.parseDouble(args[15]);
		
		SimpleTimer timer = new SimpleTimer();
		try{
			
			timer.startTimer();
			if(!benchmark) System.out.println("Reading sequences:");
			if(!benchmark) System.out.print("Reading scan region - ");
			if(benchmark) System.out.print("SeqFile: "  + scanRegionFilename + " Index: " + backgroundFilename);
			
			// reading scan region
//			BufferedReader reader = new BufferedReader(new FileReader(scanRegionFilename));
//			StringBuffer sequence = new StringBuffer();
//			String line;
//			Pattern fastaStart = Pattern.compile(">.*");
//			while((line = reader.readLine()) != null){
//				if(!fastaStart.matcher(line).matches()){
//					sequence.append(line.trim());
//				}
//			}
			// read first sequence from fasta
			SlimFastaParser parser = new SlimFastaParser(new File(scanRegionFilename));
			Vector<SimpleContigImpl> contigs = new Vector<SimpleContigImpl>();
			int counter = 0;
			while(parser.hasNextContig()){
				contigs.add(parser.parseNextContigIgnoreCase());
				counter++;
			}
			if(!benchmark) System.out.println("done " + timer.getTimeString() + " - read " + counter + " sequences");
			
			// reading background sequence
//			System.out.print("Reading background sequence - ");
//			reader = new BufferedReader(new FileReader(backgroundFilename));
//			StringBuffer backgroundSeq = new StringBuffer();
//			while((line = reader.readLine()) != null){
//				if(!fastaStart.matcher(line).matches()){
//					//background += line.trim();
//					backgroundSeq.append(line.trim());
//				}
//				else{
//					// prohibit direct concatenation of different sequences - pseudo-primers could be found
//					backgroundSeq.append(Constants.REPETITIVE_ELEMENT_CHARACTER);
//				}
//			}
//			System.out.println("done " + timer.getTimeString());
			
			//int maxWordLength = MyExtendedMath.ceil(Constants.MAX_PRIMER_LENGTH * Constants.SAFE_MISPRIMING_SEQUENCE_PRECENT_CUTOFF);
			//int maxWordLength = Constants.MAX_PRIMER_LENGTH ;
			//int maxWordLengthScan = Constants.MAX_PRIMER_LENGTH;
			//int maxWordLengthBackground = 10;
			//System.out.println("Constructing suffix trees of sequences - max word length: " + 10);
			if(!benchmark) System.out.print("Constructing scanRegion index - ");
			//DNASequenceIndex scanRegionTree = new DNASuffixTrie();
//			DNASequenceIndex scanRegionIndex = new EnhancedSuffixArrayFatOpt("");
//			//DNASequenceIndex scanRegionIndex = new EnhancedSuffixArrayIntOptRecomp("");
//			//DNASequenceIndex scanRegionIndex = new StringIndex();
//			scanRegionIndex.createIndex(sequence, sequence.length(), true);
////			scanRegionTree.createIndex(sequence.toString(), 10, true);
//			System.out.println("done " + timer.getTimeString());
			
			if(!benchmark) System.out.print("Constructing background index - ");
			//DNASequenceIndex backgroundTree = new DNASuffixTrie();
			//DNASequenceIndex backgroundIndex = new StringIndex();
			//DNASequenceIndex backgroundIndex = new DummyIndex();
//			primerDesign.dsc.indexStructures.DNASequenceIndex backgroundIndex;
//			//DNASequenceIndex backgroundIndex = new EnhancedSuffixArrayIntOptRecomp("");
//			if(scanRegionFilename.equals(backgroundFilename)){
//				backgroundIndex = DNASequenceIndexDeserializer.deserialize(backgroundFilename);
//			}
//			else{
//				System.out.print("deserializing background sequence index - ");
//				backgroundIndex = DNASequenceIndexDeserializer.deserialize(backgroundFilename);
//				//backgroundIndex.createIndex("", 0, true);
//				//backgroundIndex.createIndex(SeqTools.readFastaFile(args[1]), Integer.MAX_VALUE, true);
//				//backgroundIndex = EnhancedSuffixArrayIntOptRecomp.deserialize(backgroundFilename);
//				//backgroundTree.createIndex(backgroundSeq.toString(), backgroundSeq.length(), true);
//				//backgroundTree.createIndex(backgroundSeq.toString(), backgroundSeq.length(), true);
//			}
			
			RestrictionEnzyme enzyme = new RestrictionEnzyme(enzymeName, DNATools.createDNA(enzymeSite), enzymeFwCutPos, enzymeRevCutPos	);
			
			searchParameters = new PrimerSearchParameters(search);
			searchParameters.setNumPrimers(nbPrimers);
			searchParameters.setEnzyme(enzyme);
			if(args.length > 7){
				if(args[7].equals("Ppacificus")) searchParameters.setTargetOrganism(TargetOrganisms.Ppacificus);
				else if(args[7].equals("Celegans")) searchParameters.setTargetOrganism(TargetOrganisms.Celegans);
				else if(args[7].equals("Scerevisiae")) searchParameters.setTargetOrganism(TargetOrganisms.Scerevisiae);
				else if(args[7].equals("Dmelanogaster")) searchParameters.setTargetOrganism(TargetOrganisms.Dmelanogaster);
				else if(args[7].equals("Mmusculus")) searchParameters.setTargetOrganism(TargetOrganisms.Mmusculus);
				else if(args[7].equals("Hsapiens")) searchParameters.setTargetOrganism(TargetOrganisms.Hsapiens);
				else throw new IllegalArgumentException("Unknown organism specified!");
			}
			else searchParameters.setTargetOrganism(TargetOrganisms.Ppacificus);
			//searchParameters.setPRINT_DEBUG_LOG(true);
			searchParameters.setContigs(contigs.toArray(new SimpleContigImpl[contigs.size()]));
			//searchParameters.setScanRegionIndex(scanRegionIndex);
			//searchParameters.setBackgroundIndex(backgroundIndex);
			searchParameters.setConfigFile(PrimerDesignWebProperites.parseConfigFile(new FileInputStream("/Users/froehler/Desktop/Projects/EclipseWorkspace/PrimerDesign/deploy/3CPrimerDesign-ServletConfig-localMAC.xml")));
			
			PrimerMisprimingCheck misprimingCheck = PrimerMisprimingCheckDeserializer.deserialize(backgroundFilename, searchParameters);
			
			searchParameters.setPrimerMisprimingCheck(misprimingCheck);
			if(args.length > 8) searchParameters.setPickTaqManProbe(Boolean.parseBoolean(args[8]));
			
			if(!benchmark) System.out.println("done " + timer.getTimeString());
			else System.out.print("\t preprocessing(seq+idx):\t" + timer.getTimeString());
//			backgroundIndex = new DummyIndex();
//			backgroundIndex.createIndex("", 0, true);
//			backgroundIndex = new StringIndex();
//			backgroundIndex.createIndex(sequence.toString(), sequence.length(), true);
			
			if(benchmark) System.out.print("\tpairs: " + nbPrimers + "\t" + "TaqMan: " + searchParameters.isPickTaqManProbe() + "\t" + " dOptsorted: " + Constants.doSortedDOPTScreen + "\t");
			
			if(!benchmark) System.out.println("Scanning sequence for valid primers");			
			search.doStat = searchParameters.isComputeScanningStatistics();
			
			RestrictionSite[] optimalSites = search.naivePrimerSearch(searchParameters);
			if(!benchmark) System.out.println("Scanning sequence for valid primers took " + timer.getTimeString());
			else System.out.print("\t scan:\t" + timer.getTimeString());
						
//			System.out.println("The following restriction sites and primers were computed in the sequence");
//			//RestrictionSite rss;
//			for(int i=0; i< optimalSites.length;i++){
//				System.out.println("Valid upstream primers");
//				for(int j=0; j<optimalSites[i].getValidUpstreamPrimers().length; j++) System.out.println(optimalSites[i].getUpstreamPrimer(j).toString());
//				System.out.println("Valid downstream primers");
//				for(int j=0; j<optimalSites[i].getValidDownstreamPrimers().length; j++) System.out.println(optimalSites[i].getDownstreamPrimer(j).toString());
//				if(searchParameters.isPickTaqManProbe()){
//					System.out.println("Valid hybridization probes primers");
//					for(int j=0; j<optimalSites[i].getValidTaqManProbes().length; j++) System.out.println(optimalSites[i].getTaqManProbe(j).toString());
//				}
////				rss = optimalSites[i];
////				System.out.println("\tPosition: " + rss.getPosition());
////				
////				Primer[] primers = rss.getValidUpstreamPrimers();
////				for(int j=0; j<primers.length; j++){
////					Primer primer = primers[j];
////					System.out.println("Primer: " + primer.getSequence() + " " + primer.getDistanceToRSS() + " " + primer.getGcContent() + " " + primer.getMeltingTemp() + " " + primer.getDistanceToOptimalPrimer() + " " + primer.getSelfAlignmentScore() + " " + primer.getSelfEndAlignmentScore());
////
////				}
//			}
//			System.out.println("Displaying valid primers took " + timer.getTimeString());
			
			if(!benchmark) System.out.println("Computing best primer set");

			
			SimpleGreedyPrimerPairPicking picker = new SimpleGreedyPrimerPairPicking();
			PrimerPairSet bestPrimerPairSet = null; // = picker.pickBestPrimerSet(optimalSites, searchParameters);
			
//			AdvancedGreedyPrimerPairPicking picker = new AdvancedGreedyPrimerPairPicking();
//			PrimerPairSet bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
			
			try{
				bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
			}
			catch(EmptyResultSetException e){
				if(!benchmark) System.out.println("\n" + e.getMessage());
				if(!benchmark) System.out.println("-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!");
				if(!benchmark) if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
				if(!benchmark) if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
				if(!benchmark) if(searchParameters.isComputePickingStatistics() && searchParameters.getPickingStat() != null) System.out.println(searchParameters.getPickingStat().printStat());
				if(benchmark) System.out.println("\t FAILED!\t time:\t" + timer.getTotalTimestring());
				System.exit(1);
			}			
			
//			// while alignment scores of best pairset are unacceptable
//			while(bestPrimerPairSet.getMaxPairAlignScore() > searchParameters.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || bestPrimerPairSet.getMaxPairAlignEndScore() > searchParameters.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
//				try{
//					PrimerPair currentPair;
//					PrimerAlignmentScores scores;
//					FastVector temp = new FastVector();
//					// screen each pair in pairset for acceptable parameters
//					for(int i=0; i<bestPrimerPairSet.getPrimerPairs().size(); i++){
//						currentPair = bestPrimerPairSet.getPrimerPair(i);
//						scores = currentPair.getMaxPairElementsAlignmentScore();
//						// eval current pair, refine iff unacceptable
//						if(scores.getPairScore() > searchParameters.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() || scores.getPairEndScore() > searchParameters.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE()){
//							temp = new FastVector();
//							for(int j=0; j<optimalSites.length; j++) temp.addElement(optimalSites[j]);
//							temp = search.refineSeqRegion(currentPair.getForwardPrimer().getRestrictionSite(), temp, searchParameters, sequence.toString(), scanRegionIndex, backgroundIndex, enzyme, false);
//							assert(temp.size() == optimalSites.length);
//							for(int j=0; j<temp.size(); j++) optimalSites[j] = (RestrictionSite) temp.elementAt(j);
//						}
//					}
//					
//					bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
//				}
//				catch(Exception e){
//					break;
//				}
//				
//				//System.out.println("NO valid primer pairs can be picked - alignment values exceed your parameters!");
//			}
			if(!benchmark) System.out.println("Computing best primer set took " + timer.getTimeString());
			else System.out.print("\t picking:\t" + timer.getTimeString());
			
//			} 
			if(!benchmark){
				if(bestPrimerPairSet == null){
					System.out.println("-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!");
					if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
					if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
					if(searchParameters.isComputePickingStatistics()) System.out.println(searchParameters.getPickingStat().printStat());
				}
				else{
					System.out.println("Displaying best primer set");
					System.out.println(String.format("%10s\t%s", "Type", Primer.toFormattedStringDescription()));
	//				for(int i=0; i<bestPrimerSet.getNumberOfForwardPrimers(); i++){
	//					Primer currentPrimer = bestPrimerSet.getForwardPrimer(i);
	//					//System.out.println(currentPrimer.getRelativePosition() + "\t" + currentPrimer.getSequence() + "\t" + currentPrimer.getDistanceToRSS() + "\t" + currentPrimer.getGcContent() + "\t" + currentPrimer.getMeltingTemp() + "\t" + currentPrimer.getDistanceToOptimalPrimer() + "\t" + currentPrimer.getSelfAlignmentScore() + "\t" + currentPrimer.getSelfEndAlignmentScore());
	//					System.out.println(currentPrimer.toString());
					System.out.println(bestPrimerPairSet.toFormattedString());
					System.out.println("worst primer-  primer PA/PEA score (incl probes): " + bestPrimerPairSet.getMaxPairAlignScore() + " " + bestPrimerPairSet.getMaxPairAlignEndScore());
					System.out.println("Displaying best primer set took " + timer.getTimeString());
					
					System.out.println("### HashCode of pair set: " + bestPrimerPairSet.getHashCode() + " ###");
	
					if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
					if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
					if(searchParameters.isComputePickingStatistics()) System.out.println(searchParameters.getPickingStat().printStat());
				}
				
				System.out.println("Total runtime: " + timer.getTotalTimestring());
				System.out.println("Evaluated " + PrimerSearch.primerCount + " candidate primers in total.");
				System.out.println("Evaluated " + PrimerSearch.primerIndexCount + " candidate primers in index in total.");
			}
			else if(benchmark && bestPrimerPairSet == null){
				System.out.print("\t FAILED!\t time:\t" + timer.getTotalTimestring());
			}else if(benchmark){
				System.out.print("\t SUCCESSFUL!\t time:\t" + timer.getTotalTimestring());
			}
			// write timings to logfile
//			BufferedWriter writer = new BufferedWriter(new FileWriter("timings.log", true)); 
//			writer.write(enzymeName + "\t" + nbPrimers + "\t" + timer.getTotalTimestring() + "\n");
//			writer.close();
		}
		catch(Exception e){
			if(!benchmark){
				e.printStackTrace();
				if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
			}else{
				System.out.print("\t\t\t\t\t FAILED!\t time:\t" + timer.getTotalTimestring());
			}
		}
		System.out.println();
	}
}
