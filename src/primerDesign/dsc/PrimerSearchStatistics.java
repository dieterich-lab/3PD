/**
 * 
 */
package primerDesign.dsc;

import java.text.NumberFormat;

import primerDesign.util.PrimerSearchParameters;

/**
 * Wraps up some primer search statistics.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerSearchStatistics {
	private int primerCount = 0;
	private int primerIndexCount = 0;
	private int anotherRSSinAmpliconCount = 0;
	private int enumeratePrimersCount = 0;
	private int emptyUpstreamPrimerList = 0;
	private int nonEmptyUpstreamPrimerList = 0;
	private int emptyDownstreamPrimersList = 0;
	private int nonEmptyDownstreamPrimersList;
	private int emptyTaqManPrimersList = 0;
	private int nonEmptyTaqManPrimersList = 0;
	private int emptyUpDownTaqManList = 0;
	private int nonEmptyUpDownTaqManList = 0;
	private int primerAmbiguityCodeCount = 0;
	private int singleNucleotideRepeatPatternCount = 0;
	private int gcClampCount = 0;
	private int repetitiveElementPatternCount = 0;
	private int restrictionEnzymePatternCount = 0;
	private int invalidPrimerCount = 0;
	private int validPrimerCount = 0;
	private int noMisprimingsCount = 0;
	private int misprimingCount = 0;
	private int invalidTaqManCount = 0;
	private int validTaqManCount = 0;
	private int invalidPrimerTM = 0;
	private int invalidPrimerTMabove = 0;
	private int invalidPrimerTMbelow = 0;
	private int invalidPrimerGC = 0;
	private int invalidPrimerGCabove = 0;
	private int invalidPrimerGCbelow = 0;
	private int invalidPrimerSA = 0;
	private int invalidPrimerSEA = 0;
	private int invalidProbeTM = 0;
	private int invalidProbeTMabove = 0;
	private int invalidProbeTMbelow = 0;
	private int invalidProbeGC = 0;
	private int invalidProbeGCabove = 0;
	private int invalidProbeGCbelow = 0;
	private int invalidProbeSA = 0;
	private int invalidProbeSEA = 0;
	
	public void incPrimerCount(){
		this.primerCount++;
	}
	
	public void incPrimerIndexCount(){
		this.primerIndexCount++;
	}
	
	public void incAnotherRSSinAmpliconCount(){
		this.anotherRSSinAmpliconCount++;
	}
	
	public void incEnumeratePrimersCount(){
		this.enumeratePrimersCount++;
	}
	
	public void incEmptyUpstreamPrimersList(){
		this.emptyUpstreamPrimerList++;
	}
	
	public void incNonemptyUpstreamPrimersList(){
		this.nonEmptyUpstreamPrimerList++;
	}
	
	public void incEmptyDownstreamPrimersList(){
		this.emptyDownstreamPrimersList++;
	}
	
	public void incNonEmptyDownstreamPrimersList(){
		this.nonEmptyDownstreamPrimersList++;
	}
	
	public void incEmptyTaqManPrimersList(){
		this.emptyTaqManPrimersList++;
	}
	
	public void incNonEmptyTaqManPrimersList(){
		this.nonEmptyTaqManPrimersList++;
	}
	
	public void incEmptyUpDownTaqManList(){
		this.emptyUpDownTaqManList++;
	}
	
	public void incNonEmptyUpDpwnTaqManList(){
		this.nonEmptyUpDownTaqManList++;
	}
	
	public void incPrimerAmbiguityCode(){
		this.primerAmbiguityCodeCount++;
	}
	
	public void incSingleNucleotideRepeatPattern(){
		this.singleNucleotideRepeatPatternCount++;
	}
	
	public void incGCClampCount(){
		this.gcClampCount++;
	}
	
	public void incRepetitiveElementPattern(){
		this.repetitiveElementPatternCount++;
	}
	
	public void incRestrictionEnzymePattern(){
		this.restrictionEnzymePatternCount++;
	}
	
	public void incInvalidPrimer(){
		this.invalidPrimerCount++;
	}
	
	public void incValidPrimer(){
		this.validPrimerCount++;
	}
	
	public void incMispriming(){
		this.misprimingCount++;
	}
	
	public void incNoMispriming(){
		this.noMisprimingsCount++;
	}
	
	public void incInvalidTaqManPrimer(){
		this.invalidTaqManCount++;
	}
	
	public void incValidTaqManPrimer(){
		this.validTaqManCount++;
	}
	
	public void evaluatePrimer(Primer primer, PrimerSearchParameters params){
		if(!primer.getPrimerType().equals(PrimerTypes.forwardPrimer) && !primer.getPrimerType().equals(PrimerTypes.reversePrimer)) throw new IllegalArgumentException("This method is designed for evaluating forward and reverse primers!");
		if(primer.getMeltingTemp() < params.getMIN_TM()){
			this.invalidPrimerTM++;
			this.invalidPrimerTMbelow++;
		}
		if(primer.getMeltingTemp() > params.getMAX_TM()){
			this.invalidPrimerTM++;
			this.invalidPrimerTMabove++;
		}
		if(primer.getGcContent() < params.getMIN_GC()){
			this.invalidPrimerGC++;
			this.invalidPrimerGCbelow++;
		}
		if(primer.getGcContent() > params.getMAX_GC()){
			this.invalidPrimerGC++;
			this.invalidPrimerGCabove++;
		}
		if(primer.getSelfAlignmentScore() > params.getMAX_PRIMER_SELF_ALIGNMENT_SCORE()) this.invalidPrimerSA++;
		if(primer.getSelfEndAlignmentScore() > params.getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE()) this.invalidPrimerSEA++;
	}
	
	public void evaluateProbe(Primer primer, PrimerSearchParameters params){
		if(!primer.getPrimerType().equals(PrimerTypes.hybridizationProbe)) throw new IllegalArgumentException("This method is designed for evaluating TaqMan probe primers!");
		if(primer.getMeltingTemp() < params.getMIN_TM()){
			this.invalidProbeTM++;
			this.invalidProbeTMbelow++;
		}
		if(primer.getMeltingTemp() > params.getMAX_TM()){
			this.invalidProbeTM++;
			this.invalidProbeTMabove++;
		}
		if(primer.getGcContent() < params.getMIN_GC()){
			this.invalidProbeGC++;
			this.invalidProbeGCbelow++;
		}
		if(primer.getGcContent() > params.getMAX_GC()){
			this.invalidProbeGC++;
			this.invalidProbeGCabove++;
		}
		if(primer.getSelfAlignmentScore() > params.getMAX_PRIMER_SELF_ALIGNMENT_SCORE()) this.invalidProbeSA++;
		if(primer.getSelfEndAlignmentScore() > params.getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE()) this.invalidProbeSEA++;
	}
	
	public String getStat(){
		StringBuffer buffy = new StringBuffer();
		buffy.append("\n");
		buffy.append("Primer search statistics:\n");
		buffy.append("-------------------------\n");
		buffy.append("Enumerated primers: " + this.primerCount + "\n");
		buffy.append("Enumerated primers checked in index: " + this.primerIndexCount + "\n");
		buffy.append("Calls to 'enumeratePrimers: " + this.enumeratePrimersCount + "\n");
		buffy.append("Empty upstream primers list: " + this.emptyUpstreamPrimerList + "\n");
		buffy.append("Non-Empty upstream primers list: " + this.nonEmptyUpstreamPrimerList + "\n");
		buffy.append("Empty downstream primers list: " +  this.emptyDownstreamPrimersList + "\n");
		buffy.append("Non-Empty downstream primers list: " + this.nonEmptyDownstreamPrimersList + "\n");
		buffy.append("Empty TaqMan probes list: " + this.emptyTaqManPrimersList + "\n");
		buffy.append("Non-empty TaqMan probes list: " + this.nonEmptyTaqManPrimersList + "\n");
		buffy.append("Empty up+down+taqman list: " + this.emptyUpDownTaqManList + "\n");
		buffy.append("Non-Empty up+down+taqman list: " + this.nonEmptyUpDownTaqManList + "\n");
		buffy.append("Primers matching ambiguity pattern: " + this.primerAmbiguityCodeCount + "\n");
		buffy.append("Single nucleotide repeat count: " + this.singleNucleotideRepeatPatternCount + "\n");
		buffy.append("Repetitive element pattern count: " + this.repetitiveElementPatternCount + "\n");
		buffy.append("Restriction enzyme pattern count: " + this.restrictionEnzymePatternCount + "\n");
		buffy.append("GC clamp count: " + this.gcClampCount + "\n");
		buffy.append("Another RSS in amplicon: " + this.anotherRSSinAmpliconCount + "\n");
		buffy.append("Invalid primers: " + this.invalidPrimerCount + "\n");
		buffy.append("\tInvalid TM: " + this.invalidPrimerTM + " (below: " + this.invalidPrimerTMbelow + "/ above: " + this.invalidPrimerTMabove + ")\n");
		buffy.append("\tInvalid GC: " + this.invalidPrimerGC + " (below: " + this.invalidPrimerGCbelow + "/ above: " + this.invalidPrimerGCabove + ")\n");
		buffy.append("\tInvalid SA: " + this.invalidPrimerSA + "\n");
		buffy.append("\tInvalid SEA: " + this.invalidPrimerSEA + "\n");
		buffy.append("Valid primers: " + this.validPrimerCount + "\n");
		buffy.append("Invalid probes: " + this.invalidTaqManCount + "\n");
		buffy.append("\tInvalid TM: " + this.invalidProbeTM + " (below: " + this.invalidProbeTMbelow + "/ above: " + this.invalidProbeTMabove + ")\n");
		buffy.append("\tInvalid GC: " + this.invalidProbeGC + " (below: " + this.invalidProbeGCbelow + "/ above: " + this.invalidProbeGCabove + ")\n");
		buffy.append("\tInvalid SA: " + this.invalidProbeSA + "\n");
		buffy.append("\tInvalid SEA: " + this.invalidProbeSEA + "\n");
		buffy.append("Valid probes: " + this.validTaqManCount + "\n");
		buffy.append("Misprimings count: " + this.misprimingCount + "\n");
		buffy.append("No mispriming count: " + this.noMisprimingsCount + "\n");
		
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
		
		buffy.append("Primer search statistics:\n");
		buffy.append("-------------------------\n");
		buffy.append("Valid primers: " + format.format(this.validPrimerCount) + "\n");
		buffy.append("Invalid primers: " + format.format(this.invalidPrimerCount) + "\n");
		buffy.append("\tInvalid TM: " + format.format(this.invalidPrimerTM) + " (below: " + format.format(this.invalidPrimerTMbelow) + "/ above: " + format.format(this.invalidPrimerTMabove) + ")\n");
		buffy.append("\tInvalid GC: " + format.format(this.invalidPrimerGC) + " (below: " + format.format(this.invalidPrimerGCbelow) + "/ above: " + format.format(this.invalidPrimerGCabove) + ")\n");
		buffy.append("\tInvalid SA: " + format.format(this.invalidPrimerSA) + "\n");
		buffy.append("\tInvalid SEA: " + format.format(this.invalidPrimerSEA) + "\n");
		buffy.append("Valid probes: " + format.format(this.validTaqManCount) + "\n");
		buffy.append("Invalid probes: " + format.format(this.invalidTaqManCount) + "\n");
		buffy.append("\tInvalid TM: " + format.format(this.invalidProbeTM) + " (below: " + format.format(this.invalidProbeTMbelow) + "/ above: " + format.format(this.invalidProbeTMabove) + ")\n");
		buffy.append("\tInvalid GC: " + format.format(this.invalidProbeGC) + " (below: " + format.format(this.invalidProbeGCbelow) + "/ above: " + format.format(this.invalidProbeGCabove) + ")\n");
		buffy.append("\tInvalid SA: " + format.format(this.invalidProbeSA) + "\n");
		buffy.append("\tInvalid SEA: " + format.format(this.invalidProbeSEA) + "\n");
		buffy.append("Misprimings count: " + format.format(this.misprimingCount) + "\n");
		buffy.append("No mispriming count: " + format.format(this.noMisprimingsCount) + "\n");
		
		return buffy.toString();
	}
	
	public String toString(){
		return this.getStat();
	}
}
