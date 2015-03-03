/**
 * A simple parameter printer class.
 */
package primerDesign.web.tools;

import javax.servlet.http.HttpServletRequest;

import primerDesign.util.PrimerSearchParameters;

/**
 * @author Sebastian Fršhler
 *
 */
public class PrintParameters {
	public String print3CServletParameters(HttpServletRequest request){
		StringBuffer buffy = new StringBuffer();
		
		buffy.append("The following parameters were used during your search:\n");
		buffy.append("Number of primer pairs: " + request.getParameter("numberOfPrimerPairs") + " primer pairs\n");
		buffy.append("Target Organism: " + request.getParameter("TargetOrganism") + "\n");
		buffy.append("Restriction enzyme: " + request.getParameter("restrictionEnzyme").split(":")[0] + "\n\n");
		
		buffy.append("Primer Properties:\n");
		buffy.append("Min primer length: " + request.getParameter("minPrimerLength") + " Opt primer length: " + request.getParameter("optPrimerLength") + " Max primer length: " + request.getParameter("maxPrimerLength") + "\n");
		buffy.append("Min primer TM: " + request.getParameter("minPrimerTM") + " Opt primer TM: " + request.getParameter("optPrimerTM") + " Max primer TM: " + request.getParameter("maxPrimerTM") + "\n");
		buffy.append("Min primer GC: " + request.getParameter("minPrimerGC") + " Opt primer GC: " + request.getParameter("optPrimerGC") + " Max primer GC: " + request.getParameter("maxPrimerGC") + "\n");
		buffy.append("Max primer TM difference: " + request.getParameter("maxPrimerTMDifference") + "\n");
		buffy.append("Max primer self alignment: " + request.getParameter("saMax") + " Max primer self end alignment: " + request.getParameter("seaMax") + "\n");
		buffy.append("Max primer pair alignment: " + request.getParameter("paMax") + " Max primer pair end alignment: " + request.getParameter("peaMax") + "\n");
		
		buffy.append("\nProbe Properties:\n");
		buffy.append("Min probe length: " + request.getParameter("minProbeLength") + " Opt probeLength: " + request.getParameter("optProbeLength") + " Max probe length: " + request.getParameter("maxProbeLength") + "\n");
		buffy.append("Min probe TM: " + request.getParameter("minProbeTM") + " Opt probe TM: " + request.getParameter("optProbeTM") + " Max probe TM: " + request.getParameter("maxProbeTM") + "\n");
		buffy.append("Min probe GC: " + request.getParameter("minProbeGC") + " Opt probe GC: " + request.getParameter("optProbeGC") + " Max probe GC: " + request.getParameter("maxProbeGC") + "\n");
		buffy.append("Max probe TM difference: " + request.getParameter("minProbeTMDifference") + "\n");
		buffy.append("Max probe self alignment: " + request.getParameter("probeSaMax") + " Max probe self end alignment: " + request.getParameter("probeSeaMax") + "\n");
		buffy.append("Max probe pair alignment: " + request.getParameter("probePaMax") + " Max probe pair end alignment: " + request.getParameter("probePeaMax") + "\n");
		
		buffy.append("\nAmplicon Properties:\n");
		buffy.append("Min amplicon length: "  + request.getParameter("minAmpliconLength") + " Opt amplicon length: " + request.getParameter("optAmpliconLength") + " Max amplicon length: " + request.getParameter("maxAmpliconLength") + "\n");
		buffy.append("Safe FP amplicon length: " + request.getParameter("safeFPAmpliconLength") + "\n");
		
		buffy.append("\nWeights:\n");
		buffy.append("delta TM weight: " + request.getParameter("deltaTMWeight") + "\n");
		buffy.append("delta GC weight: " + request.getParameter("deltaGCWeight") + "\n");
		buffy.append("delta Length weight: " + request.getParameter("deltaLengthWeight") + "\n");
		buffy.append("delta Distance to RSS weight: " + request.getParameter("deltaDRSSWeight") + "\n");
		buffy.append("Self alignment weight: " + request.getParameter("deltaSAWeight") + "\n");
		buffy.append("Self end alignment weight: " + request.getParameter("deltaSEAWeight") + "\n");
		buffy.append("Pair alignment weight: " + request.getParameter("deltaPAWeight") + "\n");
		buffy.append("Pair end alignment weight: " + request.getParameter("deltaPEAWeight") + "\n");
		buffy.append("False positive weight: " + request.getParameter("FPWeight") + "\n");
		
		return buffy.toString();
	}
	
	public String print3CSearchparameters(PrimerSearchParameters params){
		StringBuffer buffy = new StringBuffer();
		
		buffy.append("The following parameters were used during your search:\n");
		buffy.append("------------------------------------------------------\n\n");
		buffy.append("Target organism: " + params.getTargetOrganism() + "\n");
		buffy.append("Number of primer pairs: " + params.getNumPrimers() + "\n");
		buffy.append("Restriction enzyme: " + params.getEnzymeName() + "\n\n");
		
		buffy.append("Primer Properties:\n");
		buffy.append("Min primer length: " + params.getMIN_PRIMER_LENGTH() + " Opt primer length: " + params.getOPT_PRIMER_LENGTH() + " Max primer length: " + params.getMAX_PRIMER_LENGTH() + "\n");
		buffy.append("Min primer TM: " + params.getMIN_TM() + " Opt primer TM: " + params.getOPT_TM() + " Max primer TM: " + params.getMAX_TM() + "\n");
		buffy.append("Min primer GC: " + params.getMIN_GC() + " Opt primer GC: " + params.getOPT_GC() + " Max primer GC: " + params.getMAX_GC() + "\n");
		buffy.append("Max primer TM difference: " + params.getMAX_PRIMER_TM_DIFFERENCE() + "\n");
		buffy.append("Max primer self alignment: " + params.getMAX_PRIMER_SELF_ALIGNMENT_SCORE() + " Max primer self end alignment: " + params.getMAX_PRIMER_SELF_END_ALIGNMENT_SCORE() + "\n");
		buffy.append("Max primer pair alignment: " + params.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() + " Max primer pair end alignment: " + params.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() + "\n");
		
		buffy.append("\nProbe Properties:\n");
		buffy.append("Min probe length: " + params.getTAQMAN_MIN_PRIMER_LENGTH() + " Opt probeLength: " + params.getTAQMAN_OPT_PRIMER_LENGTH() + " Max probe length: " + params.getTAQMAN_MAX_PRIMER_LENGTH() + "\n");
		buffy.append("Min probe TM: " + params.getTAQMAN_MIN_TM() + " Opt probe TM: " + params.getTAQMAN_OPT_TM() + " Max probe TM: " + params.getTAQMAN_MAX_TM() + "\n");
		buffy.append("Min probe GC: " + params.getTAQMAN_MIN_GC() + " Opt probe GC: " + params.getTAQMAN_OPT_GC() + " Max probe GC: " + params.getTAQMAN_MAX_GC() + "\n");
		buffy.append("Max probe TM difference: " + params.getMIN_TAQMAN_TM_DIFFERENCE() + "\n");
		buffy.append("Max probe self alignment: " + params.getMAX_TAQMAN_SELF_ALIGNMENT_SCORE() + " Max probe self end alignment: " + params.getMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE() + "\n");
		buffy.append("Max probe pair alignment: " + params.getMAX_TAQMAN_PAIR_ALIGNMENT_SCORE() + " Max probe pair end alignment: " + params.getMAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE() + "\n");
		
		buffy.append("\nAmplicon Properties:\n");
		buffy.append("Min amplicon length: "  + params.getMIN_AMPLICON_LENGTH() + " Opt amplicon length: " + params.getOPT_AMPLICON_LENGTH() + " Max amplicon length: " + params.getMAX_AMPLICON_LENGTH() + "\n");
		buffy.append("Safe FP amplicon length: " + params.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH() + "\n");
		
		buffy.append("\nWeights:\n");
		buffy.append("delta TM weight: " + params.getPRIMER_DELTA_TM_WEIGHT() + "\n");
		buffy.append("delta GC weight: " + params.getPRIMER_DELTA_GC_WEIGHT() + "\n");
		buffy.append("delta Length weight: " + params.getPRIMER_DELTA_LENGTH_WEIGHT() + "\n");
		buffy.append("delta Distance to RSS weight: " + params.getPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT() + "\n");
		buffy.append("Self alignment weight: " + params.getSELF_ALIGNMENT_WEIGHT() + "\n");
		buffy.append("Self end alignment weight: " + params.getSELF_END_ALIGNMENT_WEIGHT() + "\n");
		buffy.append("Pair alignment weight: " + params.getPAIR_ALIGNMENT_WEIGHT() + "\n");
		buffy.append("Pair end alignment weight: " + params.getPAIR_END_ALIGNMENT_WEIGHT() + "\n");
		buffy.append("False positive weight: " + params.getPRIMER_FALSE_POSITIVES_WEIGHT() + "\n");
		
		return buffy.toString();
	}
}
