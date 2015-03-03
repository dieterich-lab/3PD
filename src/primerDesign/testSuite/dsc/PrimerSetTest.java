package primerDesign.testSuite.dsc;

import junit.framework.TestCase;
import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.algo.SimpleAlignment;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.PrimerSet;
import primerDesign.dsc.PrimerTypes;
import primerDesign.util.PrimerSearchParameters;

public class PrimerSetTest extends TestCase {

	public void testAddBestPrimer() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		params.setPRIMER_CONCENTRATION(50e-9);
		params.setMONOVALENT_CATION_CONCENTRATION(50e-3);
		params.setDIVALENT_CATION_CONCENTRATION(0);
		params.setDNTP_CONCENTRATION(0);
		
		// add primer to one-element primerSet
		PrimerSet set = new PrimerSet(params);
		Primer primer = new Primer("ATGCATGCATGC", PrimerTypes.forwardPrimer, params);
		set.addForwardPrimer(primer);
		Primer[] primers = {new Primer("ATATATATATAT", PrimerTypes.forwardPrimer, params), new Primer("ATGCATGCATGC", PrimerTypes.forwardPrimer, params), new Primer("GCGCGCGCGCGC", PrimerTypes.forwardPrimer, params)};
		set.addBestPrimer(primers, PrimerTypes.forwardPrimer);
		assertEquals(primer, set.getForwardPrimers().elementAt(0));
		assertEquals(primers[0], set.getForwardPrimers().elementAt(1));
		
		// add primer to two-element primer set
		set = new PrimerSet(params);
		Primer primer1 = new Primer("ATATATATATAT", PrimerTypes.forwardPrimer, params);
		Primer primer2 = new Primer("GCGCGCGCGCGC", PrimerTypes.forwardPrimer, params);
		set.addForwardPrimer(primer1);
		set.addForwardPrimer(primer2);
		Primer[] candidates = {new Primer("ATATATATATAT", PrimerTypes.forwardPrimer, params), new Primer("ATGCATGCATGC", PrimerTypes.forwardPrimer, params), new Primer("GCGCGCGCGCGC", PrimerTypes.forwardPrimer, params)};
		set.addBestPrimer(candidates, PrimerTypes.forwardPrimer);
		assertEquals(primer1, set.getForwardPrimers().elementAt(0));
		assertEquals(primer2, set.getForwardPrimers().elementAt(1));
		assertEquals(candidates[1], set.getForwardPrimers().elementAt(2));
	}

	public void testScorePrimerTo() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		// score primer to primer set containing just itself
		Primer primer = new Primer("ATATATATATAT", PrimerTypes.forwardPrimer, params);
		PrimerSet set = new PrimerSet(params);
		set.addForwardPrimer(primer);
		PrimerAlignmentCalculation alignment = new SimpleAlignment(params.getA_t_basepair_score(), params.getG_c_basepair_score());
		PrimerAlignmentScores scores = alignment.computePairAlignment(primer.getSequence(), set.getForwardPrimer(0).getSequence());
		double expect = params.getPAIR_ALIGNMENT_WEIGHT() * Math.max(set.getMaxPairAlignmentScore(), scores.getPairScore()) + params.getPAIR_END_ALIGNMENT_WEIGHT() * Math.max(set.getMaxPairEndAlignmentScore(), scores.getPairEndScore()) + params.getSELF_ALIGNMENT_WEIGHT() * set.getMaxSelfAlignmentScore() + params.getSELF_END_ALIGNMENT_WEIGHT() * set.getMaxSelfEndAlignmentScore();
		assertEquals(expect, set.scorePrimerSetTo(primer));
		
		// score primer to primer set containing two instances of itself
		primer = new Primer("ATATATATATAT", PrimerTypes.forwardPrimer, params);
		set = new PrimerSet(params);
		set.addForwardPrimer(primer);
		set.addForwardPrimer(primer);
		scores = alignment.computePairAlignment(primer.getSequence(), set.getForwardPrimer(0).getSequence());
		expect = params.getPAIR_ALIGNMENT_WEIGHT() * Math.max(set.getMaxPairAlignmentScore(), scores.getPairScore()) + params.getPAIR_END_ALIGNMENT_WEIGHT() * Math.max(set.getMaxPairEndAlignmentScore(), scores.getPairEndScore()) + params.getSELF_ALIGNMENT_WEIGHT() * set.getMaxSelfAlignmentScore() + params.getSELF_END_ALIGNMENT_WEIGHT() * set.getMaxSelfEndAlignmentScore();
		assertEquals(expect, set.scorePrimerSetTo(primer));
	}

}
