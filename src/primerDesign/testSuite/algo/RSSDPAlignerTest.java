/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;
import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.algo.RSSDPAligner;
import primerDesign.algo.SimpleAlignment;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;

/**
 * @author froehler
 *
 */
public class RSSDPAlignerTest extends TestCase {
	static PrimerSearchParameters params = new PrimerSearchParameters();
	
	private static PrimerAlignmentCalculation reference = new SimpleAlignment(params.getA_t_basepair_score(), params.getG_c_basepair_score());
	private static RSSDPAligner aligner;
	private String primerSequence;
	private Primer primer;

	/**
	 * Test method for {@link primerDesign.algo.RSSDPAligner#computeSelfAlignment(primerDesign.dsc.Primer)}.
	 */
	public void testComputeSelfAlignment() {
		int iter = 10000;
		
		runTest(iter, PrimerTypes.forwardPrimer);
		runTest(iter, PrimerTypes.reversePrimer);
		runTest(iter, PrimerTypes.hybridizationProbe);
	}
	
	private void runTest(int iter, PrimerTypes primerType){
		for(int i=0; i<iter; i++){
			primerSequence = SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH());
			aligner = new RSSDPAligner(primerSequence, 0, params);
			primer = new Primer(primerSequence, primerType, aligner, 0, params);
			primer.setRelativePosition(0);
			assertEquals(reference.computeSelfAlignment(primerSequence).getPairScore(), aligner.computeSelfAlignment(primer).getPairScore());
			assertEquals(reference.computeSelfAlignment(primerSequence).getPairEndScore(), aligner.computeSelfAlignment(primer).getPairEndScore());
			
			aligner = new RSSDPAligner("nnn" + primerSequence + "nn", 0, params);
			primer.setRelativePosition(3);
			assertEquals(reference.computeSelfAlignment(primerSequence).getPairScore(), aligner.computeSelfAlignment(primer).getPairScore());
			assertEquals(reference.computeSelfAlignment(primerSequence).getPairEndScore(), aligner.computeSelfAlignment(primer).getPairEndScore());
		}
	}
}
