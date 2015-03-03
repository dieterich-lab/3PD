/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;
import primerDesign.algo.PrimerAlignmentCalculation;
import primerDesign.algo.SimpleAlignment;
import primerDesign.dsc.PrimerAlignmentScores;

/**
 * @author froehler
 *
 */
public class SimpleAlignmentTest extends TestCase {

	/**
	 * Test method for {@link primerDesign.algo.SimpleAlignment#computePairAlignment(java.lang.String, java.lang.String)}.
	 */
	public void testComputePairAlignment() {
		PrimerAlignmentCalculation alignment = new SimpleAlignment(2, 4);
		PrimerAlignmentScores scoresSelf;
		
		scoresSelf = alignment.computeSelfAlignment("TATA");
		assertEquals(8, scoresSelf.getPairScore());
		assertEquals(8, scoresSelf.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("GCGC");
		assertEquals(16, scoresSelf.getPairScore());
		assertEquals(16, scoresSelf.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("taat");
		assertEquals(4, scoresSelf.getPairScore());
		assertEquals(4, scoresSelf.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("GGAA");
		assertEquals(0, scoresSelf.getPairScore());
		assertEquals(0, scoresSelf.getPairEndScore());
		
//		scoresSelf = alignment.computeSelfAlignment("GGCC");
//		assertEquals(16, scoresSelf.getPairScore());
//		assertEquals(0, scoresSelf.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("ATGC");
		assertEquals(8, scoresSelf.getPairScore());
		assertEquals(8, scoresSelf.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("GGGATCGGGAT");
		assertEquals(16, scoresSelf.getPairScore());
		assertEquals(4, scoresSelf.getPairEndScore());
	}

	/**
	 * Test method for {@link primerDesign.algo.SimpleAlignment#computeSelfAlignment(java.lang.String)}.
	 */
	public void testComputeSelfAlignment() {
		PrimerAlignmentCalculation alignment = new SimpleAlignment(2, 4);
		PrimerAlignmentScores scoresSelf;
		PrimerAlignmentScores scoresPair;
		
		scoresSelf = alignment.computeSelfAlignment("TATA");
		scoresPair = alignment.computePairAlignment("TATA","TATA");
		assertEquals(8, scoresPair.getPairScore());
		assertEquals(8, scoresPair.getPairEndScore());
		assertEquals(scoresSelf.getPairScore(), scoresPair.getPairScore());
		assertEquals(scoresSelf.getPairEndScore(), scoresPair.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("GCGC");
		scoresPair = alignment.computePairAlignment("GCGC","GCGC");
		assertEquals(16, scoresPair.getPairScore());
		assertEquals(16, scoresPair.getPairEndScore());
		assertEquals(scoresSelf.getPairScore(), scoresPair.getPairScore());
		assertEquals(scoresSelf.getPairEndScore(), scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("taat","taat");
		scoresSelf = alignment.computeSelfAlignment("taat");
		assertEquals(4, scoresPair.getPairScore());
		assertEquals(4, scoresPair.getPairEndScore());
		assertEquals(scoresSelf.getPairScore(), scoresPair.getPairScore());
		assertEquals(scoresSelf.getPairEndScore(), scoresPair.getPairEndScore());
		
		scoresSelf = alignment.computeSelfAlignment("TGATCGGGAA");
		scoresPair = alignment.computePairAlignment("TGATCGGGAA","TGATCGGGAA");
		assertEquals(12, scoresPair.getPairScore());
		assertEquals(2, scoresPair.getPairEndScore());
		assertEquals(scoresSelf.getPairScore(), scoresPair.getPairScore());
		assertEquals(scoresSelf.getPairEndScore(), scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("TGAT","TGATC");
		assertEquals(8, scoresPair.getPairScore());
		assertEquals(8, scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("TGATC","TGAT");
		assertEquals(8, scoresPair.getPairScore());
		assertEquals(8, scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("TGATCGGG","TGAT");
		assertEquals(8, scoresPair.getPairScore());
		assertEquals(0, scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("GGGTGATC","TGAT");
		assertEquals(8, scoresPair.getPairScore());
		assertEquals(8, scoresPair.getPairEndScore());
		
		scoresPair = alignment.computePairAlignment("ACCCATGGCAGGTTCACCC", "ACCCATGGCAGGTTCACCC");
		assertEquals(24, scoresPair.getPairScore());
		assertEquals(8, scoresPair.getPairEndScore());
	}
}
