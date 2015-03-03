package primerDesign.testSuite.algo;

import junit.framework.TestCase;
import primerDesign.algo.KaempkePrimerAlignment;
import primerDesign.algo.PrimerAlignmentCalculation;

public class KaempkePrimerAlignmentTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	public void testDynsaAPP() {
		PrimerAlignmentCalculation alignment = new KaempkePrimerAlignment(null);
		
		alignment.computeSelfAlignment("TATA");
		assertEquals(8, alignment.getSAScore());
		assertEquals(8, alignment.getSEAScore());
		
		alignment.computeSelfAlignment("GCGC");
		assertEquals(16, alignment.getSAScore());
		assertEquals(16, alignment.getSEAScore());
		
		alignment.computeSelfAlignment("taat");
		assertEquals(4, alignment.getSAScore());
		assertEquals(4, alignment.getSEAScore());
		
		alignment.computeSelfAlignment("TGATCGGGAA");
		assertEquals(12, alignment.getSAScore());
		assertEquals(2, alignment.getSEAScore());
		
		alignment.computeSelfAlignment("ACCCATGGCAGGTTCACCC");
		assertEquals(24, alignment.getSAScore());
		assertEquals(8, alignment.getSEAScore());
	}

	public void testDynpaAPP() {
		PrimerAlignmentCalculation alignment = new KaempkePrimerAlignment(null);
		
		alignment.computeSelfAlignment("TATA");
		alignment.computePairAlignment("TATA","TATA");
		assertEquals(8, alignment.getPAScore());
		assertEquals(8, alignment.getPEAScore());
		assertEquals(alignment.getSAScore(), alignment.getPAScore());
		assertEquals(alignment.getSEAScore(), alignment.getPEAScore());
		
		alignment.computeSelfAlignment("GCGC");
		alignment.computePairAlignment("GCGC","GCGC");
		assertEquals(16, alignment.getPAScore());
		assertEquals(16, alignment.getPEAScore());
		assertEquals(alignment.getSAScore(), alignment.getPAScore());
		assertEquals(alignment.getSEAScore(), alignment.getPEAScore());
		
		alignment.computePairAlignment("taat","taat");
		alignment.computeSelfAlignment("taat");
		assertEquals(4, alignment.getPAScore());
		assertEquals(4, alignment.getPEAScore());
		assertEquals(alignment.getSAScore(), alignment.getPAScore());
		assertEquals(alignment.getSEAScore(), alignment.getPEAScore());
		
		alignment.computeSelfAlignment("TGATCGGGAA");
		alignment.computePairAlignment("TGATCGGGAA","TGATCGGGAA");
		assertEquals(12, alignment.getPAScore());
		assertEquals(2, alignment.getPEAScore());
		assertEquals(alignment.getSAScore(), alignment.getPAScore());
		assertEquals(alignment.getSEAScore(), alignment.getPEAScore());
		
		alignment.computePairAlignment("TGAT","TGATC");
		assertEquals(8, alignment.getPAScore());
		assertEquals(8, alignment.getPEAScore());
		
		alignment.computePairAlignment("TGATC","TGAT");
		assertEquals(8, alignment.getPAScore());
		assertEquals(8, alignment.getPEAScore());
		
		alignment.computePairAlignment("TGATCGGG","TGAT");
		assertEquals(8, alignment.getPAScore());
		assertEquals(0, alignment.getPEAScore());
		
		alignment.computePairAlignment("GGGTGATC","TGAT");
		assertEquals(8, alignment.getPAScore());
		assertEquals(8, alignment.getPEAScore());
		
		alignment.computePairAlignment("ACCCATGGCAGGTTCACCC", "ACCCATGGCAGGTTCACCC");
		assertEquals(24, alignment.getSAScore());
		assertEquals(8, alignment.getSEAScore());
	}
}
