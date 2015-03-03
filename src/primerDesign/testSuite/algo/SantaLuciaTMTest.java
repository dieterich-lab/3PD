package primerDesign.testSuite.algo;

import junit.framework.TestCase;
import primerDesign.algo.SantaLuciaTM;
import primerDesign.util.PrimerSearchParameters;

public class SantaLuciaTMTest extends TestCase {

//	public void testComputeTM() {
//		fail("Not yet implemented");
//	}

	public void testComputeCationCorrectedTM() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		params.setPRIMER_CONCENTRATION(50e-9);
		params.setMONOVALENT_CATION_CONCENTRATION(50e-3);
		params.setDIVALENT_CATION_CONCENTRATION(0);
		params.setDNTP_CONCENTRATION(0);
		
		double tolerance = 1.0E-06;  // the maximum difference of results between Primer3 T_m and this implementation
		
		SantaLuciaTM lucia = new SantaLuciaTM();
		// check for same result between Primer3 and this implementation of SantaLuciaTM
		assertEquals(-88.389940, lucia.computeCationCorrectedTM("ATAT", params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION()), tolerance);
		assertEquals(-22.546940, lucia.computeCationCorrectedTM("GCGC", params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION()), tolerance);
		assertEquals(38.245134, lucia.computeCationCorrectedTM("ATATATATATATATATATATATATATATATATATAT", params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION()), tolerance);
		assertEquals(90.305036, lucia.computeCationCorrectedTM("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION()), tolerance);

		// check for exception if primer length > 36 bp (formular is not valid from this length on!)
	}

}
