package primerDesign.testSuite;

import junit.framework.Test;
import junit.framework.TestSuite;
import primerDesign.testSuite.algo.LinearTimeLCPTest;
import primerDesign.testSuite.algo.PrimerMisprimingTest;
import primerDesign.testSuite.algo.PrimerPairSetAlignmentsTest;
import primerDesign.testSuite.algo.PrimerSearchTest;
import primerDesign.testSuite.algo.RestrictionSiteSearchTest;
import primerDesign.testSuite.algo.SantaLuciaTMTest;
import primerDesign.testSuite.algo.SequenceRegionAlignerTest;
import primerDesign.testSuite.algo.ShortFragmentExcluderTest;
import primerDesign.testSuite.algo.SimpleAlignmentTest;
import primerDesign.testSuite.dsc.DNASuffixTreeTest;
import primerDesign.testSuite.dsc.DNASuffixTrieWithPositionsTest;
import primerDesign.testSuite.dsc.PrimerSetTest;
import primerDesign.testSuite.dsc.SequenceRegionTest;

public class AllTests {

	public static Test suite() {
		TestSuite suite = new TestSuite("Test for primerDesign.testSuite");
		//$JUnit-BEGIN$
		
		// algo
		//suite.addTestSuite(KaempkePrimerAlignmentTest.class);
		suite.addTestSuite(LinearTimeLCPTest.class);
		suite.addTestSuite(PrimerMisprimingTest.class);
		suite.addTestSuite(PrimerPairSetAlignmentsTest.class);
		suite.addTestSuite(PrimerSearchTest.class);
		suite.addTestSuite(RestrictionSiteSearchTest.class);
		//suite.addTestSuite(RSSDPAlignerTest.class);		
		suite.addTestSuite(SantaLuciaTMTest.class);
		suite.addTestSuite(SequenceRegionAlignerTest.class);
		suite.addTestSuite(ShortFragmentExcluderTest.class);
		suite.addTestSuite(SimpleAlignmentTest.class);
		//suite.addTestSuite(SimpleGreedyPrimerPairPickingTest.class);
		
		// dsc
		suite.addTestSuite(DNASuffixTreeTest.class);
		suite.addTestSuite(DNASuffixTrieWithPositionsTest.class);
		suite.addTestSuite(PrimerSetTest.class);
		suite.addTestSuite(SequenceRegionTest.class);
		
		//$JUnit-END$
		return suite;
	}

}
