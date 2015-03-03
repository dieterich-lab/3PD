/**
 * 
 */
package primerDesign.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

import primerDesign.dsc.AlignmentType;
import primerDesign.util.PrimerSearchParameters;

/**
 * @author Sebastian Fršhler
 *
 */
public class PrimerPairSetAlignmentsTest {


	@Test
	public void testAddScanRegion() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		PrimerPairSetAlignments alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		try{
			alignments.getAlignment(0, 0, 1, 1, 0, 1, AlignmentType.forward1Forward2);
			fail("Exception expected!");
		}
		catch(ArrayIndexOutOfBoundsException e){
			
		}
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		assertEquals(0, alignments.getAlignment(0, 0, 1, 1, 0, 1, AlignmentType.forward1Forward2).getPairScore());
	}

	/**
	 * Test method for {@link primerDesign.algo.PrimerPairSetAlignments#deleteLastRegion()}.
	 */
	@Test
	public void testDeleteLastRegion() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		PrimerPairSetAlignments alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward2Probe1).getPairScore());
		
		alignments.deleteLastRegion();
		try{
			alignments.getAlignment(0, 0, 1, 1, 0, 1, AlignmentType.forward1Forward2);
			fail("Exception expected!");
		}
		catch(ArrayIndexOutOfBoundsException e){
			
		}
	}

	/**
	 * Test method for {@link primerDesign.algo.PrimerPairSetAlignments#getAlignment(int, int, int, int, int, int, java.lang.Enum)}.
	 */
	@Test
	public void testGetAlignment() {
		PrimerSearchParameters params = new PrimerSearchParameters();
		PrimerPairSetAlignments alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		alignments.addScanRegion("ATGC".toCharArray(), "ATGC".toCharArray(), params);
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(8,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward2Probe1).getPairScore());
		
		alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("AAAA".toCharArray(), "TTTT".toCharArray(), params);
		alignments.addScanRegion("GGGG".toCharArray(), "CCCC".toCharArray(), params);
		
		assertEquals(0,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(0,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward2Probe1).getPairScore());
		
		alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("ATTA".toCharArray(), "CGGC".toCharArray(), params);
		alignments.addScanRegion("GCCG".toCharArray(), "TAAT".toCharArray(), params);
		
		assertEquals(0,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(16,alignments.getAlignment(0, 0, 4, 1, 0, 4, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(0,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(2,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 1, 3, 1, 1, 3, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward2Probe1).getPairScore());
		
		assertEquals(0,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(2,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 1, 2, 1, 1, 2, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 2, 1, 1, 2, 1, AlignmentType.forward2Probe1).getPairScore());
		
		alignments = new PrimerPairSetAlignments(2);
		alignments.addScanRegion("ATTAATTA".toCharArray(), "CGGCCGG".toCharArray(), params);
		alignments.addScanRegion("GCCG".toCharArray(), "TAATAT".toCharArray(), params);
		
		assertEquals(0,alignments.getAlignment(0, 0, 8, 1, 0, 4, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 8, 1, 0, 6, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(8,alignments.getAlignment(0, 0, 8, 1, 0, 6, AlignmentType.forward1Probe2).getPairEndScore());
		assertEquals(16,alignments.getAlignment(0, 0, 4, 1, 0, 7, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(12,alignments.getAlignment(0, 0, 4, 1, 0, 7, AlignmentType.forward2Probe1).getPairEndScore());
		
		assertEquals(0,alignments.getAlignment(0, 2, 4, 1, 1, 3, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(0,alignments.getAlignment(0, 2, 2, 1, 2, 2, AlignmentType.forward1Forward2).getPairScore());
		assertEquals(6,alignments.getAlignment(0, 2, 4, 1, 1, 5, AlignmentType.forward1Probe2).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 2, 4, 1, 1, 5, AlignmentType.forward1Probe2).getPairEndScore());
		assertEquals(12,alignments.getAlignment(0, 0, 3, 1, 1, 5, AlignmentType.forward2Probe1).getPairScore());
		assertEquals(4,alignments.getAlignment(0, 0, 3, 1, 1, 5, AlignmentType.forward2Probe1).getPairEndScore());		
	}
}
