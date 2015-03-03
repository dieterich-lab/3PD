/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;

import org.junit.Test;

import primerDesign.algo.SequenceRegionAligner;
import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.SequenceRegionAlignment;

/**
 * @author froehler
 *
 */
public class SequenceRegionAlignerTest extends TestCase{

	/**
	 * Test method for {@link primerDesign.algo.SequenceRegionAligner#alignSequenceRegions(char[], char[])}.
	 */
	@Test
	public void testAlignSequenceRegions() {
		char[] first = "ATGC".toCharArray();
		char[] second = "GCAT".toCharArray();
		SequenceRegionAlignment alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		PrimerAlignmentScores scores;
		
		// test the whole alignment
		assertEquals(12, alignment.getGlobalAlignmentValues(0, 4, 0, 4).getPairScore(), 0.001);
		
		// test the 'center-box'-alignment
		assertEquals(6, alignment.getGlobalAlignmentValues(1, 2, 1, 2).getPairScore(), 0.001);
		
		// test each of the four corner alignments
		assertEquals(4, alignment.getGlobalAlignmentValues(0, 2, 2, 2).getPairScore(), 0.001);
		
		assertEquals(0, alignment.getGlobalAlignmentValues(0, 2, 0, 2).getPairScore(), 0.001);
		
		assertEquals(0, alignment.getGlobalAlignmentValues(2, 2, 2, 2).getPairScore(), 0.001);
		
		assertEquals(8, alignment.getGlobalAlignmentValues(2, 2, 0, 2).getPairScore(), 0.001);
		
		// test single alignments
		assertEquals(2, alignment.getGlobalAlignmentValues(0, 1, 3, 1).getPairScore(), 0.001);
		
		assertEquals(0, alignment.getGlobalAlignmentValues(0, 1, 0, 1).getPairScore(), 0.001);
		
		assertEquals(0, alignment.getGlobalAlignmentValues(3, 1, 3, 1).getPairScore(), 0.001);
		
		assertEquals(4, alignment.getGlobalAlignmentValues(3, 1, 0, 1).getPairScore(), 0.001);
		
		first = "CTTGACTATTAACTCACTTGTAGT".toCharArray();
		second = "TCCCTAATCTTTCATGATCTA".toCharArray();
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		scores = alignment.getGlobalAlignmentValues(0, first.length, 0, second.length);
		assertEquals(22, scores.getPairScore());
		assertEquals(8, scores.getPairEndScore());
		
		first = "TCCCTAATCTTTCATGATCTA".toCharArray();
		second = "CTTGACTATTAACTCACTTGTAGT".toCharArray();
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		scores = alignment.getGlobalAlignmentValues(0, first.length, 0, second.length);
		assertEquals(22, scores.getPairScore());
		assertEquals(8, scores.getPairEndScore());
		
		first = "aaaaaaaaaTCCCTAATCTTTCATGATCTAttttttttt".toCharArray();
		second = "aaaaaaaaaCTTGACTATTAACTCACTTGTAGTttttttttt".toCharArray();
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		scores = alignment.getGlobalAlignmentValues(9, 21, 9, 24);
		assertEquals(22, scores.getPairScore());
		assertEquals(8, scores.getPairEndScore());
		
		first = "TCCCTAATCTTTCATGATCTAaaaaaaaaa".toCharArray();
		second = "aaaaaaaaaCTTGACTATTAACTCACTTGTAGTaaaaaaaaa".toCharArray();
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		scores = alignment.getGlobalAlignmentValues(0, 21, 9, 24);
		assertEquals(22, scores.getPairScore());
		assertEquals(8, scores.getPairEndScore());
		
		first = "aaaaaaaaaTCCCTAATCTTTCATGATCTA".toCharArray();
		second = "aaaaaaaaaCTTGACTATTAACTCACTTGTAGTttttttttt".toCharArray();
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second, 2, 4);
		scores = alignment.getGlobalAlignmentValues(9, 21, 9, 24);
		assertEquals(22, scores.getPairScore());
		assertEquals(8, scores.getPairEndScore());
	}
}
