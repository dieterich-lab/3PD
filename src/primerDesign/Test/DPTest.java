/**
 * 
 */
package primerDesign.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * @author froehler
 *
 */
public class DPTest {

	/**
	 * Test method for {@link shortReadMapping.test.DP#getDP(char[], int, int, char[], int, int)}.
	 */
	@Test
	public void testGetDP() {
		char[] sequence = "atgc".toCharArray();
		char[] query = "atgc".toCharArray();
		assertEquals(0, DP.getDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		assertEquals(0, DP.getDP(sequence, 1, sequence.length - 1, query, 1, query.length - 1));
		assertEquals(0, DP.getDP(sequence, 2, sequence.length - 1, query, 2, query.length - 1));
		assertEquals(0, DP.getDP(sequence, 3, sequence.length - 1, query, 3, query.length - 1));
		
		assertEquals(0, DP.getDP(sequence, 0, sequence.length - 2, query, 0, query.length - 2));
		assertEquals(0, DP.getDP(sequence, 0, sequence.length - 3, query, 0, query.length - 3));
		assertEquals(0, DP.getDP(sequence, 0, sequence.length - 4, query, 0, query.length - 4));
		
		sequence = "aaaaaaaaaatttttttttt".toCharArray();
		query = "aaattt".toCharArray();
		
		assertEquals(14, DP.getDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		assertEquals(0, DP.getDP(sequence, 7, sequence.length - 8, query, 0, query.length - 1));
		assertEquals(1, DP.getDP(sequence, 6, sequence.length - 9, query, 0, query.length - 1));
		assertEquals(1, DP.getDP(sequence, 8, sequence.length - 7, query, 0, query.length - 1));
		
		sequence = "tttttttttt".toCharArray();
		query = "aaaaaaaaaa".toCharArray();
		assertEquals(10, DP.getDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "".toCharArray();
		query = "aaaaaaaaaa".toCharArray();
		
		assertEquals(10, DP.getDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "aaaaaaaaaa".toCharArray();
		query = "".toCharArray();
		
		assertEquals(10, DP.getDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
	}
	
	@Test
	public void testGetRevcompDP(){
		char[] sequence = "atgc".toCharArray();
		char[] query = "gcat".toCharArray();
		assertEquals(0, DP.getRevcompDP(query, 0, query.length - 1, sequence, 0, sequence.length - 1));
		assertFalse(DP.getRevcompDP(sequence, 0, sequence.length - 1, sequence, 0, sequence.length - 1) == 0);
		assertTrue(DP.getRevcompDP(sequence, 0, sequence.length - 1, sequence, 0, sequence.length - 1) == 4);
		
		assertEquals(0, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		assertEquals(0, DP.getRevcompDP(sequence, 1, sequence.length - 1, query, 0, query.length - 2));
		assertEquals(0, DP.getRevcompDP(sequence, 2, sequence.length - 1, query, 0, query.length - 3));
		assertEquals(0, DP.getRevcompDP(sequence, 3, sequence.length - 1, query, 0, query.length - 4));
		
		assertEquals(0, DP.getRevcompDP(sequence, 0, sequence.length - 2, query, 1, query.length - 1));
		assertEquals(0, DP.getRevcompDP(sequence, 0, sequence.length - 3, query, 2, query.length - 1));
		assertEquals(0, DP.getRevcompDP(sequence, 0, sequence.length - 4, query, 3, query.length - 1));
		
		sequence = "aaaaaaaaaatttttttttt".toCharArray();
		query = "aaattt".toCharArray();
		
		assertEquals(14, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		assertEquals(0, DP.getRevcompDP(sequence, 7, sequence.length - 8, query, 0, query.length - 1));
		assertEquals(1, DP.getRevcompDP(sequence, 6, sequence.length - 9, query, 0, query.length - 1));
		assertEquals(1, DP.getRevcompDP(sequence, 8, sequence.length - 7, query, 0, query.length - 1));
		
		sequence = "tttttttttt".toCharArray();
		query = "aaaaaaaaaa".toCharArray();
		assertEquals(0, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "tttttttttt".toCharArray();
		query = "tttttttttt".toCharArray();
		assertEquals(10, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "tttttttttt".toCharArray();
		query = "gggggggggg".toCharArray();
		assertEquals(10, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "tttttttttt".toCharArray();
		query = "cccccccccc".toCharArray();
		assertEquals(10, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "".toCharArray();
		query = "aaaaaaaaaa".toCharArray();
		
		assertEquals(10, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
		
		sequence = "aaaaaaaaaa".toCharArray();
		query = "".toCharArray();
		
		assertEquals(10, DP.getRevcompDP(sequence, 0, sequence.length - 1, query, 0, query.length - 1));
	}
}
