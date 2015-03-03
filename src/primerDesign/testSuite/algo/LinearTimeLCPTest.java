/**
 * 
 */
package primerDesign.testSuite.algo;

import static org.junit.Assert.assertArrayEquals;

import java.io.File;
import java.io.IOException;

import junit.framework.TestCase;

import org.junit.Test;

import primerDesign.algo.KoAluruSuffixSort;
import primerDesign.algo.LinearTimeLCP;
import primerDesign.util.SlimFastaParser;

/**
 * This unit test checks the proper implementation of the speedup method "LinearTimeLCP".
 * 
 * This method should lead to the same result as the non-linear time method proposed by Kurtz et.al. in their enhanced suffix arrays paper.
 * 
 * @author froehler
 *
 */
public class LinearTimeLCPTest extends TestCase{

	/**
	 * Test method for {@link primerDesign.algo.LinearTimeLCP#getLCP(char[], int[])}.
	 * @throws IOException 
	 */
	@Test
	public void testGetLCP() throws IOException {
		//SlimFastaParser parser = new SlimFastaParser(new File("Testsequenzen/Ppa/BAC-Ppa50-C09.fa"));
		SlimFastaParser parser = new SlimFastaParser(new File("Testsequenzen/Ppa/Ppacificus-unmasked.fa"));
		
		char[] text = parser.parseNextContigIgnoreCase().getSequence();
		//text = "acaaacatat".toCharArray();
		KoAluruSuffixSort sort = new KoAluruSuffixSort();
		int[] suffixArray = sort.getSuffixArray(text);
		int[] linearTimeLcp = LinearTimeLCP.getLCP(text, suffixArray);
		int[] nonLinearTimeLcp = computeLCPTable(suffixArray, text);
		
		assertTrue(linearTimeLcp.length == nonLinearTimeLcp.length);
		assertArrayEquals(nonLinearTimeLcp, linearTimeLcp);
	}
	
	private final int[] computeLCPTable(int[] suftab, char[] sequence) {
		int[] lcpTab = new int[suftab.length];
		int lcp;
		int stop = suftab.length;
		for (int i = 0; i < stop; i++) {
			if(i == 0) lcpTab[i] = 0;
			else{
				lcp = computeLcp(suftab, i, i-1, sequence, sequence.length);
				lcpTab[i] = lcp;
			}
		}
		return lcpTab;
	}
	
	private final int computeLcp(int[] suftab, int i, int j, char[] sequence, int sequenceLength) {
		int idxI = suftab[i];
		int idxJ = suftab[j];
		int end = sequenceLength - Math.max(idxI, idxJ);
		int lcp = 0;
		for(int k = 0; k < end; k++) {
			if(sequence[idxI + k] == sequence[idxJ + k]) lcp++;
			else break;
		}
		return lcp;
	}

}
