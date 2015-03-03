/**
 * 
 */
package primerDesign.algo;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

/**
 * @author Sebastian Fršhler
 *
 */
public class KoAluruSuffixSortTest {

	/**
	 * Test method for {@link primerDesign.algo.KoAluruSuffixSort#getSuffixArray(char[])}.
	 */
	@Test
	public void testGetSuffixArray() {
		char[] text;
		int[] expect;
		int[] result;
		
		KoAluruSuffixSort sorter = new KoAluruSuffixSort();
		
		text = "abc".toCharArray();
		expect = new int[]{3,0,1,2};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "cba".toCharArray();
		expect = new int[]{3,2,1,0};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "acb".toCharArray();
		expect = new int[]{3,0,2,1};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "abcd".toCharArray();
		expect = new int[]{4,0,1,2,3};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "dcba".toCharArray();
		expect = new int[]{4,3,2,1,0};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "adcb".toCharArray();
		expect = new int[]{4, 0,3,2,1};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		text = "adbc".toCharArray();
		expect = new int[]{4, 0,2,3,1};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		// majority: type-S suffix word
		text = "mississippi".toCharArray();
		expect = new int[]{11, 10,7,4,1,0,9,8,6,3,5,2};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
		
		// majority: type-L suffix word
		text = "acaaacatat".toCharArray();
		expect = new int[]{10, 2, 3, 0, 4, 8, 6, 1, 5, 9, 7};
		result = sorter.getSuffixArray(text);
		
		for(int i=0; i<expect.length; i++) assertEquals(expect[i], result[i]);
	}
}
