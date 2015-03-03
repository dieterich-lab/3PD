package primerDesign.Test;

import java.util.Arrays;

import junit.framework.TestCase;

public class EnhancedSuffixArrayIntOldTest extends TestCase {

	public void testComputeCombinedChildTable() {
		EnhancedSuffixArrayIntOld array = new EnhancedSuffixArrayIntOld("acaaacatat");
		array.createIndex(Integer.MAX_VALUE, true);
		int[] expect = array.getChildTab();
		
		array.initChildTab();
		array.computeCombinedChildTable();
		int[] real = array.getChildTab();
		
		assertTrue(Arrays.equals(expect, real));
	}

}
