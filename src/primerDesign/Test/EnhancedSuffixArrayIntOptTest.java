package primerDesign.Test;

import java.util.Arrays;

import junit.framework.TestCase;
import cern.colt.list.IntArrayList;

public class EnhancedSuffixArrayIntOptTest extends TestCase {

	public void testQuicksortChildmap() {
		EnhancedSuffixArrayIntOpt array = new EnhancedSuffixArrayIntOpt("abc");
		IntArrayList map = new IntArrayList();
		map.add(0); map.add(1);
		map.add(10); map.add(8);
		map.add(6); map.add(7);
		map.add(4); map.add(4);
		array.setChildMap(map);
		array.quicksortChildmap(0, map.size()-1);
		IntArrayList result = array.getChildMap();
		
		IntArrayList expect = new IntArrayList();
		expect.add(0); expect.add(1);
		expect.add(4); expect.add(4);
		expect.add(6); expect.add(7);
		expect.add(10); expect.add(8);
		
		assertEquals(expect, result);
	}
	
	public void testChildTabComp(){
		EnhancedSuffixArrayIntOpt array = new EnhancedSuffixArrayIntOpt("acaaacatat");
		array.createIndex(Integer.MAX_VALUE, true);
		array.initChildTab();
		array.computeChildNextIndex();
		byte[] next = array.getChildTab();
		byte empty = Byte.MIN_VALUE;
		byte[] expect = new byte[]{-122, empty, -126, empty, empty, empty, -126, empty, -126, empty, empty};
		assertTrue(Arrays.equals(next, expect));
		
		array.computeChildUpDown();
		byte[] all = array.getChildTab();
		expect = new byte[]{-122, empty, -126, empty, -127, Byte.MAX_VALUE, -126, empty, -126, empty, empty};
		assertTrue(Arrays.equals(all, expect));
		
		array.initChildTab();
		array.computeCombinedChildTable();
		all = array.getChildTab();
		assertTrue(Arrays.equals(all, expect));
	}
	
	public void testSearch(){
		EnhancedSuffixArrayIntOpt array = new EnhancedSuffixArrayIntOpt("acaaacatat");
		array.createIndex(Integer.MAX_VALUE, true);
		array.findMatchPositions("ACA");
	}

}
