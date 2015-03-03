package primerDesign.testSuite.dsc;


import java.util.Arrays;
import java.util.HashSet;

import junit.framework.TestCase;
import primerDesign.dsc.DNASuffixTrieWithPositions;

public class DNASuffixTrieWithPositionsTest extends TestCase {

	public void testSearchMatchPositionsInIndex() {
		DNASuffixTrieWithPositions tree = new DNASuffixTrieWithPositions();
		tree.createIndex("ATGCNATGCN", 5, true);
		Integer[] matches = tree.searchMatchPositionsInIndex("ATGCN");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{0,5}, matches));
		assertTrue(Arrays.equals(new Integer[]{0}, tree.searchMatchPositionsInIndex("ATGCNATGCN")));
		assertTrue(Arrays.equals(new Integer[]{1}, tree.searchMatchPositionsInIndex("TGCNATGCN")));
		assertTrue(Arrays.equals(new Integer[]{2}, tree.searchMatchPositionsInIndex("GCNATGCN")));
		assertTrue(Arrays.equals(new Integer[]{3}, tree.searchMatchPositionsInIndex("CNATGCN")));
		assertTrue(Arrays.equals(new Integer[]{4}, tree.searchMatchPositionsInIndex("NATGCN")));
		matches = tree.searchMatchPositionsInIndex("ATGCN");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{0,5}, matches));
		matches = tree.searchMatchPositionsInIndex("TGCN");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{1,6}, matches));
		matches = tree.searchMatchPositionsInIndex("GCN");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{2,7}, matches));
		matches = tree.searchMatchPositionsInIndex("CN");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{3,8}, matches));
		matches = tree.searchMatchPositionsInIndex("N");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{4,9}, matches));
		matches = tree.searchMatchPositionsInIndex("GCNAT");
		Arrays.sort(matches);
		assertTrue(Arrays.equals(new Integer[]{2}, matches));
	}

	public void testSearchNbMatchesInIndex() {
		DNASuffixTrieWithPositions tree = new DNASuffixTrieWithPositions();
		tree.createIndex("ATGCNATGCN", 5, true);
		assertEquals(2, tree.searchNbMatchesInIndex("ATGCN"));
		assertEquals(2, tree.searchNbMatchesInIndex("TGCN"));
		assertEquals(2, tree.searchNbMatchesInIndex("GCN"));
		assertEquals(2, tree.searchNbMatchesInIndex("CN"));
		assertEquals(2, tree.searchNbMatchesInIndex("N"));
		
		// string longer than max word size
		assertEquals(1, tree.searchNbMatchesInIndex("ATGCNATGCN"));
		// string longer than index string!
		assertEquals(0, tree.searchNbMatchesInIndex("ATGCNATGCNN"));
	}

	public void testCompareMatchSets(){
		HashSet<Integer> first = new HashSet<Integer>();
		first.add(0);
		first.add(5);
		HashSet<Integer> second = new HashSet<Integer>();
		second.add(4);
		HashSet<Integer> expected = new HashSet<Integer>();
		expected.add(0);
		assertEquals(expected, DNASuffixTrieWithPositions.compareMatchSets(first, second, 4));
		second.add(10);
		assertEquals(expected, DNASuffixTrieWithPositions.compareMatchSets(first, second, 4));
		second.add(9);
		expected.add(5);
		assertEquals(expected, DNASuffixTrieWithPositions.compareMatchSets(first, second, 4));
	}
}
