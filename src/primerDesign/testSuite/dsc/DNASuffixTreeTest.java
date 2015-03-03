package primerDesign.testSuite.dsc;

import junit.framework.TestCase;
import primerDesign.dsc.DNASuffixTrie;
import primerDesign.util.Constants;
import cern.colt.list.ObjectArrayList;

public class DNASuffixTreeTest extends TestCase {

	public void testDNASuffixTree_newObjectArrayList() {
		ObjectArrayList alphabet = new ObjectArrayList();
		alphabet.add((char) 'A');
		alphabet.add((char) 'T');
		alphabet.add((char) 'G');
		alphabet.add((char) 'C');
		alphabet.add((char) Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase().charAt(0));
		
		DNASuffixTrie tree = new DNASuffixTrie(alphabet);
		assertEquals(new DNASuffixTrie(), tree);
	}

//
//	public void testCreateIndex() {
//		fail("Not yet implemented");
//	}
//
	public void testSearchNbMatchesInIndex() {
		ObjectArrayList alphabet = new ObjectArrayList();
		alphabet.add((char) 'A');
		alphabet.add((char) 'T');
		alphabet.add((char) 'G');
		alphabet.add((char) 'C');
		alphabet.add((char) Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase().charAt(0));
		
		DNASuffixTrie tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgc", 4, true);
		assertEquals(0, tree.searchNbMatchesInIndex("n"));
		
		tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgcn", 5, true);
		assertEquals(1, tree.searchNbMatchesInIndex("atgcn"));
		assertEquals(1, tree.searchNbMatchesInIndex("tgcn"));
		assertEquals(1, tree.searchNbMatchesInIndex("gcn"));
		assertEquals(1, tree.searchNbMatchesInIndex("cn"));
		assertEquals(1, tree.searchNbMatchesInIndex("t"));
		assertEquals(1, tree.searchNbMatchesInIndex("g"));
		assertEquals(1, tree.searchNbMatchesInIndex("c"));
		assertEquals(1, tree.searchNbMatchesInIndex("n"));
		
		tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgcnatgcn", 5, true);
		assertEquals(2, tree.searchNbMatchesInIndex("atgcn"));
		assertEquals(2, tree.searchNbMatchesInIndex("tgcn"));
		assertEquals(2, tree.searchNbMatchesInIndex("gcn"));
		assertEquals(2, tree.searchNbMatchesInIndex("cn"));
		assertEquals(2, tree.searchNbMatchesInIndex("t"));
		assertEquals(2, tree.searchNbMatchesInIndex("g"));
		assertEquals(2, tree.searchNbMatchesInIndex("c"));
		assertEquals(2, tree.searchNbMatchesInIndex("n"));
		
		// search string length > maxWordLength of index
		tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgcnatgcnatgcnatgcn", 5, true);
		assertEquals(3, tree.searchNbMatchesInIndex("atgcnatgcn"));
		assertEquals(2, tree.searchNbMatchesInIndex("atgcnatgcnatgcn"));
		assertEquals(1, tree.searchNbMatchesInIndex("atgcnatgcnatgcnatgcn"));
		
		// search string length < maxWordLength of index
		tree = new DNASuffixTrie(alphabet);
		try{
			tree.createIndex("atgcn", 30, true);
			assertEquals(1, tree.searchNbMatchesInIndex("atgcn"));
			fail("Illegal argument exception is expected here!");
		}
		catch(IllegalArgumentException e){
			;
		}
	}

	public void testSearchMatchPositionsInIndex() {
		ObjectArrayList alphabet = new ObjectArrayList();
		alphabet.add((char) 'A');
		alphabet.add((char) 'T');
		alphabet.add((char) 'G');
		alphabet.add((char) 'C');
		alphabet.add((char) Constants.REPETITIVE_ELEMENT_CHARACTER.toUpperCase().charAt(0));
		
		DNASuffixTrie tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgcn", 5, true);
		Integer[] result1 = {0};
		Integer[] realResult = tree.searchMatchPositionsInIndex("atgcn");
		for(int i=0; i < Math.max(result1.length, realResult.length); i++) assertEquals(result1[i], realResult[i]);

		Integer[] result2 = {1};
		realResult = tree.searchMatchPositionsInIndex("tgcn");
		for(int i=0; i < Math.max(result2.length, realResult.length); i++) assertEquals(result2[i], realResult[i]);
		
		Integer[] result3 = {2};
		realResult = tree.searchMatchPositionsInIndex("gcn");
		for(int i=0; i < Math.max(result3.length, realResult.length); i++) assertEquals(result3[i], realResult[i]);
		
		Integer[] result4 = {3};
		realResult = tree.searchMatchPositionsInIndex("cn");
		for(int i=0; i < Math.max(result4.length, realResult.length); i++) assertEquals(result4[i], realResult[i]);

		Integer[] result5 = {4};
		realResult = tree.searchMatchPositionsInIndex("n");
		for(int i=0; i < Math.max(result5.length, realResult.length); i++) assertEquals(result5[i], realResult[i]);

		Integer[] result6 = {1};
		realResult = tree.searchMatchPositionsInIndex("tgc");
		for(int i=0; i < Math.max(result6.length, realResult.length); i++) assertEquals(result6[i], realResult[i]);
		
		tree = new DNASuffixTrie(alphabet);
		tree.createIndex("atgcnatgcn", 10, true);
		Integer[] result7 = {0};
		realResult = tree.searchMatchPositionsInIndex("atgcnatgcn");
		for(int i=0; i < Math.max(result7.length, realResult.length); i++) assertEquals(result7[i], realResult[i]);
		
		Integer[] result8 = {0,5};
		realResult = tree.searchMatchPositionsInIndex("atgcn");
		for(int i=0; i < Math.max(result8.length, realResult.length); i++) assertEquals(result8[i], realResult[i]);
	}

}
