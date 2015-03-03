package primerDesign.Test;

import junit.framework.TestCase;
import cern.colt.list.ObjectArrayList;

public class EnhancedSuffixArrayTest extends TestCase {

	public void testCreateIndex() {
		EnhancedSuffixArray array = new EnhancedSuffixArray();
		String sequence = "acaaacatat";
		array.createIndex(sequence);
		int[] suftab = new int[sequence.length()+1];
		int[] lcptab = new int[sequence.length()+1];
		int[] childtab = new int[sequence.length()+1];
		
		suftab[0] = 2;
		suftab[1] = 3;
		suftab[2] = 0;
		suftab[3] = 4;
		suftab[4] = 8;
		suftab[5] = 6;
		suftab[6] = 1;
		suftab[7] = 5;
		suftab[8] = 9;
		suftab[9] = 7;
		suftab[10] = Integer.MAX_VALUE;
		
		lcptab[0] = 0;
		lcptab[1] = 2;
		lcptab[2] = 1;
		lcptab[3] = 3;
		lcptab[4] = 1;
		lcptab[5] = 2;
		lcptab[6] = 0;
		lcptab[7] = 2;
		lcptab[8] = 0;
		lcptab[9] = 1;
		lcptab[10] = 0;
		
		childtab[0] = 6;
		childtab[1] = 1; 
		childtab[2] = 4;
		childtab[3] = 3;
		childtab[4] = 5;
		childtab[5] = 2;
		childtab[6] = 8;
		childtab[7] = 7;
		childtab[8] = 10;
		childtab[9] = 9;
		childtab[10] = 0;
		for(int i=0; i<suftab.length; i++){
			for(int j=0; j<4; j++){
				assertEquals(suftab[i], array.getSufTab(i));
				assertEquals(lcptab[i], array.getLcpTab(i));
				assertEquals(childtab[i], array.getChildTab(i));
			}
		}
	}

	public void testGetChildIntervals() {
		ObjectArrayList expected = new ObjectArrayList();
		expected.add(new Integer[]{0,5});
		expected.add(new Integer[]{6,7});
		expected.add(new Integer[]{8,9});
		expected.add(new Integer[]{10,10});
		
		EnhancedSuffixArray array = new EnhancedSuffixArray();
		String sequence = "acaaacatat";
		array.createIndex(sequence);
		ObjectArrayList intervals = array.getChildIntervals(0, sequence.length());
		ObjectArrayList actual = array.getChildIntervalsValues(intervals);
		for(int i=0; i<actual.size(); i++){
			assertEquals((Integer) (((Integer[]) (expected.get(i)))[0]), ((Integer[]) actual.get(i))[0]);
			assertEquals((Integer) (((Integer[]) (expected.get(i)))[1]), ((Integer[]) actual.get(i))[1]);
		}
	}

	public void testGetChildInterval() {
		EnhancedSuffixArray array = new EnhancedSuffixArray();
		String sequence = "acaaacatat";
		array.createIndex(sequence);
		int[] actual = array.intervalToArray(array.getChildInterval(0, sequence.length(), 'a'));
		int[] expect = new int[]{0,5};
		for(int i=0; i<actual.length; i++){
			assertEquals(expect[i], actual[i]);
		}
	}

	public void testFindMatchPositions() {
		EnhancedSuffixArray array = new EnhancedSuffixArray();
		String sequence = "acaaacatat";
		array.createIndex(sequence);
		int[] actual = array.findMatchPositions("aca");
		int[] expect = new int[]{0,4};
		for(int i=0; i<actual.length; i++){
			assertEquals(expect[i], actual[i]);
		}
		
		array = new EnhancedSuffixArray();
		 sequence = "acaaacatat";
		array.createIndex(sequence);
		actual = array.findMatchPositions("acat");
		expect = new int[]{4};
		for(int i=0; i<actual.length; i++){
			assertEquals(expect[i], actual[i]);
		}
	}
}
