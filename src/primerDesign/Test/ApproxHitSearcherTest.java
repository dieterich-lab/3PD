/**
 * 
 */
package primerDesign.Test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.dsc.indexStructures.esa.EnhancedSuffixArrayFatOpt;
import primerDesign.util.SeqTools;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class ApproxHitSearcherTest {

	/**
	 * Test method for {@link primerDesign.Test.ApproxHitSearcher#findHits(java.lang.String, primerDesign.dsc.indexStructures.DNASequenceIndex, int, int, int, int, int)}.
	 */
	@Test
	public void testFindHits() {
		DNASequenceIndex index = new EnhancedSuffixArrayFatOpt("TGTAATCACATGAAGAAGTACTTGG", "testname");
		index.createIndex();
		ObjectArrayList hits = ApproxHitSearcher.findHits("TGTAATCACATGAAGAAGTACTTGG", index, 11, 1, 5, 2, 0);
		ObjectArrayList expect = new ObjectArrayList();
		expect.add(new IndexHitImpl(index.getContig()[0], 0, true));
		assertEquals(expect, hits);
		
		hits = ApproxHitSearcher.findHits(SeqTools.complementDNA("TGTAATCACATGAAGAAGTACTTGG".toCharArray()), index, 11, 1, 5, 2, 0);
		expect = new ObjectArrayList();
		expect.add(new IndexHitImpl(index.getContig()[0], 0, false));
		assertEquals(expect, hits);
		
		String query = "TGTAATCACATGAAGAAGTACTTGG";
		String pad = "AAAAAAAAA";
		index = new EnhancedSuffixArrayFatOpt(query.substring(14) + pad + query + pad + SeqTools.complementDNA(query.toCharArray()) + pad + SeqTools.complementDNA(query.substring(14).toCharArray()), "testname");
		index.createIndex();
		hits = ApproxHitSearcher.findHits(query, index, 11, 1, 5, 2, 10);
		expect = new ObjectArrayList();
		expect.add(new IndexHitImpl(index.getContig()[0], 20, true));
		expect.add(new IndexHitImpl(index.getContig()[0], 54, false));
		assertEquals(expect, hits);
		
		index = new EnhancedSuffixArrayFatOpt(query + pad + SeqTools.complementDNA(query.toCharArray()), "testname");
		index.createIndex();
		hits = ApproxHitSearcher.findHits(query, index, 11, 1, 5, 2, 10);
		expect = new ObjectArrayList();
		expect.add(new IndexHitImpl(index.getContig()[0], 0, true));
		expect.add(new IndexHitImpl(index.getContig()[0], 34, false));
		assertEquals(expect, hits);
		
		index = new EnhancedSuffixArrayFatOpt(query + pad + query.substring(14), "testname");
		index.createIndex();
		hits = ApproxHitSearcher.findHits(query, index, 11, 1, 5, 2, 8);
		expect = new ObjectArrayList();
		expect.add(new IndexHitImpl(index.getContig()[0], 0, true));
		expect.add(new IndexHitImpl(index.getContig()[0], 20, true));
		assertEquals(expect, hits);
	}
}
