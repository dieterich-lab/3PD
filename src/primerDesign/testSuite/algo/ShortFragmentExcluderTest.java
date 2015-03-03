/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.algo.RestrictionFragmentSizeFilter;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.SequenceRegion;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleContigImpl;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class ShortFragmentExcluderTest extends TestCase {

	/**
	 * Test method for {@link primerDesign.algo.RestrictionFragmentSizeFilter#excludeSites(cern.colt.list.ObjectArrayList)}.
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 */
	public void testExcludeSites() throws IllegalAlphabetException, IllegalSymbolException {
		PrimerSearchParameters params = new PrimerSearchParameters();
		params.setMIN_RESTRICTION_FRAGMENT_LENGTH(1000);
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("Dummy", DNATools.createDNA("atgc"), 0, 0);
		SequenceRegion region = new SequenceRegion(new SimpleContigImpl("Dummy"),0,3000);
		
		RestrictionSite site1 = new RestrictionSite(0, enzyme, params);
		site1.setSequenceRegion(region);
		RestrictionSite site2 = new RestrictionSite(1001, enzyme, params);
		site2.setSequenceRegion(region);
		RestrictionSite site3 = new RestrictionSite(2000, enzyme, params);
		site3.setSequenceRegion(region);
		RestrictionSite site4 = new RestrictionSite(2999, enzyme, params);
		site4.setSequenceRegion(region);
		RestrictionSite site5 = new RestrictionSite(4000, enzyme, params);
		site5.setSequenceRegion(region);
		
		region.addRestrictionSite(site1);
		region.addRestrictionSite(site2);
		region.addRestrictionSite(site3);
		region.addRestrictionSite(site4);
		region.addRestrictionSite(site5);
		
		ObjectArrayList sites = new ObjectArrayList();
		sites.add(site1);
		sites.add(site2);
		sites.add(site3);
		sites.add(site4);
		sites.add(site5);
		
		sites = RestrictionFragmentSizeFilter.excludeSitesMinLength(sites, params.getMIN_RESTRICTION_FRAGMENT_LENGTH());
		
		assertTrue(sites.size() == 3);
		assertEquals((RestrictionSite) sites.get(0), site1);
		assertEquals((RestrictionSite) sites.get(1), site2);
		assertEquals((RestrictionSite) sites.get(2), site5);
	}
}
