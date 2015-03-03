package primerDesign.testSuite.dsc;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.SequenceRegion;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleContigImpl;

public class SequenceRegionTest extends TestCase {

	public void testGetRestrictionSitesIterator() throws IllegalAlphabetException, IllegalSymbolException {
		// init
		SequenceRegion region = new SequenceRegion(new SimpleContigImpl("TEST"),0,100);
		assertFalse(region.getRestrictionSitesIterator(false).hasNext());
		
		//add sites, check if sites are returned ordered w.r.t sequence region mean
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		PrimerSearchParameters params = new PrimerSearchParameters();
		RestrictionSite site1 = new RestrictionSite(10, enzyme, params);
		site1.setDistanceToIntervalMean(40);
		region.addRestrictionSite(site1);
		
		RestrictionSite site2 = new RestrictionSite(45, enzyme, params);
		site2.setDistanceToIntervalMean(5);
		region.addRestrictionSite(site2);
		
		RestrictionSite site3 = new RestrictionSite(90, enzyme, params);
		site3.setDistanceToIntervalMean(40);
		region.addRestrictionSite(site3);
		
		// check auto-sort upon iterator call iterator
		assertEquals(site2, region.getRestrictionSitesIterator(false).next());
		
		// check explicit sort call
		region.sortRestrictionSites();
		assertEquals(site2, region.getRestrictionSitesIterator(false).next());
	}
}
