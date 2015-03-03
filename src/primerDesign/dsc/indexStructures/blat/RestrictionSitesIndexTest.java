/**
 * 
 */
package primerDesign.dsc.indexStructures.blat;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * Implements unit tests for a restrictino sites index.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RestrictionSitesIndexTest {

	/**
	 * Test method for {@link primerDesign.dsc.indexStructures.blat.RestrictionSitesIndex#addSites(java.lang.String, java.lang.String, java.lang.String, int[])}.
	 */
	@Test
	public void testAddSites() {
		
	}

	/**
	 * Test method for {@link primerDesign.dsc.indexStructures.blat.RestrictionSitesIndex#getSites(java.lang.String, java.lang.String, java.lang.String)}.
	 */
	@Test
	public void testGetSites() {
		
	}

	/**
	 * Test method for {@link primerDesign.dsc.indexStructures.blat.RestrictionSitesIndex#isSufficientlyClose(int, int, java.lang.String, java.lang.String, java.lang.String)}.
	 */
	@Test
	public void testIsSufficientlyClose() {
		RestrictionSitesIndex index = new RestrictionSitesIndex();
		index.addSites("Enzyme", "Organism", "Chromosome", new int[]{0, 1000, 10000});
		assertFalse(index.isSufficientlyClose(20000, 1000, "Enzyme", "Organism", "Chromosome"));
		assertTrue(index.isSufficientlyClose(2000, 1000, "Enzyme", "Organism", "Chromosome"));
		assertFalse(index.isSufficientlyClose(2001, 1000, "Enzyme", "Organism", "Chromosome"));
	}

}
