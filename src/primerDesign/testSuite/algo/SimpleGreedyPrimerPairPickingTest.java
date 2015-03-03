/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.junit.Test;

import primerDesign.algo.PrimerSearch;
import primerDesign.algo.SimpleGreedyPrimerPairPicking;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.SequenceRegion;
import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.DummyIndex;
import primerDesign.util.PrimerSearchParameters;

/**
 * @author froehler
 *
 */
public class SimpleGreedyPrimerPairPickingTest extends TestCase{

	/**
	 * Test method for {@link primerDesign.algo.SimpleGreedyPrimerPairPicking#pickBestPrimerSet(primerDesign.dsc.RestrictionSite[], primerDesign.util.PrimerSearchParameters)}.
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 */
	@Test
	public void testPickBestPrimerSet() throws IllegalAlphabetException, IllegalSymbolException {
		// init method
		PrimerSearchParameters searchParams = new PrimerSearchParameters();
		RestrictionEnzyme enzyme = new RestrictionEnzyme("Dummy", DNATools.createDNA("CCCC"), 0, 0);
		searchParams.setEnzyme(enzyme);
		searchParams.setNumPrimers(2);
		
		RestrictionSite[] optimalSites = new RestrictionSite[2];
		
		RestrictionSite site1 = new RestrictionSite(5000, enzyme, searchParams);
		RestrictionSite site2 = new RestrictionSite(15000, enzyme, searchParams);
		
		SequenceRegion region1 = new SequenceRegion(null, 0,9999);
		SequenceRegion region2 = new SequenceRegion(null, 10000, 20000);
		
		region1.addRestrictionSite(site1);
		region2.addRestrictionSite(site2);
		
		optimalSites[0] = site1;
		optimalSites[1] = site2;
		
		site1.setDistanceToIntervalMean(0);
		site2.setDistanceToIntervalMean(0);
		
		site1.setForwardScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGATTTTTTTTTTT".toCharArray());
		site1.setReverseScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGATTTTTTTTTTT".toCharArray());
		site1.setProbeScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGATTTTTTTTT".toCharArray());
		
		site2.setForwardScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGATTTTTTTTTTT".toCharArray());
		site2.setReverseScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGATTTTTTTTTTT".toCharArray());
		site2.setProbeScanSequence("TTTTTTTTTTAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGATTTTTTTTT".toCharArray());
		
		PrimerSearch search = new PrimerSearch();
		searchParams.setPrimerSearch(search);
		
		primerDesign.dsc.DNASequenceIndex scanRegionIndex = new DummyIndex();
		DNASequenceIndex backgroundIndex = new DummyIndex();
		//searchParams.setScanRegionIndex(scanRegionIndex);
		//searchParams.setBackgroundIndex(backgroundIndex);
		
		site1.setValidUpstreamPrimers(search.enumeratePrimers(new String(site1.getForwardScanSequence()), PrimerTypes.forwardPrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), 4900, true, 100, site1, enzyme, searchParams));
		site1.setValidDownstreamPrimers(search.enumeratePrimers(new String(site1.getReverseScanSequence()), PrimerTypes.reversePrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), 5100, true, 100, site1, enzyme, searchParams));
		site1.setValidTaqManProbes(search.enumeratePrimers(new String(site1.getProbeScanSequence()), PrimerTypes.hybridizationProbe, searchParams.getTAQMAN_MIN_PRIMER_LENGTH(), searchParams.getTAQMAN_MAX_PRIMER_LENGTH(), 5000, true, 0, site1, enzyme, searchParams));
		
		site2.setValidUpstreamPrimers(search.enumeratePrimers(new String(site2.getForwardScanSequence()), PrimerTypes.forwardPrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), 14900, true, 100, site2, enzyme, searchParams));
		site2.setValidDownstreamPrimers(search.enumeratePrimers(new String(site2.getReverseScanSequence()), PrimerTypes.reversePrimer, searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH(), 15100, true, 100, site2, enzyme, searchParams));
		site2.setValidTaqManProbes(search.enumeratePrimers(new String(site2.getProbeScanSequence()), PrimerTypes.hybridizationProbe, searchParams.getTAQMAN_MIN_PRIMER_LENGTH(), searchParams.getTAQMAN_MAX_PRIMER_LENGTH(), 15000, true, 0, site2, enzyme, searchParams));

		Primer upstream = new Primer("AGAGAGAGAGAGAGAGAGAGA", PrimerTypes.forwardPrimer, 10, searchParams);
		Primer downstream = new Primer("AGAGAGAGAGAGAGAGAGAGA", PrimerTypes.reversePrimer, 10, searchParams);
		Primer probe = new Primer("AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA", PrimerTypes.hybridizationProbe, 10, searchParams);
		
		assertEquals(upstream.getSequence(), site1.getValidUpstreamPrimers()[0].getSequence());
		assertEquals(downstream.getSequence(), site1.getValidDownstreamPrimers()[0].getSequence());
		assertEquals(probe.getSequence(), site1.getValidTaqManProbes()[0].getSequence());
		
		assertEquals(upstream.getSequence(), site2.getValidUpstreamPrimers()[0].getSequence());
		assertEquals(downstream.getSequence(), site2.getValidDownstreamPrimers()[0].getSequence());
		assertEquals(probe.getSequence(), site2.getValidTaqManProbes()[0].getSequence());
		
		PrimerPair targetPair1 = new PrimerPair(upstream, downstream, probe, searchParams);
		PrimerPair targetPair2 = new PrimerPair(upstream, downstream, probe, searchParams);
		
		site1.enumeratePrimerPairs();
		site2.enumeratePrimerPairs();
		
		// pick best primer pairs
		System.out.println("Picking from: " + optimalSites[0].getNumberOfValidPrimerPairs() + " x " + optimalSites[1].getNumberOfValidPrimerPairs() + " Primer pair sets");
		SimpleGreedyPrimerPairPicking picker = new SimpleGreedyPrimerPairPicking();
		PrimerPairSet optimalSet = picker.pickBestPrimerSet(optimalSites, searchParams);
		
		assertEquals(targetPair1.getForwardPrimer().getSequence(), optimalSet.getPrimerPair(0).getForwardPrimer().getSequence());
		
		assertEquals(targetPair2.getForwardPrimer().getSequence(), optimalSet.getPrimerPair(1).getForwardPrimer().getSequence());
		
		System.out.println(optimalSet.toString());
	}
}
