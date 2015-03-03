package primerDesign.testSuite.algo;

import java.util.Vector;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.algo.RestrictionSiteSearch;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.NoRestrictionSitesFoundException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleContigImpl;
import weka.core.FastVector;

public class RestrictionSiteSearchTest extends TestCase {

	public void testGetMostHomogenuousRSSs() throws IllegalAlphabetException, IllegalSymbolException {
		PrimerSearchParameters params = new PrimerSearchParameters();
		SimpleContig[] contigs;
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		
		// no site
		String sequence = "aaaaaaaa";

		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		RestrictionSiteSearch search = new RestrictionSiteSearch();
		try{
			assertEquals(new FastVector(), search.getMostHomogenuousRSSs(enzyme, 1, params));
		}
		catch(NoRestrictionSitesFoundException e){
			;
		}
		
		// uncentred sites in each interval
		// one site in whole sequence
		sequence = "aaagaattcaaa";
		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		//              |4
		
		search = new RestrictionSiteSearch();
		FastVector hits = search.getMostHomogenuousRSSs( enzyme, 1, params);
		RestrictionSite expectSite = new RestrictionSite(4, enzyme, params);
		assertEquals(expectSite.getPosition(), ((RestrictionSite) hits.elementAt(0)).getPosition());
		
		// one site in each interval
		//           int1       |    int2
		sequence = "aaagaattcaaaaaagaattcaaa";	
		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		//              |4          |16
		hits = search.getMostHomogenuousRSSs( enzyme, 2, params);
		RestrictionSite expectSite1 = new RestrictionSite(4, enzyme, params);
		RestrictionSite expectSite2 = new RestrictionSite(16, enzyme, params);
		assertEquals(expectSite1.getPosition(), ((RestrictionSite) hits.elementAt(0)).getPosition());
		assertEquals(expectSite2.getPosition(), ((RestrictionSite) hits.elementAt(1)).getPosition());
		
		// more than one site in each interval
		// create left site of sequence
		sequence = enzyme.getRecognitionSite().seqString();
		for(int i=0; i<params.getMIN_RESTRICTION_FRAGMENT_LENGTH(); i++) sequence += "A";
		sequence += enzyme.getRecognitionSite().seqString();
		for(int i=0; i<params.getMIN_RESTRICTION_FRAGMENT_LENGTH(); i++) sequence += "A";
		sequence += enzyme.getRecognitionSite().seqString();
		// insert separating sequence
		sequence += "TTT";
		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		
		// create right site of sequence
		sequence += enzyme.getRecognitionSite().seqString();
		for(int i=0; i<params.getMIN_RESTRICTION_FRAGMENT_LENGTH(); i++) sequence += "A";
		sequence += enzyme.getRecognitionSite().seqString();
		for(int i=0; i<params.getMIN_RESTRICTION_FRAGMENT_LENGTH(); i++) sequence += "A";
		sequence += enzyme.getRecognitionSite().seqString();
		
		//                          |16                                 |52
		hits = search.getMostHomogenuousRSSs( enzyme, 2, params);
		expectSite1 = new RestrictionSite(1007,enzyme, params);
		expectSite2 = new RestrictionSite(3028,enzyme, params);
		assertEquals(expectSite1.getPosition(), ((RestrictionSite) hits.elementAt(0)).getPosition());
		assertEquals(expectSite2.getPosition(), ((RestrictionSite) hits.elementAt(1)).getPosition());
		
		// sites centred w.r.t each interval
		// one site centred in whole sequence
		sequence = "aaaagaattc";
		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		hits = search.getMostHomogenuousRSSs( enzyme, 1, params);
		expectSite = new RestrictionSite(5,enzyme, params);
		assertEquals(expectSite.getPosition(), ((RestrictionSite) hits.elementAt(0)).getPosition());
		
		// one site centred in each interval
		      //      int1       |    int2
		sequence = "aaaagaattcaaaagaattc";
		contigs = new SimpleContigImpl[1];
		contigs[0] = new SimpleContigImpl("TEST", sequence.toCharArray());
		params.setContigs(contigs);
		hits = search.getMostHomogenuousRSSs(enzyme, 2, params);
		expectSite1 = new RestrictionSite(5, enzyme, params);
		expectSite2 = new RestrictionSite(15, enzyme, params);
		assertEquals(expectSite1.getPosition(), ((RestrictionSite) hits.elementAt(0)).getPosition());
		assertEquals(expectSite2.getPosition(), ((RestrictionSite) hits.elementAt(1)).getPosition());
	}
}
