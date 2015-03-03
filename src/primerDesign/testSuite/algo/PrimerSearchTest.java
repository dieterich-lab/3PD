package primerDesign.testSuite.algo;

import java.util.regex.Pattern;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.Test.EnhancedSuffixArrayFatOpt;
import primerDesign.algo.PrimerSearch;
import primerDesign.dsc.DNASequenceIndex;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;

public class PrimerSearchTest extends TestCase {

	public void testEnumeratePrimers() throws IllegalAlphabetException, IllegalSymbolException {
		PrimerSearchParameters params = new PrimerSearchParameters();
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);

		// exclude primer with masked basepairs check
		String sequence = "AAAAAAAAAANAAAAAAAAAAAAAAAAAAA";
		
		DNASequenceIndex index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(sequence, sequence.length(), true);
		
		primerDesign.dsc.indexStructures.DNASequenceIndex bgIndex = new EnhancedSuffixArrayFatOpt(sequence, "test");
		bgIndex.createIndex();
		
		PrimerSearch search = new PrimerSearch();
		search.setEnzymePatterns(Pattern.compile(".*" + enzyme.getRecognitionSite().seqString().toUpperCase() + ".*"), Pattern.compile(".*" + SeqTools.revcompDNA(enzyme.getRecognitionSite().seqString().toUpperCase().toCharArray()) + ".*"));
		params.setEnzyme(enzyme);
		
		RestrictionSite site = new RestrictionSite(params);
		
		Primer[] result = search.enumeratePrimers(sequence, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		
		// invalid primer check - according to constraints specified in Constants! (e.g.: extreme %GC, extreme T_m,...)
		// GC min
		sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		
		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(sequence, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(sequence, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		search = new PrimerSearch();
		result = search.enumeratePrimers(sequence, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		// GC max
		sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(sequence, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(sequence, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		search = new PrimerSearch();
		result = search.enumeratePrimers(sequence, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		// sa, sea max
		sequence = "GCGCGCGCGCGCGCGGCGCGCGCGCGCGCG";

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(sequence, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(sequence, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		search = new PrimerSearch();
		result = search.enumeratePrimers(sequence, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		
		// primer mispriming checks
		// another hit in search region
		String searchRegion = "ATATATATATATATAGCGCGCGCGCGCGCGATATATATATATATAGCGCGCGCGCGCGCG";
		String background = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(background, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		
		// another hit outside of search region
		searchRegion = "ATATATATATATATAGCGCGCGCGCGCGCGATATATATATATATAGCGCGCGCGCGCGCG";
		background = "ATATATATATATATAGCGCGCGCGCGCGCG";

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(background, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
		
		// NO other hit, neither within, nor outside of scan region
		searchRegion = "TTTTTTTTCGGTTTGCGGGCTTTTTTTTTTT";
		background = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(sequence);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(background, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(6, result.length);
		
		// fp amplicon checks
		int intraSequence = params.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH() - searchRegion.length();// - params.getMIN_AMPLICON_LENGTH()/2 - searchRegion.length();
		searchRegion = "TTTTCGGTTTGCGGGCTT";
		// fp amplicon > threshold
		String bg = searchRegion;
		// "set distance of FP match to "one above safe distance""
		for(int i=0; i< intraSequence+1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(searchRegion);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(1, result.length);
		
		// fp amplicon == threshold -> reject primer since constraint: distance > threshold!
		bg = searchRegion;
		// "set distance of FP match to "safe distance""
		for(int i=0; i< intraSequence; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(searchRegion);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(1, result.length);
		
		// fp amplicon < threshold
		bg = searchRegion;
		// "set distance of FP match to "one below safe distance""
		for(int i=0; i< intraSequence-1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		// insert false positive hit into background sequence
		bg += "AAAA" + searchRegion;
		for(int i=0; i< intraSequence-1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();

		index = new primerDesign.dsc.EnhancedSuffixArrayFatOpt(searchRegion);
		index.createIndex(searchRegion, sequence.length(), true);
		
		bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		
		site = new RestrictionSite(params);
		
		result = search.enumeratePrimers(searchRegion, PrimerTypes.forwardPrimer, params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH(), 0 , true, 45, site, enzyme, params);
		assertEquals(0, result.length);
	}

//	public void testNaivePrimerSearch() {
	// RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
//		fail("Not yet implemented");
//	}
//
//	public void testPickBestPrimerSet() {
	// RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
//		fail("Not yet implemented");
//	}

}
