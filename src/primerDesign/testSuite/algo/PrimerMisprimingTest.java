package primerDesign.testSuite.algo;

import junit.framework.TestCase;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.Test.EnhancedSuffixArrayFatOpt;
import primerDesign.algo.PrimerMispriming;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.util.PrimerSearchParameters;
import cern.colt.list.ObjectArrayList;

/**
 * Tests whether a primer mispriming has a safe distance to the next restriction site.
 * 
 * Since the background sequence is supposed to contain the primer scan sequence,
 * a mispriming in this context occurs iff a second hit of the primer end in the
 * background sequence is found and iff this hit is "sufficiently close" to the
 * next restriction site.
 * 
 * @author froehler
 *
 */
public class PrimerMisprimingTest extends TestCase {

	public void testHasSafeDistanceToNextRSS() throws IllegalAlphabetException, IllegalSymbolException {
		PrimerSearchParameters params = new PrimerSearchParameters();
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		params.setEnzyme(enzyme);
		String searchRegion = "TTTTTTTTCGGTTTGCGGGC";
		int intraSequence = params.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH() - searchRegion.length() - enzyme.getRecognitionSite().length();
			
		// "set distance of FP match to "safe distance +1""
		String bg = searchRegion;
		for(int i=0; i< intraSequence+1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		// insert false positive mispriming
		for(int i=0; i< intraSequence+1; i++) bg += "A";
		bg += searchRegion;		
		for(int i=0; i< intraSequence+1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		primerDesign.dsc.indexStructures.DNASequenceIndex bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		//params.setBackgroundIndex(bgIndex);
		
		ObjectArrayList hits = bgIndex.findHitPositions(searchRegion.substring(searchRegion.length() - params.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
		
		assertTrue(PrimerMispriming.hasSafeDistanceToNextRSS(new Primer(searchRegion,PrimerTypes.forwardPrimer, params), enzyme, hits, params));
		
		// "set distance of FP match to "safe distance""
		bg = searchRegion;
		for(int i=0; i< intraSequence; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		// insert false positive mispriming
		for(int i=0; i< intraSequence; i++) bg += "A";
		bg += searchRegion;	
		for(int i=0; i< intraSequence; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		//params.setBackgroundIndex(bgIndex);
		
		hits = bgIndex.findHitPositions(searchRegion.substring(searchRegion.length() - params.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
		
		assertFalse(PrimerMispriming.hasSafeDistanceToNextRSS(new Primer(searchRegion,PrimerTypes.forwardPrimer, params), enzyme, hits, params));
		
		// "set distance of FP match to "safe distance -1""
		bg = searchRegion;
		for(int i=0; i< intraSequence-1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		// insert false positive mispriming
		for(int i=0; i< intraSequence-1; i++) bg += "A";
		bg += searchRegion;	
		for(int i=0; i< intraSequence-1; i++) bg += "A";
		bg += enzyme.getRecognitionSite().seqString();
		
		bgIndex = new EnhancedSuffixArrayFatOpt(bg, "test");
		bgIndex.createIndex();
		//params.setBackgroundIndex(bgIndex);
		
		hits = bgIndex.findHitPositions(searchRegion.substring(searchRegion.length() - params.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
		
		assertFalse(PrimerMispriming.hasSafeDistanceToNextRSS(new Primer(searchRegion,PrimerTypes.forwardPrimer, params), enzyme, hits, params));
	}
}
