/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.algo.PrimerSearch;
import primerDesign.algo.SimpleGreedyPrimerPairPicking;
import primerDesign.dsc.DNASequenceIndex;
import primerDesign.dsc.EnhancedSuffixArrayFatOpt;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.indexStructures.rmi.MultiSeqESAFatOpt;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SlimFastaParser;

/**
 * @author froehler
 *
 */
public class ConstrainedPrimerSearch {
	public static void main(String[] args) throws IOException, NumberFormatException, IllegalAlphabetException, IllegalSymbolException{
		PrimerSearch search = new PrimerSearch();
		PrimerSearchParameters params = new PrimerSearchParameters();
		
		// parse scan sequence
		SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
		String sequence = new String(parser.parseNextContigIgnoreCase().getSequence());
		params.setSequence(sequence);
		
		// create scan region index
		DNASequenceIndex scanRegionIndex = new EnhancedSuffixArrayFatOpt("");
		scanRegionIndex.createIndex(sequence.toString(), sequence.length(), true);
		params.setScanRegionIndex(scanRegionIndex);
		
		// deserialize background index
		params.setBackgroundIndex(MultiSeqESAFatOpt.deserialize(new File(args[1])));
		
		// set num primer pairs
		params.setNumPrimers(Integer.parseInt(args[2]));
		
		// set enzyme
		RestrictionEnzyme enzyme = new RestrictionEnzyme(args[3], DNATools.createDNA(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]));
		
		// get optimal restriction sites
		RestrictionSite[] optimalSites = search.naivePrimerSearch(params);
		
		// post-process optimal sites
		
		
		// scan for valid pairs on constrained set
		SimpleGreedyPrimerPairPicking picker = new SimpleGreedyPrimerPairPicking();
		PrimerPairSet bestPrimerPairSet = null; 
		
		try{
			bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, params);
		}
		catch(EmptyResultSetException e){
			System.out.println("\n" + e.getMessage());
			System.out.println("-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!");
			System.exit(1);
		}
		
		// output best set
		System.out.println("Displaying best primer set");
		System.out.println("Type\t" + Primer.toStringDescription());
		System.out.println(bestPrimerPairSet.toString());
		System.out.println("worst primer-  primer PA/PEA score (incl probes): " + bestPrimerPairSet.getMaxPairAlignScore() + " " + bestPrimerPairSet.getMaxPairAlignEndScore());
	}	
}
