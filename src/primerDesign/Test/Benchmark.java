/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.util.Vector;
import java.util.regex.Pattern;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;

import primerDesign.algo.PrimerSearch;
import primerDesign.algo.SimpleGreedyPrimerPairPicking;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.indexStructures.TargetOrganisms;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheckDeserializer;
import primerDesign.util.Constants;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import primerDesign.web.PrimerDesignWebProperites;

/**
 * @author froehler
 *
 */
public class Benchmark {
	
	private static final int iter = 10;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// eval early vs late mismatch verification and greedy vs exhaustive screening and sorted vs unsorted search
		
		// EARLY MM check
		Constants.doEarlyMMScan = true;
		System.out.println("STATUS: EARLY MISMATCH SCAN");
			// GREEDY PICK
			Constants.doGreedySitesScreening = true;
			System.out.println("STATUS: GREEDY SITES SCREENING");
				// SORTED SCREEN
				Constants.doSortedDOPTScreen = true;
				System.out.println("STATUS: SORTED dOPT SCREEN");
				for(int i=0; i<iter; i++){
					doRun();
					Runtime.getRuntime().gc();
				}
				// UNSORTED SCREEN
				Constants.doSortedDOPTScreen = false;
				System.out.println("STATUS: UNSORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
			// EXHAUSTIVE PICK
			Constants.doGreedySitesScreening = false;
			System.out.println("STATUS: EXHAUSTIVE SITES SCREENING");
				// SORTED SCREEN
				Constants.doSortedDOPTScreen = true;
				System.out.println("STATUS: SORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
				// UNSORTED SCREEN
				Constants.doSortedDOPTScreen = false;
				System.out.println("STATUS: UNSORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
		// LATE MM check
		Constants.doEarlyMMScan = false;
		System.out.println("STATUS: LATE MISMATCH SCAN");
			// GREEDY PICK
			Constants.doGreedySitesScreening = true;
			System.out.println("STATUS: GREEDY SITES SCREENING");
				// SORTED SCREEN
				Constants.doSortedDOPTScreen = true;
				System.out.println("STATUS: SORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
				// UNSORTED SCREEN
				Constants.doSortedDOPTScreen = false;
				System.out.println("STATUS: UNSORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
			// EXHAUSTIVE PICK
			Constants.doGreedySitesScreening = false;
			System.out.println("STATUS: EXHAUSTIVE SITES SCREENING");
				// SORTED SCREEN
				Constants.doSortedDOPTScreen = true;
				System.out.println("STATUS: SORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
				// UNSORTED SCREEN
				Constants.doSortedDOPTScreen = false;
				System.out.println("STATUS: UNSORTED dOPT SCREEN");
				for(int i=0; i<iter; i++) doRun();
	}
	
	public static void doRun(){
		PrimerSearch search = new PrimerSearch();
		PrimerSearchParameters searchParameters = null;
		
		String scanRegionFilename = "/Users/froehler/Documents/Projects/HOX cluster analysis/experiments/000052/EcoRI-6pairs-SeqForMuPleX.result";
		String backgroundFilename = "/Users/froehler/SEQUENCE_INDICES/Ppacificus-unmasked.fa.MultiSeqESAFatOpt.esaidx";
		int nbPrimers = 6;
		String enzymeName = "EcoRI";
		String enzymeSite = "gaattc";
		search.setEnzymePatterns(Pattern.compile(".*" + enzymeSite.toUpperCase() + ".*"), Pattern.compile(".*" + SeqTools.revcompDNA(enzymeSite.toUpperCase().toCharArray()) + ".*"));
		int enzymeFwCutPos = 0;
		int enzymeRevCutPos = 0;
		
		try{
			SimpleTimer timer = new SimpleTimer();
			timer.startTimer();
			//System.out.println("Reading sequences:");
			//System.out.print("Reading scan region - ");
		
			// read first sequence from fasta
			SlimFastaParser parser = new SlimFastaParser(new File(scanRegionFilename));
			Vector<SimpleContigImpl> contigs = new Vector<SimpleContigImpl>();
			int counter = 0;
			while(parser.hasNextContig()){
				contigs.add(parser.parseNextContigIgnoreCase());
				counter++;
			}
			//System.out.println("done " + timer.getTimeString() + " - read " + counter + " sequences");
			
		
			//System.out.print("Constructing scanRegion index - ");

			
			//System.out.print("Constructing background index - ");
		
			
			RestrictionEnzyme enzyme = new RestrictionEnzyme(enzymeName, DNATools.createDNA(enzymeSite), enzymeFwCutPos, enzymeRevCutPos	);
			
			searchParameters = new PrimerSearchParameters(search);
			searchParameters.setNumPrimers(nbPrimers);
			searchParameters.setEnzyme(enzyme);
			searchParameters.setTargetOrganism(TargetOrganisms.Ppacificus);
			searchParameters.setContigs(contigs.toArray(new SimpleContigImpl[contigs.size()]));
			searchParameters.setConfigFile(PrimerDesignWebProperites.parseConfigFile(new FileInputStream("/Users/froehler/Desktop/Projects/EclipseWorkspace/PrimerDesign/deploy/3CPrimerDesign-ServletConfig.xml")));
			
			PrimerMisprimingCheck misprimingCheck = PrimerMisprimingCheckDeserializer.deserialize(backgroundFilename, searchParameters);
			
			searchParameters.setPrimerMisprimingCheck(misprimingCheck);
			
//			System.out.println("done " + timer.getTimeString());
//			System.out.println("Scanning sequence for valid primers");			
			search.setDoStat(searchParameters.isComputeScanningStatistics());
			timer.getTimeString();
			RestrictionSite[] optimalSites = search.naivePrimerSearch(searchParameters);
			System.out.print(timer.getTimeString() + "\t");						

			SimpleGreedyPrimerPairPicking picker = new SimpleGreedyPrimerPairPicking();
			PrimerPairSet bestPrimerPairSet = null; // = picker.pickBestPrimerSet(optimalSites, searchParameters);
			
			try{
				bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
			}
			catch(EmptyResultSetException e){
				//System.out.println("\n" + e.getMessage());
				System.out.println("-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!");
				//if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
				if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
				if(searchParameters.isComputePickingStatistics() && searchParameters.getPickingStat() != null) System.out.println(searchParameters.getPickingStat().printStat());
				//System.exit(1);
			}			
			
//			
			System.out.print(timer.getTimeString() + "\t");
			
//			} 
			if(bestPrimerPairSet == null){
				System.out.println("-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!");
			//	if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
//				if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
//				if(searchParameters.isComputePickingStatistics()) System.out.println(searchParameters.getPickingStat().printStat());
			}
			else{
//				System.out.println("Displaying best primer set");
//				System.out.println(String.format("%10s\t%s", "Type", Primer.toFormattedStringDescription()));
//				for(int i=0; i<bestPrimerSet.getNumberOfForwardPrimers(); i++){
//					Primer currentPrimer = bestPrimerSet.getForwardPrimer(i);
//					//System.out.println(currentPrimer.getRelativePosition() + "\t" + currentPrimer.getSequence() + "\t" + currentPrimer.getDistanceToRSS() + "\t" + currentPrimer.getGcContent() + "\t" + currentPrimer.getMeltingTemp() + "\t" + currentPrimer.getDistanceToOptimalPrimer() + "\t" + currentPrimer.getSelfAlignmentScore() + "\t" + currentPrimer.getSelfEndAlignmentScore());
//					System.out.println(currentPrimer.toString());
				//System.out.println(bestPrimerPairSet.toFormattedString());
				System.out.print(bestPrimerPairSet.getMaxPairAlignScore() + " " + bestPrimerPairSet.getMaxPairAlignEndScore() + "\t");
				//System.out.println("Displaying best primer set took " + timer.getTimeString());
				
				System.out.print(bestPrimerPairSet.getHashCode() + " \t");

			//	if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
//				if(searchParameters.isComputeScanningStatistics() && searchParameters.isComputePickingStatistics()) System.out.println();
//				if(searchParameters.isComputePickingStatistics()) System.out.println(searchParameters.getPickingStat().printStat());
			}
			
			System.out.println(timer.getTotalTimestring());
//			System.out.println("Evaluated " + PrimerSearch.primerCount + " candidate primers in total.");
//			System.out.println("Evaluated " + PrimerSearch.primerIndexCount + " candidate primers in index in total.");
			
			// write timings to logfile
//			BufferedWriter writer = new BufferedWriter(new FileWriter("timings.log", true)); 
//			writer.write(enzymeName + "\t" + nbPrimers + "\t" + timer.getTotalTimestring() + "\n");
//			writer.close();
		}
		catch(Exception e){
			e.printStackTrace();
		//	if(searchParameters.isComputeScanningStatistics()) System.out.println(search.stat.getStat());
		}
	}

}
