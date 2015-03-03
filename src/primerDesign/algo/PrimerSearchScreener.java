package primerDesign.algo;

import java.io.File;
import java.util.ArrayList;
import java.util.regex.Pattern;

import org.biojava.bio.molbio.RestrictionEnzyme;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerPairPickingStatistics;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.indexStructures.TargetOrganisms;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheckDeserializer;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.RestrictionEnzymeListParser;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;

public class PrimerSearchScreener {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 4){
			System.out.println("Usage: PrimerSearchScreener <Target Region Sequence> <Background Region Index> <# Primer Pairs> <Enzymes File>");
			System.out.println();
			System.exit(1);
		}
		PrimerSearch search = new PrimerSearch();
		
		String scanRegionFilename = args[0];
		String backgroundFilename = args[1];
		int nbPrimers = Integer.parseInt(args[2]);
		String enzymesFilename = args[3];
		Enum targetOrganism = TargetOrganisms.Ppacificus;

		try{
			SimpleTimer timer = new SimpleTimer();
			timer.startTimer();
			
			System.err.println("Reading sequences:");
			System.err.print("Reading scan region - ");			
			String sequence = null;
			SlimFastaParser parser = new SlimFastaParser(new File(scanRegionFilename));
			if(parser.hasNextContig()) sequence = new String(parser.parseNextContigIgnoreCase().getSequence());			
			System.err.println("done " + timer.getTimeString());
			
//			System.err.print("Constructing scanRegion index - ");
//			DNASequenceIndex scanRegionIndex = new EnhancedSuffixArrayFatOpt("");
//			scanRegionIndex.createIndex(sequence, sequence.length(), true);
//			System.err.println("done " + timer.getTimeString());
			
//			System.err.print("Constructing background index - ");
//			System.err.print("deserializing background sequence index - ");
//			MultiSeqESAFatOpt backgroundIndex = MultiSeqESAFatOpt.deserialize(new File(backgroundFilename));
//			System.err.println("done " + timer.getTimeString());
			
			ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File(enzymesFilename));
			RestrictionEnzyme currentEnzyme;
			PrimerPairSet finalBestPrimerPairSet = null;
			int enzymeIndex = -1;
			
			double currentBestHomogenityScore = Double.MAX_VALUE;
			double currentdOptPrimerPairSet = Double.MAX_VALUE;
			int pa_max = Integer.MAX_VALUE;
			int pea_max = Integer.MAX_VALUE;
			
			PrimerPairPickingStatistics stat = new PrimerPairPickingStatistics();
			SimpleGreedyPrimerPairPicking picker = null;
			
			for(int k=0; k<enzymes.size(); k++){
				currentEnzyme = enzymes.get(k);
				PrimerPairSet bestPrimerPairSet;
				PrimerSearchParameters searchParameters;

				try{
					
					search.setEnzymePatterns(Pattern.compile(".*" + currentEnzyme.getRecognitionSite().seqString().toUpperCase() + ".*"), Pattern.compile(".*" + SeqTools.revcompDNA(currentEnzyme.getRecognitionSite().seqString().toUpperCase().toCharArray()) + ".*"));
					
					System.err.println("Scanning sequence for valid primers using enzyme: " + currentEnzyme.getName());
					
					searchParameters = new PrimerSearchParameters(search);
					searchParameters.setContigs(new SimpleContigImpl[]{new SimpleContigImpl("TEST", sequence.toCharArray())});
					searchParameters.setNumPrimers(nbPrimers);
					searchParameters.setEnzyme(currentEnzyme);
					searchParameters.setTargetOrganism(targetOrganism);
//					searchParameters.setScanRegionIndex(scanRegionIndex);
//					searchParameters.setBackgroundIndex(backgroundIndex);
					
					PrimerMisprimingCheck misprimingCheck = PrimerMisprimingCheckDeserializer.deserialize(backgroundFilename, searchParameters);
					searchParameters.setPrimerMisprimingCheck(misprimingCheck);
					
					RestrictionSite[] optimalSites = search.naivePrimerSearch(searchParameters);
					System.err.println("Scanning sequence for valid primers took " + timer.getTimeString());
											
					System.err.println("Computing best primer set");
					//PrimerPickingAlgorithm picker = new SimpleGreedyPrimerPicking();
					//PrimerSet bestPrimerSet = picker.pickBestPrimerSet(optimalSites);
					
					picker = new SimpleGreedyPrimerPairPicking();
					bestPrimerPairSet = null;
				
					
					bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
					stat.combineStats(picker.getStat());
				}
				catch(Exception e){
					System.err.println("Exception thrown for restiction enzyme: " + enzymes.get(k) + "\n" + e.toString());
					if(picker != null) stat.combineStats(picker.getStat());
					continue;
				}
				
	//			AdvacedGreedyPrimerPairPicking picker = new AdvacedGreedyPrimerPairPicking();
	//			PrimerPairSet bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParameters);
			
			
				System.err.println("Computing best primer set took " + timer.getTimeString());
				
				if(finalBestPrimerPairSet != null && bestPrimerPairSet != null){
					System.err.println("Current primer pair set to be evaluated:\n" + bestPrimerPairSet.toString() + "\n");
				}	
				if(bestPrimerPairSet != null && bestPrimerPairSet.getHomogenityScore() * searchParameters.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + bestPrimerPairSet.getAvgDistOptPrimerPair() * searchParameters.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)bestPrimerPairSet.getMaxPairAlignScore())/searchParameters.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParameters.getPAIR_ALIGNMENT_WEIGHT() + ((double)bestPrimerPairSet.getMaxPairAlignEndScore())/searchParameters.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParameters.getPAIR_END_ALIGNMENT_WEIGHT() < currentBestHomogenityScore * searchParameters.getPRIMER_PAIR_HOMOGENITY_WEIGHT() + currentdOptPrimerPairSet * searchParameters.getPRIMER_PAIR_DOPT_WEIGHT() + ((double)pa_max)/searchParameters.getMAX_PRIMER_PAIR_ALIGNMENT_SCORE() * searchParameters.getPAIR_ALIGNMENT_WEIGHT() + ((double)pea_max)/searchParameters.getMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE() * searchParameters.getPAIR_END_ALIGNMENT_WEIGHT()){
						finalBestPrimerPairSet = bestPrimerPairSet;
						currentBestHomogenityScore = bestPrimerPairSet.getHomogenityScore();
						currentdOptPrimerPairSet = bestPrimerPairSet.getAvgDistOptPrimerPair();
						pa_max = bestPrimerPairSet.getMaxPairAlignScore();
						pea_max = bestPrimerPairSet.getMaxPairAlignEndScore();
						enzymeIndex = k;
				}
				if(finalBestPrimerPairSet != null){
					System.err.println("Current best set:");
					System.err.println("Enzyme: " + enzymes.get(enzymeIndex));
					System.err.println(finalBestPrimerPairSet.toString());
					System.err.println();
				}
			}
			
			System.out.println("Screened setups for " + enzymes.size() + " restriction enzymes.");
			if(finalBestPrimerPairSet == null){
				System.out.println("-- No valid primer pair set can be picked at all for any enzyme - check your selection parameters/ modify enzyme list!\n");
			}
			else{
				System.out.println("The following setup was found to be best:");
				System.out.println("Enzyme: " + enzymes.get(enzymeIndex));
				System.out.println("Displaying best primer set");
				System.out.println("Type\t" + Primer.toStringDescription());
				System.out.println(finalBestPrimerPairSet.toString());
				System.out.println("Displaying best primer set took " + timer.getTimeString());
			}
			System.out.println("Total runtime: " + timer.getTotalTimestring());
			System.out.println("Evaluated " + search.getPrimerCount() + " candidate primers in total.");
			System.out.println("Evaluated " + search.getPrimerIndexCount() + " candidate primers in index in total.");		
			//System.out.println(stat.printStat());
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
}
