package primerDesign.algo;

import java.util.Iterator;
import java.util.TreeSet;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.indexStructures.IndexHit;
import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import cern.colt.list.ObjectArrayList;

/**
 * Encapsulates tools to scan for whether a mispriming of a primer affects the 3C-qPCR experiment.
 * 
 * A mispriming of a primer is defined to affect a 3C-qPCR experiment iff the distance to the next
 * restriction site of the specified restriction enzyme is below a used-defined threshold.
 * In such a case, a false positive amplicon could be generated which can not easily discriminated
 * from the desired target amplicon.
 * Only if the FP amplicon is above a specific length, it can safely be discriminated from the real qPCR amplicon!
 * This means that the next restriction site is more distant than this threshold w.r.t the mispriming position.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerMispriming {
	private static RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
	
	/**
	 * Computes whether a false positive match can produce an amplicon which is capable of disturbing qPCR.
	 * 
	 * All false positive amplicons of length > Constants.SAFE_FALSE_POSITIVE_AMPLICON_LENGTH are assumed to
	 * have no influence on the qPCR experiment
	 * 
	 * @param primer the primer to check
	 * @param enzyme the restriction enzyme
	 * @param hits a list of hits of primer primer in the background sequence
	 * @param searchParams the primer search parameters
	 * 
	 * @return true iff length of FP amplicon is > Threshold (s.o.)
	 */
	public static boolean hasSafeDistanceToNextRSS(Primer primer, RestrictionEnzyme enzyme, ObjectArrayList hits, PrimerSearchParameters searchParams){
		int offset = searchParams.getSAFE_FALSE_POSITIVE_AMPLICON_LENGTH();
		//Integer[] positions = backgroundIndex.searchMatchPositionsInIndex(seq);
		
//		PrimerMispriming.positions = new ObjectArrayList();
//		if(!primerType.equals(PrimerTypes.hybridizationProbe)){
//			// hybridization probes are NOT tested for misprimings since flourescence is only emmited when those probes are cleaved by a polymerase! -> distancxe to next RSS is NOT important
//			Integer[] matches = backgroundIndex.searchMatchPositionsInIndex(sequence);
//			if(matches.length > Constants.MAX_PRIMER_MISPRIMING_CUTOFF){
//				for(int i=0; i<matches.length; i++) PrimerMispriming.positions.add(matches[i]);
//	//			int stop = backgroundIndex.getSequence().length();
//	//			int match = -1;
//	//			//Integer[] matches = backgroundIndex.searchMatchPositionsInIndex(sequence.toUpperCase());
//	//			//if(matches.length != 0){
//	//			//	Arrays.sort(matches);
//	//			//	for(int i=0; i<matches.length; i++) PrimerMispriming.positions.add(matches[i]);
//	//			//}
//	//			while(match<stop){
//	//				match = backgroundIndex.getSequence().indexOf(sequence.toUpperCase(), match +1);
//	//				// if primer has no mispriming in background seq
//	//				if(match == -1) break;
//	//				// else store position
//	//				else PrimerMispriming.positions.add(match);
//	//			}
//			}
//		}
		ObjectArrayList positions;
		//RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
		
		if(primer.getPrimerType().equals(PrimerTypes.hybridizationProbe)) return true;
		else{
			positions = hits;
		}
		
		synchronized(mapper){
			//mapper.addEnzyme(enzyme);
			RestrictionEnzymeManager.register(enzyme, new TreeSet());
			mapper.clearEnzymes();
			mapper.addEnzyme(enzyme);
		}
		
		// if no mispriming can be found XOR one 'mispriming' (original primer match!) can be found and the index includes the primer scan region 
		//if(PrimerMispriming.positions.size() == Constants.MAX_PRIMER_MISPRIMING_CUTOFF && !backgroundIndex.includesScanRegion()) return true;
		//else if(PrimerMispriming.positions.size() == Constants.MAX_PRIMER_MISPRIMING_CUTOFF +1 && backgroundIndex.includesScanRegion()) return true; 
		
		// a background sequence is supposed to contain the target genomic region. Therefore one 'mispriming' (the required priming) is expected!
		//if(PrimerMispriming.positions.size() <= searchParams.getMAX_PRIMER_MISPRIMING_CUTOFF() + 1) return true;
		
		// else scan each mispriming of forward and reverse primers
		boolean result = true;
		IndexHitImpl mispriming;
		// allow for one 'mispriming' (the desired priming position)
		// 2Do: check whether the first position which is skipped is indeed the real, desired priming position!!!
		
		// whether the one real and true positive binding has already been found
		boolean truePositivePrimingFound = false;
		String commonChars;
		char[] primerSeq = primer.getSequence().toCharArray();
		double tm;
		int start;
		int end;
		String scanRegion;
		int position;
		for(int i=0; i< positions.size(); i++){
			//truePositivePrimingFound = false;
			
			// the start position of a mispriming is the last 'PRIMER_END_MISMATCH_SCAN_LENGTH' basepairs of a primer
			mispriming = (IndexHitImpl) positions.getQuick(i);
			position = mispriming.getPosition();
			
			// check whether primer mispriming is sufficiently strong to promote a false positive amplification
			// (if the melting temperature of a mispriming is 'sufficienly' lower than the original TM, the mispriming should not take place in reality)
//			if(mispriming.isForwardHit()){
//				start = position;
//				end = Math.min(mispriming.getContigLength() - 1, position + primer.getLength() - 1);
//			}
//			else{
//				start = position;
//				end = Math.min(mispriming.getContigLength() - 1, position + primer.getLength() - 1);
//			}
			
			// this positioning requires hit positions to be left-based for forward as well as for reverse hits!!!
			// the position of the forward hit corresponds to the start of the last 'PRIMER_END_MISMATCH_SCAN_LENGTH()' bp of the primer -> needs to be shifted to extract whole primer sequence
			// the position of the reverse hit correpsonds to the end of the reverse primer and need NOT to be shifted!
			if(mispriming.isForwardHit()){
				// oligonucleotide primers are single stranded, therefore, a mispriming only affects the sequence at the 3' end of the mismatch
				start = Math.max(0, position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH() - primer.getLength());
				end = position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH() - 1;
				commonChars = SeqTools.getCommonChars(primerSeq, mispriming.getSubSequence(start, end));
			}
			else{
				// old, erroneous whole-hit sequence extraction
				start = Math.max(0, position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH() - primer.getLength());
				end = Math.min(position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH() - 1, mispriming.getContigLength() - 1);
				commonChars = SeqTools.getCommonChars(primerSeq, SeqTools.revcompDNA(mispriming.getSubSequence(start, end)).toCharArray());
				
				// new, correct whole-hit sequence extraction
//				start = position;
//				end = Math.min(mispriming.getContigLength() - 1, position + primer.getSequenceLength() - 1);
//				commonChars = SeqTools.getCommonChars(primerSeq, mispriming.getSubSequence(start, end));
			}
//			start = position;
//			end = Math.min(mispriming.getContigLength() - 1, position + primer.getLength() - 1);
			// computes the common chars between two sequences, the subsequence of a mispriming is automatically reverse-complemented in case it is a reverse hit!
			//commonChars = SeqTools.getCommonChars(primerSeq, mispriming.getSubSequence(start, end));
			
			if(commonChars.length() >= 2){
				tm = Constants.PRIMER_TM_CALC_METHOD.computeMonoAndDivalentCationCorrectedTM(commonChars, searchParams.getPRIMER_CONCENTRATION(), searchParams.getMONOVALENT_CATION_CONCENTRATION(), searchParams.getDIVALENT_CATION_CONCENTRATION(), searchParams.getDNTP_CONCENTRATION());
				// if the mispriming temperature differs significantly from the real priming temperature, do NOT further verify this mispriming but consider it as being NOT relevant!
				if(Math.abs(primer.getMeltingTemp() - tm) >= searchParams.getMIN_MISPRIMING_TM_DIFFERENCE()){
					if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Skip ckeck: Mispriming TM diff > threshold: " + tm + "/(" + primer.getMeltingTemp() + ") lengths: " + commonChars.length() + "(" + primer.getLength() + ")");
					continue;
				}
				if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Verify: Mispriming TM diff < threshold: " + tm + "/(" + primer.getMeltingTemp() + ") lengths: " + commonChars.length() + "(" + primer.getLength() + ") " + primer.getSequence());
			}
			else{
				if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Mispriming of length: " + commonChars.length() + " - SKIP");
				continue;
			}
			
			// oligonucleotide primers are single-stranded, therefore, only the 3' part of a mispriming needs to be scanned for another RSS
			if(mispriming.isForwardHit()){
				start = Math.max(0, position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH() - primer.getLength());
				end = Math.min(start + offset, mispriming.getContigLength());
				scanRegion = mispriming.getContig().getSubsequence(start, end);
			}
			else{
				end = Math.min(position + searchParams.getPRIMER_END_MISMATCH_SCAN_LENGTH(), mispriming.getContigLength()); 
				start = Math.max(0, end - offset);
				scanRegion = mispriming.getContig().getSubsequence(start, end);
			}
			//int position = (Integer) PrimerMispriming.positions.get(i);
			//start = Math.max(0, position - offset);
			//int end = Math.min(position + offset, backgroundIndex.getSequence().length());
			//end = Math.min(position + offset, mispriming.getContigLength() - 1);
//			int length;
//			if(start == 0) length = Math.min(backgroundIndex.getSequence().length(), offset +1); 
//			else length = Math.min(backgroundIndex.getSequence().length(), 2*offset +1);
			//String scanRegion;
			//scanRegion = backgroundIndex.getSequence().substring(start, end);
			
			//String scanRegion = mispriming.getContig().getSubsequence(start, end);
			
//			// Only restriction sites in the direction of amplification can lead to false positive amplicons!!!
//			String scanRegion;
//			if(primer.getPrimerType().equals(PrimerTypes.forwardPrimer)){
//				scanRegion = mispriming.getContig().getSubsequence(position, end);
//			}
//			else if(primer.getPrimerType().equals(PrimerTypes.reversePrimer)){
//				scanRegion = mispriming.getContig().getSubsequence(start, position);
//			}
//			else throw new IllegalStateException("Unhandled case!");
			
//			if(primerType.equals(PrimerTypes.forwardPrimer)) scanRegion = backgroundIndex.getSequence().substring(start, start + length);
//			else if(primerType.equals(PrimerTypes.reversePrimer)) scanRegion = backgroundIndex.getSequence().substring(start - length, start);
//			else throw new IllegalArgumentException("Unsupported primer type!");
			Sequence dnaSeq = null;
			Sequence newSeq;
			try{
				dnaSeq = DNATools.createDNASequence(scanRegion, "");
			}catch(IllegalSymbolException e){
				e.printStackTrace();
			}
			synchronized (mapper) {
				newSeq = mapper.annotate(dnaSeq);
			}
			
			Iterator iter = newSeq.features();
			// skip MAX_PRIMER_MISPRIMING_CUTOFF 'acceptable' misprimings
			for(int j=0; j < searchParams.getMAX_PRIMER_MISPRIMING_CUTOFF(); j++)	iter.next();
			// check remaining misprimings
//			if(iter.hasNext() &! backgroundIndex.includesScanRegion()){
//				// if mispriming scan region does not include primer scan region, NO mispriming match is allowed but at least one was found
//				result = false;
//				break;
//			}else if(iter.hasNext() && backgroundIndex.includesScanRegion()){
//				iter.next();
//				if(iter.hasNext()){
//					// if mispriming scan region includes primer scan region, exactly one pseudo-mispriming match is allowed
//					result =  false;
//					break;
//				}else{
//					result = true;
//				}
//			}
			//if(searchParams.isPRINT_DEBUG_LOG()) 
				//System.err.print("Mispriming - forward?: " + mispriming.isForwardHit() + " RSSs found in vicinity: " + newSeq.countFeatures());
			
			if(iter.hasNext() && !truePositivePrimingFound){
				truePositivePrimingFound = true;
				iter.next();
				if(iter.hasNext()){
					//System.err.println(" vote: UNSAFE");
					if(searchParams.isPRINT_DEBUG_LOG()){
						System.err.println("\tUNSAFE Mispriming(s) at " + mispriming.getContigName() + " at position: " + position);
						for(int j=0; j<hits.size(); j++){
							System.err.println("\t" + ((IndexHit)hits.get(j)).getContigName() + " " + ((IndexHit)hits.get(j)).getPosition());
						}
					}
					return false;
				}
				else{
					//System.err.println();
				}
			}
			else if(iter.hasNext()){
				//System.err.println(" vote: UNSAFE");
				if(searchParams.isPRINT_DEBUG_LOG()){
					System.err.println("\tUNSAFE Mispriming(s) at " + mispriming.getContigName() + " at position: " + position);
					for(int j=0; j<hits.size(); j++){
						System.err.println("\t" + ((IndexHit)hits.get(j)).getContigName() + " " + ((IndexHit)hits.get(j)).getPosition());
					}
				}
				return false;
			}
			else{
				//System.err.println();
				result = true;
				if(searchParams.isPRINT_DEBUG_LOG()){
					System.err.println("SAFE mispriming at: " + mispriming.getContigName() + " - " + mispriming.getPosition());
				}
			}
		}
		//System.err.println("Mispriming vote: SAFE");
		if(searchParams.isPRINT_DEBUG_LOG() && positions.size() == 0) System.err.println("\tSAFE Mispriming");
		if(searchParams.isPRINT_DEBUG_LOG()) System.err.println("Total result (checking " + hits.size() + " candidate misprimings) - is safe?: " + result);
		return result;
	}
}
