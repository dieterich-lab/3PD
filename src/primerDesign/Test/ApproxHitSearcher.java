/**
 * 
 */
package primerDesign.Test;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.IndexHit;
import primerDesign.util.PrimerSearchParameters;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class ApproxHitSearcher {
	/**
	 * Implements an approximate search for hits of query 'sequence' against index 'index'.
	 * 
	 * A query is split up into word of size 'wordsize' and seed hits of those smaller words are queried in the index.
	 * Therafter, foreach hit a DP of its 3' end (of length 'threePrimeDPLength') is computed. 
	 * If this alignment value is <= threshold 'threePrimeDPThreshold', a DP for the whole hit region is computed. 
	 * Iff the DP value for this whole alignment is below threshold 'wholeDPThreshold',
	 * a hit is reported (added to the list of hits which is returned by this method).
	 * 
	 * This search strategy is supposed to be analogous to a BLAST search in a pre-existing index.
	 * 
	 * @param sequence the sequence to search hits for
	 * @param index the index to search hits in
	 * @param wordsize the wordsize to search seed-hits with
	 * @param stepsize the stepsize to shift seed queries of size 'wordsize' over the query sequence
	 * @param threePrimeDPLength the length of the 3' DP for verification
	 * @param threePrimeDPThreshold the threshold for 3' DP verification of candidate hits
	 * @param wholeDPThreshold thw threshold for whole-length DP verification of filtered candidate hits
	 * @return
	 */
	public static ObjectArrayList findHits(String sequence, DNASequenceIndex index, int wordsize, int stepsize, int threePrimeDPLength, int threePrimeDPThreshold, int wholeDPThreshold){
		char[] seq = sequence.toCharArray();
		int seqlength = seq.length;
		if(seqlength < wordsize) throw new IllegalArgumentException("Sequence length must be >= wordsize!");
		else if(seqlength < threePrimeDPLength) throw new IllegalArgumentException("Sequence length must be >= threePrimeDPLength!");
		else if(seqlength < threePrimeDPThreshold || seqlength < wholeDPThreshold) throw new IllegalArgumentException("Sequence length must be >= DP thresholds!");
		ObjectArrayList result = new ObjectArrayList();
		ObjectArrayList hits;
		String currentQuery;
		IndexHit currentHit;
		int k;
		int start;
		int stop;
		
		// determine seed candidate hits
		for(int i=sequence.length() - wordsize; i>=0; i-=stepsize){
			
			currentQuery = sequence.substring(i, i+wordsize);
			
			hits = index.findHitPositions(currentQuery);
			
			for(int j=0; j<hits.size(); j++){
				currentHit = (IndexHit)hits.getQuick(j);
				
				if(currentHit.isForwardHit()){
					start = Math.max(0, currentHit.getPosition() - i + seqlength - 1 - threePrimeDPLength + 1);
					stop = Math.min(currentHit.getPosition() - i + seqlength - 1, currentHit.getContigLength() - 1);
					// if the three prime end of a 'full hit' is not within the sequence limits, this cannot be a real full hit!
					if(start > currentHit.getContigLength() - 1 || stop > currentHit.getContigLength() - 1) k = Integer.MAX_VALUE;
					else{
						k = DP.getDP(seq, seqlength-1-threePrimeDPLength+1, seqlength-1, currentHit.getContigSequence(), start, stop);
					}
				}
				else{
					start = Math.max(0, currentHit.getPosition() - i + seqlength - 1 - threePrimeDPLength + 1);
					stop = Math.min(currentHit.getPosition() - i + seqlength - 1, currentHit.getContigLength() - 1);
					// if the three prime end of a 'full hit' is not within the sequence limits, this cannot be a real full hit!
					if(start < 0 || stop < 0) k = Integer.MAX_VALUE;
					else{
						k = DP.getCompDP(seq, seqlength-1-threePrimeDPLength+1, seqlength-1, currentHit.getContigSequence(), start, stop);
					}
				}
				
				if(k > threePrimeDPThreshold){
					// this subquery cannot lead to a match, discard subquery, do not refine and test!
					hits.remove(j);
					//System.err.println("Removing hit " + currentHit.getPosition() + " " + currentHit.isForwardHit());
					if(j-1 < -1){
						j=-1; // counter has to be >= 0 in next iteration
					}else{
						j--;
					}
				}
			}
			
			// shift hit positions
			int newPos;
			for(int j=0; j<hits.size(); j++){
				currentHit = (IndexHit)hits.getQuick(j);
				newPos = currentHit.getPosition() - i;

				if(newPos >= 0){
					currentHit.setPosition(newPos);
				}
				//else throw new IllegalStateException("Invalid state!");
				//else hits.remove(j); // hits 'extend beyond the 5' (left) end end of the sequence'
			}
			
			// append hits to result
			result = appendUnique(result, hits);
			
			// if there at most one hit (the real TP one), all more specific versions of the query can lead to max the same number of hits! -> prevent needless verification of TP hit
			if(hits.size() <= 1) break;
		}
		
		// do full-length DP verification of candidate hits already having acceptable 3' DPs
		for(int i=0; i<result.size();i++){
			currentHit = (IndexHit)result.getQuick(i);
			
			if(currentHit.isForwardHit()){
				start = Math.max(0, currentHit.getPosition());
				stop = Math.max(currentHit.getPosition() + seqlength -1, currentHit.getContigLength() - 1);
				// if the three prime end of a 'full hit' is not within the sequence limits, this cannot be a real full hit!
				if(start > currentHit.getContigLength() - 1 || stop > currentHit.getContigLength() - 1) throw new IllegalStateException("Illegal state!");  // those cases should be filtered out in the previous loop!!!
				else{
					k = DP.getDP(seq, 0, seqlength-1, currentHit.getContigSequence(), start, stop);
				}
			}
			else{
				start = Math.max(0, currentHit.getPosition());
				stop = Math.min(currentHit.getPosition() + seqlength - 1, currentHit.getContigLength() - 1);
				// if the three prime end of a 'full hit' is not within the sequence limits, this cannot be a real full hit!
				if(start < 0 || stop < 0) throw new IllegalStateException("Illegal state!");
				else{
					k = DP.getCompDP(seq, 0, seqlength-1, currentHit.getContigSequence(), start, stop);
				}
			}
			
			if(k > wholeDPThreshold){
				// this subquery cannot lead to a match, discard subquery, do not refine and test!
				result.remove(i);
				//System.err.println("Removing hit " + currentHit.getPosition() + " " + currentHit.isForwardHit());
				if(i-1 < -1){
					i=-1; // counter has to be >= 0 in next iteration
				}else{
					i--;
				}
			}
		}
		
		return result;
	}
	
	private static ObjectArrayList appendUnique(ObjectArrayList list, ObjectArrayList candidates){
		for(int i=0; i<candidates.size(); i++){
			if(!list.contains((IndexHit)candidates.getQuick(i), true)){
				list.add(candidates.getQuick(i));
			}
		}
		return list;
	}
	
	public static void main(String[] args){
		String query = "TCGATTGTCATT";
		PrimerSearchParameters params = new PrimerSearchParameters();
		DNASequenceIndex index = EnhancedSuffixArrayFatOpt.deserialize("/Users/froehler/SEQUENCE_INDICES/Ppacificus-unmasked.fa.ESAFatOpt.esaidx");
		//MultiSeqESAFatOpt index = MultiSeqESAFatOpt.deserialize(new File("/export/Sebastian/PrimerDesign/Testsequenzen/BAC-Ppa50-C09.fa.MultiSeqESAFatOpt.esaidx"));
		int hits = ApproxHitSearcher.findHits(query, index, params.getSEED_WORD_SIZE(), params.getSEED_STEP_SIZE(), params.getTHREE_PRIME_DP_LENGTH(), params.getTHREE_PRIME_DP_THRESHOLD(), params.getWHOLE_DP_THRESHOLD()).size();
		System.out.println("Found " + hits + " hits for query " + query);
	}
}
