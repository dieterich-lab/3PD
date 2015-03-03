package primerDesign.Test;

import primerDesign.algo.SimpleAlignment;
import primerDesign.dsc.SequenceRegionAlignment;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;

/**
 * This class implements an alignment of a whole sequence region of a dna strand.
 * 
 * It features extracting sub-alignments of subsequences.
 * 
 * @author sebastianfrohler
 *
 */
public class SequenceRegionAligner {

    /**
     * Aligns two sequence regions computing both: global and local alignment.
     * 
     * @param first the first sequence region to align
     * @param second the second sequence region to align
     * @return the local and glocal sequence alignment of the above sequence regions
     */
	public static SequenceRegionAlignment alignSequenceRegions(char[] first, char[] second, int atMatchCost, int gcMatchCost){
		// compute alignments
		int firstLength = first.length;
		int secondLength = second.length;
		int increment;
		int[][] localAlignment = new int[firstLength + 1][secondLength + 1];    // matrix is auto-inited with 0s
		for(int i=1; i<=firstLength; i++){
			for(int j=1; j<=secondLength; j++){
				increment = compare(first[i-1], second[secondLength - 1 - j+1], atMatchCost, gcMatchCost);
				// compute global alignment
				localAlignment[i][j] = localAlignment[i-1][j-1] + increment;
				
			}
		}
		
		
		// return alignments object
		return new SequenceRegionAlignment(localAlignment);
	}
	
	private static int compare(char first, char second, int atMatchCost, int gcMatchCost){
		if(first == SeqTools.revcompChar(second)){
			if(first == 'a' || first == 'A' || first == 't' || first == 'T') return atMatchCost;
			else if(first == 'g' || first == 'G' || first == 'c' || first == 'C') return gcMatchCost;
			else throw new IllegalStateException("Unhandled case!");
		}
		else return 0;
	}
	
	private static int min(int a, int b, int c){
		return Math.min(a, Math.min(b, c));
	}
	
	public static void main(String[] args){		
		char[] first = "CTCGTAAGTAGGATGAAAAC".toCharArray();
		char[] second = "CTCGTAAGTAGGATGAAAAC".toCharArray();
		String f = new String(first);
		String s = new String(second);
		
		SimpleAlignment aligner = new SimpleAlignment(2,4);
		
		int iter = 1000000;
		SimpleTimer timer = new SimpleTimer();
		
		System.out.print("Aligning sequence region of length 20");
		SequenceRegionAlignment alignment = SequenceRegionAligner.alignSequenceRegions(first, second,2, 4);
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (new) alignment values");
		for(int i=0; i<iter; i++){
			alignment.getGlobalAlignmentValues(0, 20, 0, 20);
		}
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (old) alignment values");
		for(int i=0; i<iter; i++){
			aligner.computeSelfAlignment(f);
		}
		System.out.println(" - took " + timer.getTimeString());
		
		first = "CTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAAC".toCharArray();
		second = "CTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAAC".toCharArray();
		f = new String(first).substring(20);
		
		System.out.print("Aligning sequence region of length 40");
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second,2, 4);
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (new) alignment values");
		for(int i=0; i<iter; i++){
			alignment.getGlobalAlignmentValues(0, 20, 0, 20);
		}
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (old) alignment values");
		for(int i=0; i<iter; i++){
			aligner.computeSelfAlignment(f);
		}
		System.out.println(" - took " + timer.getTimeString());
		
		first = "CTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAAC".toCharArray();
		second = "CTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAACCTCGTAAGTAGGATGAAAAC".toCharArray();
		f = new String(first).substring(60);
		
		System.out.print("Aligning sequence region of length 80");
		alignment = SequenceRegionAligner.alignSequenceRegions(first, second,2 ,4);
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (new) alignment values");
		for(int i=0; i<iter; i++){
			alignment.getGlobalAlignmentValues(0, 20, 0, 20);
		}
		System.out.println(" - took " + timer.getTimeString());
		System.out.print("\tGetting " + iter + " (old) alignment values");
		for(int i=0; i<iter; i++){
			aligner.computeSelfAlignment(f);
		}
		System.out.println(" - took " + timer.getTimeString());
	}
}
