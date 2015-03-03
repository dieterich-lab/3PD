/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.PrimerTypes;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;

/**
 * Implements a dynamic programming based restriction site fragment alignment.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RSSDPAligner {
	
	private int initGapPenalty = 0; // the penalty for an initial gap
	private String forward;
	private int startPos;
//	private String reverse;
//	private String probe;
	private int[][] forwardForward;  // alignment of the whole forward primer scan region
//	private int[][] reverseReverse;  // alignment of the whole reverse primer scan region
//	private int[][] probeProbe;  // alignment of the whole probe primer scan region
//	private int[][] forwardReverse;  // alignment of the whole forward primer scan region with the whole reverse primer scan region
//	private int[][] forwardProbe;  // alignment of the whole forward primer scan region with the whole probe primer scan region
//	private int[][] reverseProbe;  // alignment of the whole reverse primer scan region with the whole probe primer scan region
	private PrimerSearchParameters searchParams;

	/**
	 * Initializes a new restriction site scan regions alignment.
	 * 
	 * @param forward the forward primer scan region in 5'->3' direction
	 * @param startPos the start position of this forward region
	 * @param searchParams the 3PD search parameters
	 */
	public RSSDPAligner(String forward, int startPos, PrimerSearchParameters searchParams){ // , String reverse, String probe){
		this.forward = forward.toUpperCase();
		this.startPos = startPos;
		this.searchParams = searchParams;
//		this.reverse = reverse;
//		this.probe = probe;
		// init dp matrices
		this.forwardForward = computeDPAlignmentMatrix(this.forward, SeqTools.reverseDNA(this.forward), this.forwardForward);
//		this.reverseReverse = computeDPAlignmentMatrix(reverse, SeqTools.reverseDNA(reverse), this.reverseReverse);
//		this.probeProbe = computeDPAlignmentMatrix(probe, SeqTools.reverseDNA(probe), this.probeProbe);
//		this.forwardReverse = computeDPAlignmentMatrix(forward, reverse, this.forwardReverse);
//		this.forwardProbe = computeDPAlignmentMatrix(forward, probe, this.forwardProbe);
//		this.reverseProbe = computeDPAlignmentMatrix(reverse, SeqTools.reverseDNA(probe), this.reverseProbe);
	}
	
	/**
	 * Computes the alignment matrix of region1 and region2 using dynamic programming.
	 * 
	 * @param region1 the first region of the alignment
	 * @param region2 the second region of the alignment
	 * @param matrix the alignment matrix to store computed values in
	 */
	private int[][] computeDPAlignmentMatrix(String region1, String region2, int[][] matrix){
		int length1 = region1.length();
		int length2 = region2.length();
		matrix = new int[length1+1][length2+1];
		for(int i=0; i<=length1; i++) matrix[i][0] = initGapPenalty;
		for(int i=0; i<=length2; i++) matrix[0][i] = initGapPenalty;
		for(int i=1; i<=length1; i++){
			for(int j=1; j<=length2; j++){
				matrix[i][j] = matrix[i-1][j-1] + characterScoring(region1.charAt(i-1), region2.charAt(j-1));
			}
		}
		return matrix;
	}
	
	/**
	 * This implements the character pair scoring function 's(x_i,y_i)'.
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 * 
	 * @param xi the character in sequence 1
	 * @param yi the character in sequence 2
	 * 
	 * @return the score of this character pair
	 */
	private int characterScoring(char xi, char yi){
		if((xi == 'A' && yi == 'T') || (xi == 'T' && yi == 'A')){
			return this.searchParams.getA_t_basepair_score();
		}
		else if((xi == 'G' && yi == 'C') || (xi == 'C' && yi == 'G')){
			return this.searchParams.getG_c_basepair_score();
		}
		else{
			return 0;
		}
	}
	
	public PrimerAlignmentScores computeSelfAlignment(Primer primer){
		
		final int start = Math.abs(primer.getRelativePosition() - this.startPos) + 1;
		final int end = Math.abs(primer.getRelativePosition() - this.startPos) + primer.getLength();
		
		int sa_max = 0;
		int sea_max = 0;
		int pa_current;
		int reverseStart;
		int reverseEnd;
		for(int i=0; i<primer.getLength(); i++){
			// update pa_max
			reverseStart = this.forward.length() - start + 1;
			reverseEnd = reverseStart - primer.getLength() + 1;
			pa_current = this.forwardForward[end-i][reverseStart] - this.forwardForward[start-1][reverseEnd+i-1];
			if(pa_current > sa_max) sa_max = pa_current;
			pa_current = this.forwardForward[end][reverseStart-i] - this.forwardForward[start+i-1][reverseEnd-1];
			if(pa_current > sa_max) sa_max = pa_current;
		}
		
		// update pea_max
		int current_value;
		int prev_value;
		int candidate;
		boolean cont = true;
		
		reverseEnd = this.forward.length() - end + 1;
			
		if(start == 0) sea_max = this.forwardForward[end][reverseEnd];
		else sea_max = this.forwardForward[end][reverseEnd] - this.forwardForward[end-1][reverseEnd-1];
		for(int i=1; i<primer.getLength() && cont; i++){
			prev_value = this.forwardForward[end][reverseEnd+i];
			for(int j=1; j<=i && cont; j++){
				current_value = this.forwardForward[end-j][reverseEnd+i-j];
				if(current_value < prev_value){
					candidate = this.forwardForward[end][reverseEnd+i] - this.forwardForward[end-j-1][reverseEnd+i-j-1];
					if(j==i && candidate > sea_max) sea_max = candidate;
					prev_value = current_value;
				}else break; //cont = false;
			}
		}
		
		return new PrimerAlignmentScores(sa_max, sea_max);
	}
	
	public String toString(int[][] matrix){
		StringBuffer buffer = new StringBuffer();
		
		buffer.append("\t\t");
		for(int i=this.forward.length()-1; i>=0; i--){
			buffer.append(this.forward.charAt(i));
			if(i>0) buffer.append("\t");
		}
		buffer.append("\n");
		
		for(int i=0; i<matrix.length; i++){
			if(i>0) buffer.append(this.forward.charAt(i-1) + "\t");
			else buffer.append("\t");
			for(int j=0; j<matrix[i].length; j++){
				buffer.append(matrix[i][j]);
				if(j<matrix[i].length-1) buffer.append("\t");
			}
			buffer.append("\n");
		}
		return buffer.toString();
	}
	
	public String toString(){
		return this.toString(this.forwardForward);
	}
	
	public static void main(String[] args){
//		RSSDPAligner aligner = new RSSDPAligner("ctccaccacacgccaaacaacatgctcaatatgctcttcaagtagggtatcacggaattgg", 57582); //, "ATAT", "ATAT");
//		System.out.println(aligner.toString());
//		Primer primer = new Primer("TCAATATGCTCTTCAAGTAGGG", PrimerTypes.reversePrimer, aligner, 57557);
//		primer.setRelativePosition(57557);
//		System.out.println(aligner.computeSelfAlignment(primer).getPairScore());
//		System.out.println(aligner.computeSelfAlignment(primer).getPairEndScore());
		PrimerSearchParameters searchParams = new PrimerSearchParameters();
		RSSDPAligner aligner = new RSSDPAligner("ggatatg", 0, searchParams); //, "ATAT", "ATAT");
		System.out.println(aligner.toString());
		Primer primer = new Primer("atat", PrimerTypes.forwardPrimer, aligner, 2, searchParams);
		primer.setRelativePosition(2);
		System.out.println(aligner.computeSelfAlignment(primer).getPairScore());
		System.out.println(aligner.computeSelfAlignment(primer).getPairEndScore());
	}
}
