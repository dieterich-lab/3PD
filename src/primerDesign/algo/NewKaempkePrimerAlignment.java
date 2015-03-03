/**
 * 
 */
package primerDesign.algo;

import java.util.HashMap;

import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.util.Constants;

/**
 * THIS IS A DEADEND CLASS - NOT TO BE FURTHER IMPLEMENTED!!!
 * 
 * @author Sebastian Fršhler
 *
 */
public class NewKaempkePrimerAlignment {
	//private OpenIntIntHashMap skValues = new OpenIntIntHashMap();
	private HashMap<String, Integer> skValues = new HashMap<String, Integer>();


	/**
	 * This implements the DYNSA-APP algorithm ('primer self alignments').
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 *
	 *@param primr the primer for which the self alignment is to be computed in 5'->3' direction
	 */
	public PrimerAlignmentScores computeSelfAlignment(String primr){
		String primer = primr.toUpperCase();
		
			int maxScore_sa = 0;
			
			// initialization
			String w2 = primer.substring(0,2);
			//String w1 = primer.substring(0,2); //(new StringBuilder(w2)).reverse().toString();
			
			int sa_current = alignmentScoring(w2, w2).getPairScore(); //w1, w2).getPairScore();
			if(sa_current > maxScore_sa) maxScore_sa = sa_current;
			
			// iteration
			for(int n=1; n<primer.length()-1;n++){
				String wn = primer.substring(n, n+2);
				//String wn_rev = primer.substring(n, n+2); //(new StringBuilder(wn)).reverse().toString();
	
				sa_current = alignmentScoring(wn, wn).getPairScore(); //wn_rev, wn).getPairScore();
				if(sa_current > maxScore_sa) maxScore_sa = sa_current;
				
				for(int i=1; i<= Math.min(n, Constants.MAX_PRIMER_LENGTH - 2); i++){
					//StringBuilder sequence = new StringBuilder(primer.substring(n-i, n+2));
					String seq = primer.substring(n-i, n+2); //sequence.toString();
					//String seq_rev = primer.substring(n-i, n+2); //sequence.reverse().toString();
					sa_current = alignmentScoring(seq, seq).getPairScore();  //seq_rev, seq).getPairScore();
					
					if(sa_current > maxScore_sa) maxScore_sa = sa_current;
				}
			}
			
			this.saScore = maxScore_sa;
			this.seaScore = this.alignmentScoring(primer, primer).getPairEndScore(); //(new StringBuilder(primer)).reverse().toString(), primer).getPairEndScore();
			
		return new PrimerAlignmentScores(this.saScore, this.seaScore);
	}
	
	// DYNPA-APP2
	
	/**
	 * This implements the alignment scoring function 'S(x,y)'.
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 * 
	 * @param forward sequence  in 5'->3' direction
	 * @param k the value of k used in dynamic programming
	 * 
	 * @return the maximum alignment value of the two sequences
	 */
	private int alignmentScoring(String forward, int k){
		int result = -1;
		if(k <= -1) result =  skValues.get(forward.substring(1, forward.length()));
		else if(k == 0) result =  skValues.get(forward.substring(1, forward.length()-1)) + 2 * characterScoring(forward.charAt(0), forward.charAt(forward.length()-1));
		else if(k >= 1) result =  skValues.get(forward.substring(0, forward.length()-1));
		else throw new IllegalStateException("Unhandled case!");
		
		this.skValues.put(forward, result);
		
		return result;
	}
	
	/**
	 * This implements the character pair scoring function 's(x_i,y_i)'.
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 * 
	 * @param xi the character in sequence 1
	 * @param yi the character in sequence 2
	 * @return the score of this character pair
	 */
	private int characterScoring(char xi, char yi){
		if((xi == 'A' && yi == 'T') || (xi == 'T' && yi == 'A')){
			return Constants.a_t_basepair_score;
		}
		else if((xi == 'G' && yi == 'C') || (xi == 'C' && yi == 'G')){
			return Constants.g_c_basepair_score;
		}
		else{
			return 0;
		}
	}
}
