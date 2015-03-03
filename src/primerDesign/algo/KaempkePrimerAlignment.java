package primerDesign.algo;

import java.util.HashMap;

import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;

/**
 * Alignment between two PCR-primers.
 * 
 * @author Sebastian Fršhler
 *
 */
public class KaempkePrimerAlignment implements PrimerAlignmentCalculation {
	private int saScore;
	private int seaScore;
	private int paScore;
	private int peaScore;
	private static HashMap<String, Integer> skValues_sa_hash = new HashMap<String, Integer>();
	private static HashMap<String, Integer> skValues_sea_hash = new HashMap<String, Integer>();
	private static HashMap<String, Integer> dynsa_APP_hash_sa = new HashMap<String, Integer>(); 
	private static HashMap<String, Integer> dynsa_APP_hash_sea = new HashMap<String, Integer>();
	private static HashMap<String, Integer> dynpa_APP_hash_pa = new HashMap<String, Integer>(); 
	private static HashMap<String, Integer> dynpa_APP_hash_pea = new HashMap<String, Integer>();
	
	private PrimerSearchParameters params;
	
	public KaempkePrimerAlignment(PrimerSearchParameters params){
		this.params = params;
	}
	
	/**
	 * Returns the sa (self alignment) score of a  primer.
	 * 
	 * @return the sa score of a primer
	 */
	public int getSAScore(){
		return saScore;
	}
	
	/**
	 * Returns the sea (self end alignment) score of a primer.
	 * 
	 * @return the sea score of a primer
	 */
	public int getSEAScore(){
		return seaScore;
	}
	
	/**
	 * Returns the pa (pair alignment) score of two primers.
	 * 
	 * @return the pa score of two primers
	 */
	public int getPAScore(){
		return paScore;
	}
	
	/**
	 * Returns the pea (pair end alignment) score of two primers.
	 * 
	 * @return the pea score of two primers
	 */
	public int getPEAScore(){
		return peaScore;
	}
	
	/**
	 * This implements the DYNSA-APP algorithm ('primer self alignments').
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 *
	 * @param primr the primer for which the self alignment is to be computed in 5'->3' direction
	 */
	public PrimerAlignmentScores computeSelfAlignment(String primr){
		String primer = primr.toUpperCase();
		
		PrimerAlignmentScores alignment = null;
		
		// if alignment values for this primer have already been computed, re-use these
		if(this.params.isUseIndices() && KaempkePrimerAlignment.dynsa_APP_hash_sa.containsKey(primer) && KaempkePrimerAlignment.dynsa_APP_hash_sea.containsKey(primer)){
			this.saScore = KaempkePrimerAlignment.dynsa_APP_hash_sa.get(primer);
			this.seaScore = KaempkePrimerAlignment.dynsa_APP_hash_sea.get(primer);
			alignment = new PrimerAlignmentScores(this.saScore, this.seaScore);
		}else{
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
				
				for(int i=1; i<= Math.min(n, this.params.getMAX_PRIMER_LENGTH() - 2); i++){
					//StringBuilder sequence = new StringBuilder(primer.substring(n-i, n+2));
					String seq = primer.substring(n-i, n+2); //sequence.toString();
					//String seq_rev = primer.substring(n-i, n+2); //sequence.reverse().toString();
					sa_current = alignmentScoring(seq, seq).getPairScore();  //seq_rev, seq).getPairScore();
					
					if(sa_current > maxScore_sa) maxScore_sa = sa_current;
				}
			}
			
			this.saScore = maxScore_sa;
			this.seaScore = this.alignmentScoring(primer, primer).getPairEndScore(); //(new StringBuilder(primer)).reverse().toString(), primer).getPairEndScore();
			
			alignment = new PrimerAlignmentScores(this.saScore, this.seaScore);
			
			if(this.params.isUseIndices()){
				KaempkePrimerAlignment.dynsa_APP_hash_sa.put(primr, this.saScore);
				KaempkePrimerAlignment.dynsa_APP_hash_sea.put(primr, this.seaScore);
			}
		}
		return alignment;
	}
	
	/**
	 * This implements the DYNPA-APP algorithm ('primer1-primer2 alignments').
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 *
	 * @param primr1 the first primer for which the pair alignment is to be computed in 5'->3' direction
	 * @param primr2 the second primer for which the pair alignment is to be computed in 5'->3' direction
	 */
	public PrimerAlignmentScores computePairAlignment(String primr1, String primr2){
		String primer1 = primr1.toUpperCase();
		String primer2 = primr2.toUpperCase();
		
		PrimerAlignmentScores alignment = null;
		
		// if alignment values for this primer pair have already been computed, re-use these
		if(this.params.isUseIndices() && KaempkePrimerAlignment.dynpa_APP_hash_pa.containsKey(primer1 + primer2) && KaempkePrimerAlignment.dynpa_APP_hash_pea.containsKey(primer1 + primer2)){
			this.paScore = KaempkePrimerAlignment.dynpa_APP_hash_pa.get(primer1 + primer2);
			this.peaScore = KaempkePrimerAlignment.dynpa_APP_hash_pea.get(primer1 + primer2);
			alignment = new PrimerAlignmentScores(this.paScore, this.peaScore);
		}else if(this.params.isUseIndices() && KaempkePrimerAlignment.dynpa_APP_hash_pa.containsKey(primer2 + primer1) && KaempkePrimerAlignment.dynpa_APP_hash_pea.containsKey(primer2 + primer1)){
			this.paScore = KaempkePrimerAlignment.dynpa_APP_hash_pa.get(primer2 + primer1);
			this.peaScore = KaempkePrimerAlignment.dynpa_APP_hash_pea.get(primer2 + primer1);
			alignment = new PrimerAlignmentScores(this.paScore, this.peaScore);
		}else{
			int max_Score_pa = 0;
			String w = primer1.substring(0, 2); //(new StringBuilder(primer1.toUpperCase().substring(0, 2))).reverse().toString();
			String v = primer2.toUpperCase().substring(0, 2);
			
			// initialization
			int pa_current = alignmentScoring(w, v).getPairScore();
			if(pa_current > max_Score_pa) max_Score_pa = pa_current;
			
			// iteration
			int N = primer1.length() - 1; // w
			int M = primer2.length() - 1; // v
			for(int k=1; k<=Math.max(N, M)-1; k++){
				if(k <= N-1){
					for(int i=0; i<= Math.min(k, this.params.getMAX_PRIMER_LENGTH() -2); i++){
						String w_current = primer1.substring(k-i, k+2); //(new StringBuilder(primer1.substring(k-i, k+2))).reverse().toString();
						for(int j=0; j<=Math.min(k, M); j++){
							String q_current = primer2.substring(0, j+1);
							pa_current = alignmentScoring(w_current, q_current).getPairScore();
							
							if(pa_current > max_Score_pa) max_Score_pa = pa_current;
						}
					}
				}
				if(k <= M-1){
					for(int i=0; i<=Math.min(k, this.params.getMAX_PRIMER_LENGTH() -2); i++){
						String v_current = primer2.substring(k-i, k+2);
						for(int j=0; j<=Math.min(k+1, N); j++){
							String w_current = primer1.substring(0,j+1); //(new StringBuilder(primer1.substring(0,j+1))).reverse().toString();
							pa_current = alignmentScoring(w_current, v_current).getPairScore();
							
							if(pa_current > max_Score_pa) max_Score_pa = pa_current;
						}
					}
				}
			}
			this.paScore = max_Score_pa;
			this.peaScore = this.alignmentScoring(primer1, primer2).getPairEndScore(); //(new StringBuilder(primer1)).reverse().toString(), primer2).getPairEndScore();
			
			alignment = new PrimerAlignmentScores(this.paScore, this.peaScore);
			
			if(this.params.isUseIndices()){
				KaempkePrimerAlignment.dynpa_APP_hash_pa.put(primer1 + primer2, this.paScore);
				KaempkePrimerAlignment.dynpa_APP_hash_pea.put(primer1 + primer2, this.peaScore);
			}
		}
		return alignment;
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
			return this.params.getA_t_basepair_score();
		}
		else if((xi == 'G' && yi == 'C') || (xi == 'C' && yi == 'G')){
			return this.params.getG_c_basepair_score();
		}
		else{
			return 0;
		}
	}
	
	/**
	 * This implements the alignment scoring function 'S(x,y)'.
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 * 
	 * @param x sequence 1 in 5'->3' direction
	 * @param y sequence 2 in 5'->3' direction
	 * 
	 * @return the maximum alignment value of the two sequences
	 */
	private PrimerAlignmentScores alignmentScoring(String x, String y){
		
		PrimerAlignmentScores scores = null;
		
		int max_sa = 0;
		int max_sea = 0;
		
		// retrieve pre-computed value iff exists
		if(this.params.isUseIndices() && skValues_sa_hash.containsKey(x + y)){
			max_sa = skValues_sa_hash.get(x + y);
			max_sea = skValues_sea_hash.get(x + y);
			scores = new PrimerAlignmentScores(max_sa, max_sea);
		}else if(this.params.isUseIndices() && skValues_sa_hash.containsKey(y + x)){
			max_sa = skValues_sa_hash.get(y + x);
			max_sea = skValues_sea_hash.get(y + x);
			scores = new PrimerAlignmentScores(max_sa, max_sea);
		}else{
			//for(int k=(y.length()-1); k>-x.length();k--){
			int lengthX = x.length();
			int lengthY = y.length();
			int length = Math.min(x.length(), y.length());

			for(int i=0; i<length; i++){
				int sum_sa_up = 0;
				int sum_sea = 0;
				int sum_sa_down = 0;
				boolean sum_sea_finished = false;
				
//				for(int i=x.length()-1; i>= 0; i--){
//					// 'virtually enlarge sequence y to fit for alignment'
//					if(i+k >= 0 && i+k < y.length()){
//						int score = characterScoring(x.charAt(i),y.charAt(i+k));
//						sum_sa += score;
//						if(!sum_sea_finished && k >= 0){
//							if(score == 0){
//								sum_sea_finished = true;
//							}else{
//								sum_sea += score;
//							}
//						}
//					}
//				}
				int score_up = 0;
				int score_down = 0;
				for(int k=0; k<=i; k++){
					score_up = characterScoring(x.charAt(lengthX-1-k), y.charAt(lengthY-1-i+k));
					sum_sa_up += score_up;
					if(!sum_sea_finished){
						if(score_up == 0){
							sum_sea_finished = true;
						}else{
							sum_sea += score_up;
						}
					}
					if(i<length-1){
						score_down = characterScoring(x.charAt(0+k), y.charAt(0+i-k));
						sum_sa_down += score_down;
					}
				}
				if(sum_sa_up > max_sa){
					max_sa = sum_sa_up;
				}
				if(sum_sea > max_sea){
					max_sea = sum_sea;
				}
				if(sum_sa_down > max_sa){
					max_sa = sum_sa_down;
				}
			}
			if(this.params.isUseIndices()){
				skValues_sa_hash.put(x + y, max_sa);
				skValues_sea_hash.put(x + y, max_sea);
			}
			scores = new PrimerAlignmentScores(max_sa, max_sea);
		}		
		return scores;
	}
	
	public PrimerSearchParameters getParams(){
		return this.params;
	}
	
	/**
	 * Returns the maximum value of an integer array.
	 * 
	 * @param array the integer array
	 * 
	 * @return the maximum of the respective integer array
	 */
	private int max(int[] array){
		int max = 0;
		for(int i=0; i<array.length; i++){
			if(array[i] > max){
				max = array[i];
			}
		}
		return max;
	}
	
	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
		String primer = args[0];
		int number = Integer.parseInt(args[1]);
		
		PrimerSearchParameters params = new PrimerSearchParameters();
		PrimerAlignmentCalculation alignment = new KaempkePrimerAlignment(params);
		String randomPrimer;
		
		System.out.print("Performing " + number + " primer alignments of primer " + primer);
		for(int i=0; i<number; i++){
			randomPrimer = SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH());
			alignment.computeSelfAlignment(randomPrimer);
		}
		System.out.println(" - done in " + timer.getTimeString());
	}
}
