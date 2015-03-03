/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleTimer;

/**
 * @author Sebastian Fršhler
 *
 */
public class SimpleAlignment implements PrimerAlignmentCalculation {
	private PrimerAlignmentScores scores;
	private String p1;
	private int at_score;
	private int gc_score;
	private int mismatch_score = 0;
	
	public SimpleAlignment(int at_score, int gc_score){
		this.at_score = at_score;
		this.gc_score = gc_score;
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#computePairAlignment(java.lang.String, java.lang.String)
	 */
	public PrimerAlignmentScores computePairAlignment(String primer1, String primer2) {
		scores = alignmentScoring(new StringBuilder(primer1).reverse().toString().toUpperCase(), primer2.toUpperCase());
		try{
			throw new Exception("This method is inefficient and tupposed NOT to be used anymore!");
		}
		catch(Exception e ){
			e.printStackTrace();
			System.exit(1);
		}
		return scores;
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#computeSelfAlignment(java.lang.String)
	 */
	public PrimerAlignmentScores computeSelfAlignment(String primer) {
		p1 = primer.toUpperCase();
		scores = computePairAlignment(p1, p1);
		return scores;
	}
	
	/**
	 * This implements the alignment scoring function 'S(x,y)'.
	 * 
	 * Reference: Kaempke et.al.: Efficient primer design algorithms, Bioinformatics 17 (3) 2001, pp. 214-225
	 * 
	 * @param x sequence 1 in 3'->5' direction
	 * @param y sequence 2 in 5'->3' direction
	 * 
	 * @return the maximum alignment value of the two sequences
	 */
	private PrimerAlignmentScores alignmentScoring(String x, String y){
		
		int max_sa = 0;
		int max_sea = 0;
		
		//System.out.println("New alignment: " + x + " " + y);
		for(int k=(y.length()-1); k>-x.length();k--){
			int sum_sa = 0;
			int sum_sea = 0;
			boolean sum_sea_finished = false;
			
			for(int i=x.length()-1; i>= 0; i--){
				// 'virtually enlarge sequence y to fit for alignment'
				if(i+k >= 0 && i+k < y.length()){
					int score = characterScoring(x.charAt(i),y.charAt(i+k));
					//System.out.println("Comparing: " + x.charAt(i) + " " + y.charAt(i+k) + " score: " + score);
					sum_sa += score;
					if(!sum_sea_finished && k >= 0){
						if(score <= 0){
							sum_sea_finished = true;
							//sum_sea += 0;
						}else{
							sum_sea += score;
						}
					}
				}
			}
			if(sum_sa > max_sa){
				max_sa = sum_sa;
			}
			if(sum_sea > max_sea){
				max_sea = sum_sea;
			}
			//System.out.println("next round: " + max_sa + " " + max_sea);
		}
		return new PrimerAlignmentScores(max_sa, max_sea);
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
			return this.at_score;
		}
		else if((xi == 'G' && yi == 'C') || (xi == 'C' && yi == 'G')){
			return this.gc_score;
		}
		else{
			return this.mismatch_score;
		}
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#getPAScore()
	 */
	public int getPAScore() {
		return this.scores.getPairScore();
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#getPEAScore()
	 */
	public int getPEAScore() {
		return this.scores.getPairEndScore();
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#getSAScore()
	 */
	public int getSAScore() {
		return this.scores.getPairScore();
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerAlignmentCalculation#getSEAScore()
	 */
	public int getSEAScore() {
		return this.scores.getPairEndScore();
	}

	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
//		String primer = args[0];
//		int number = Integer.parseInt(args[1]);
		
		PrimerSearchParameters searchParams = new PrimerSearchParameters();
		PrimerAlignmentCalculation alignment = new SimpleAlignment(2, 4);
		
		PrimerAlignmentScores scores = alignment.computePairAlignment("WWWWATWWGCA", "WWWWTGWWCAT");
		System.out.println(scores.getPairScore() + " " + scores.getPairEndScore());
		
		scores = alignment.computePairAlignment("CATTATGGGTGGTATGTTGG", "CATTATGGGTGGTATGTTGG");
		System.out.println(scores.getPairScore() + " " + scores.getPairEndScore());
		
		scores = alignment.computePairAlignment("GGATTGATAATGTAATAGG", "GGATTGATAATGTAATAGG");
		System.out.println(scores.getPairScore() + " " + scores.getPairEndScore());
		
//		String randomPrimer;
//		
//		System.out.print("Performing " + number + " primer alignments of primer " + primer);
//		for(int i=0; i<number; i++){
//			randomPrimer = SeqTools.getRandomPrimerSequence(searchParams.getMIN_PRIMER_LENGTH(), searchParams.getMAX_PRIMER_LENGTH());
//			alignment.computeSelfAlignment(randomPrimer);
//		}
//		System.out.println(" - done in " + timer.getTimeString());
	}
}
