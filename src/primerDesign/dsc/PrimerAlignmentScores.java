package primerDesign.dsc;

/**
 * Encapsulates an alignment (the scores of an alignment).
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerAlignmentScores {
	private int pairScore;
	private int pairEndScore;
	
	/**
	 * Initializes a new Primer alignment score.
	 * 
	 * @param pairScore the score for primer pair alignment
	 * @param pairEndScore the score for primer pair-end alignment
	 */
	public PrimerAlignmentScores(int pairScore, int pairEndScore){
		this.pairScore = pairScore;
		this.pairEndScore = pairEndScore;
	}

	/**
	 * Returns the pair-end-score of this alignment.
	 * 
	 * @return the pair-end-score of this alignment
	 */
	public int getPairEndScore() {
		return pairEndScore;
	}

	/**
	 * Returns the pair-score of this alignment.
	 * 
	 * @return the pair-score of this alignment
	 */
	public int getPairScore() {
		return pairScore;
	}
	
	/**
	 * Returns the optimal alignment scores.
	 * 
	 * @return the optimal alignment scores
	 */
	public static PrimerAlignmentScores getOptimalAlignmentScores(){
		return new PrimerAlignmentScores(0,0);
	}
	
	/**
	 * Returns the maximum scores from a list.
	 * 
	 * @param list the list
	 * @return the maximum scores from a list
	 */
	public static PrimerAlignmentScores getMaxScores(PrimerAlignmentScores[] list){
		int sa = -1;
		int sea = -1;
		for(PrimerAlignmentScores scores : list){
			if(sa < scores.getPairScore()) sa = scores.getPairScore();
			if(sea < scores.getPairEndScore()) sea = scores.getPairEndScore();
		}
		return new PrimerAlignmentScores(sa,sea);
	}
}
