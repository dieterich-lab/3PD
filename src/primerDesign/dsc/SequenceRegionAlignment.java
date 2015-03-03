package primerDesign.dsc;


/**
 * Encapsulates a sequence region alignment and the associated methods to extract subalignments.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SequenceRegionAlignment {
	private int[][] globalAlignment;
	private static final int SEA_LENGTH = 5;
	
	/**
	 * Initializes a sequence region alignment.
	 * 
	 * @param globalAlignment the global alignment matrix
	 */
	public SequenceRegionAlignment(int[][] globalAlignment){
		if(globalAlignment.length < 1 || globalAlignment[0].length < 1){
			throw new IllegalArgumentException("Global alignment matrix must have size >= 1x1!");
		}
		this.globalAlignment = globalAlignment;
	}
	
	/**
	 * Returns the maximum global alignment value of the endings of an alignment.
	 * 
	 * An alignment ending is defined to be the last column and the last row of an alignment matrix.
	 * 
	 * @param positionFirst the start position of the first sequence in the alignment matrix
	 * @param lengthFirst the length of the first sequence
	 * @param positionSecond the start position o the second sequence in the alignment matrix
	 * @param lengthSecond the length of the second sequence
	 * 
	 * @return the max alignment value of the two sequences extracted from the alignment of a larger region
	 */
	public PrimerAlignmentScores getGlobalAlignmentValues(int positionFirst, int lengthFirst, int positionSecond, int lengthSecond){
		if(positionFirst + lengthFirst >= this.globalAlignment.length || positionSecond + lengthSecond >= this.globalAlignment[0].length) throw new ArrayIndexOutOfBoundsException();
		return getAlignmentValue(this.globalAlignment, positionFirst, lengthFirst, positionSecond, lengthSecond);
	}
	
	/**
	 * Returns the maximum alignment value of the endings of an alignment.
	 * 
	 * An alignment ending is defined to be the last column and the last row of an alignment matrix.
	 * 
	 * @param matrix the alignment matrix to extract the alignment value from
	 * @param positionFirst the start position of the first sequence in the alignment matrix
	 * @param lengthFirst the length of the first sequence
	 * @param positionSecond the start position o the second sequence in the alignment matrix
	 * @param lengthSecond the length of the second sequence
	 *
	 * @return the max alignment value of the two sequences extracted from the alignment of a larger region
	 */
	private PrimerAlignmentScores getAlignmentValue(int[][] matrix, int positionFirst, int lengthFirst, int positionSecond, int lengthSecond){
		assert(positionFirst >= 0 && positionFirst < matrix.length - 1 && positionSecond >= 0 && positionSecond < matrix[0].length - 1);
		// retrieve alignment value from local int[][]:
		// in case this is an end alignment, only scan for the maximum of the lowermost row in the alignment
		// else also scan for the maximum of the rightmost column in the alignment			
//		int sa_max = Integer.MIN_VALUE;
//		int sea_max = Integer.MIN_VALUE;
		
//		final int LAST_LINE = 1 + positionFirst + lengthFirst - 1; // the last line relevant to this alignment
//		final int RIGHTMOST_COLUMN = this.columnsLength - 1 - positionSecond; // the rightmost column relevant to this alignment
//		int subtractor = 0; // the value to subtract from the current alignment value (the first value left or above the current alignment which is NOT relevant to this alignment)
//		int value1 = Integer.MIN_VALUE;
//		int value2 = Integer.MIN_VALUE;
//		int shift = 0;
//		int corrector = 0; // a term to correct for non-quadratic matrices
//		final int COLUMN_BEFORE_FIRST_COLUMN = RIGHTMOST_COLUMN - lengthSecond;
//		// scan lower triangular matrix (pa and pea)
//		for(int i=COLUMN_BEFORE_FIRST_COLUMN + 1 ; i<=RIGHTMOST_COLUMN; i++){
//			if(LAST_LINE - 1 - shift + corrector < 0 || COLUMN_BEFORE_FIRST_COLUMN + corrector < 0){
//				subtractor = 0;
//			}
//			else{
//				subtractor = matrix[LAST_LINE - 1 - shift + corrector][COLUMN_BEFORE_FIRST_COLUMN + corrector];
//			}
//			// calculate pa
//			value1 = matrix[LAST_LINE][i] - subtractor;
//			if(value1 > sa_max) sa_max = value1;
//			// as long as both self-end-alignments completely overlap...
//			if(false && LAST_LINE - (LAST_LINE - 1 - shift - 1) + i - COLUMN_BEFORE_FIRST_COLUMN <= SEA_LENGTH + 1){
//				sea_max = sa_max;
//			}
//			else{
//				subtractor = (LAST_LINE - SEA_LENGTH >= 0 && i-SEA_LENGTH >= 0) ? matrix[LAST_LINE - SEA_LENGTH][i-SEA_LENGTH] : 0;
//				value1 = matrix[LAST_LINE][i] - subtractor;
//				subtractor = matrix[Math.max(LAST_LINE - 1 - shift + corrector, 0)][Math.max(COLUMN_BEFORE_FIRST_COLUMN + corrector, 0)];
//				assert(subtractor >= 0);
//
//				value2 = matrix[Math.min(LAST_LINE - 1 - shift + corrector + SEA_LENGTH, this.rowsLength -1)][Math.min(COLUMN_BEFORE_FIRST_COLUMN + corrector + SEA_LENGTH, this.columnsLength - 1)] - subtractor;
//
//				// self end alignments are only meaningful iff both ends align at least to a certain extent
//				if(value1 > sea_max || value2 > sea_max){
//					sea_max = (value1 > 0 && value2 > 0) ? Math.max(value1, value2) : 0;
//				}
//			}
//			if(shift >= lengthFirst) corrector++;
//			shift++;
//		}
//		shift = 0;
//		corrector = 0;
//		// scan upper triangular matrix (pa)
//		for(int i=LAST_LINE; i>=positionFirst+1 && shift <= lengthSecond; i--){
//			if(positionFirst - 1  < 0 || COLUMN_BEFORE_FIRST_COLUMN + shift < 0){
//				subtractor = 0;
//			}
//			else{
//				subtractor = matrix[positionFirst - 1][COLUMN_BEFORE_FIRST_COLUMN + shift++];
//			}
//			value1 = matrix[i][RIGHTMOST_COLUMN] - subtractor;
//			if(sa_max < value1) sa_max = value1;
//		}
//		return new PrimerAlignmentScores(sa_max, sea_max);
		
		final int LMC = this.globalAlignment[0].length - positionSecond - lengthSecond; // the leftmost column of the alignment sub-matrix
		final int RMC = this.globalAlignment[0].length - 1 - positionSecond; // the rightmost column of the alignment sub-matrix
		final int TMR = 1 + positionFirst; // the topmost row of the alignment sub-matrix
		final int LWMR = 1 + positionFirst + lengthFirst - 1; // the lowermost row of the alignment sub-matrix
		
		int subPosV = 0; // the vertical shift of the subtractor
		int subPosH = 0; // the horizontal shift of the subtractor
		int valPosV = 0; // the vertical shift of the alignment value
		int valPosH = 0; // the horizontal shift of the alignment value
		
		final int END_ALIGNMENT_LENGTH = Math.min(lengthFirst, lengthSecond);
		final int breakH = LMC + END_ALIGNMENT_LENGTH;
		final int breakV = LWMR - END_ALIGNMENT_LENGTH;
		
		boolean isQuadraticAlignment = RMC - LMC == LWMR - TMR;
		
		int sa_final = -1;
		int sea_final = -1;
		
		int candidateSA = -1;
		int candidateSEA = -1;
		
		// traverse lower triangular matrix, compute SEA and SA
		for(int i=0; i<END_ALIGNMENT_LENGTH; i++){
			// check global pair alignment value
			candidateSA = this.globalAlignment[LWMR][LMC + i] - this.globalAlignment[LWMR - 1 - i][LMC - 1];
			if(candidateSA > sa_final) sa_final = candidateSA;
			
			// if there exist two non-completely-overlapping pair end alignments, scan each of them
			if(i >= SEA_LENGTH){
				// check first global pair end alignment value
				//candidateSEA = this.globalAlignment[LWMR][LMC + i] - ((LWMR - SEA_LENGTH >= 0 && LMC + i - SEA_LENGTH >= 0) ? this.globalAlignment[LWMR - SEA_LENGTH][LMC + i - SEA_LENGTH] : 0);
				candidateSEA = this.globalAlignment[LWMR][LMC + i] - ((LWMR - SEA_LENGTH + 1 >= 0 && LMC + i - SEA_LENGTH + 1 >= 0) ? this.globalAlignment[LWMR - SEA_LENGTH + 1][LMC + i - SEA_LENGTH + 1] : 0);
				if(candidateSEA > sea_final) sea_final = candidateSEA;
				
				// check second global alignment value
				//candidateSEA = this.globalAlignment[LWMR - 1 - i + SEA_LENGTH][LMC - 1 + SEA_LENGTH] - this.globalAlignment[LWMR - 1 - i][LMC - 1];
				candidateSEA = this.globalAlignment[LWMR - 1 - i + SEA_LENGTH][LMC - 1 + SEA_LENGTH] - this.globalAlignment[LWMR - i][LMC];
				if(candidateSEA > sea_final) sea_final = candidateSEA;
			}
			// if both end alignments completely overlap
			else if(candidateSA > sea_final){
				sea_final = candidateSA;
			}
		}
		
		// traverse upper triangular matrix, compute SA
		if(!isQuadraticAlignment){
			if(lengthFirst > lengthSecond){
				subPosH = LMC - 1;
				subPosV = breakV;
				valPosH = RMC;
				valPosV = LWMR - 1;
			}else if(lengthSecond > lengthFirst){
				subPosH = LMC;
				subPosV = TMR - 1;
				valPosH = breakH;
				valPosV = LWMR;
			}
		}else{
			subPosH = LMC - 1;
			subPosV = TMR - 1;
			valPosH = RMC;
			valPosV = LWMR;
		}
			
		final int REST_LENGTH = Math.max(lengthFirst, lengthSecond);
		for(int i=0; i<REST_LENGTH; i++){
			candidateSA = this.globalAlignment[valPosV][valPosH] - this.globalAlignment[subPosV][subPosH];
			if(candidateSA > sa_final) sa_final = candidateSA;
			
			if(subPosV > TMR - 1) subPosV--;
			else subPosH++;
			
			if(valPosH < RMC) valPosH++;
			else valPosV--;
		}
		
		return new PrimerAlignmentScores(sa_final, sea_final);
	}
	
	public String toString(){
		StringBuffer buffy = new StringBuffer();
		for(int i=0; i<this.globalAlignment.length; i++){
			for(int j=0; j<this.globalAlignment[i].length; j++){
				buffy.append(this.globalAlignment[i][j]);
				if(j<this.globalAlignment[i].length - 1) buffy.append("\t");
			}
			buffy.append("\n");
		}
		return buffy.toString();
	}
}
