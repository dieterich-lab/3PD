/**
 * 
 */
package primerDesign.Test;

import primerDesign.util.Constants;

/**
 * Implements the amplicon TM Calc method as suggested by Burpo et.al.: Biochemistry (2001)
 * @author froehler
 *
 */
public class AmpliconTMCalc {
	
	public static double getAmpliconTM(String sequence){
		double cationConcentration = Constants.MONOVALENT_CATION_CONCENTRATION + 120 * (Math.sqrt(Constants.DIVALENT_CATION_CONCENTRATION-Constants.DNTP_CONCENTRATION));
		return 0.41 * getGCContent(sequence) + 16.6 * Math.log10(cationConcentration) - 675 / sequence.length() + 81.5;
	}

	private static double getGCContent(String sequence){
		double gcContent = 0;
		int maskedCharacters = 0; 
		char maskingChar = Constants.REPETITIVE_ELEMENT_CHARACTER.toCharArray()[0];
		for(int i=0; i < sequence.length(); i++){
			char base = sequence.charAt(i);
			if (base == 'G' || base == 'C'){
				gcContent++;
			}
			else if(base == maskingChar){
				maskedCharacters++;
			}
		}
		return gcContent/ (sequence.length() - maskedCharacters);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println(args[0]);
		System.out.println("%GC: " + getGCContent(args[0]));
		System.out.println("TM: " + getAmpliconTM(args[0]));
	}
}
