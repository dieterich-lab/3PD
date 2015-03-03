/**
 * 
 */
package primerDesign.Test;

import primerDesign.algo.MWGSimplePrimerTMCalc;
import primerDesign.algo.PrimerMeltingTempCalculation;
import primerDesign.algo.SantaLuciaTM;
import primerDesign.util.Constants;
import primerDesign.util.SeqTools;

/**
 * @author froehler
 *
 */
public class PrimerTMCalcComparison {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length == 0) System.out.println("Usage: PrimerTMCalcComparison <#random primers for TM comparison");
		int number = Integer.parseInt(args[0]);
		PrimerMeltingTempCalculation mwg = new MWGSimplePrimerTMCalc();
		PrimerMeltingTempCalculation lucia = new SantaLuciaTM();
		String primer;
		String separator = "\t";
		
		for(int i=0; i<number; i++){
			primer = SeqTools.getRandomPrimerSequence(Constants.MIN_PRIMER_LENGTH, Constants.MAX_PRIMER_LENGTH);
			System.out.println(primer + separator + mwg.computeTM(primer) + separator + lucia.computeTM(primer));
		}
	}

}
