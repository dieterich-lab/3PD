/**
 * 
 */
package primerDesign.Test;

import primerDesign.algo.SantaLuciaTM;
import primerDesign.util.PrimerSearchParameters;

/**
 * @author froehler
 *
 */
public class TmMinMaxPrinter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int start = 10;
		int end = 30;
		String lowSeq = "";
		String highSeq = "";
		
		for(int i=0; i<start; i++){
			lowSeq += "A";
			highSeq += "G";
		}
		
		SantaLuciaTM calc = new SantaLuciaTM();
		PrimerSearchParameters params = new PrimerSearchParameters();
		System.out.println("Iter\tlow\thigh");
		for(int i=start; i<=end; i++){
			System.out.print(i + ": " + calc.computeMonoAndDivalentCationCorrectedTM(lowSeq, params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION(), params.getDIVALENT_CATION_CONCENTRATION(), params.getDNTP_CONCENTRATION()));
			System.out.print("\t");
			System.out.println(calc.computeMonoAndDivalentCationCorrectedTM(highSeq, params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION(), params.getDIVALENT_CATION_CONCENTRATION(), params.getDNTP_CONCENTRATION()));
			
			lowSeq += "A";
			highSeq += "G";
		}

	}

}
