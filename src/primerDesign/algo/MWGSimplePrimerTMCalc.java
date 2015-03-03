/**
 * 
 */
package primerDesign.algo;

/**
 * This class implements a simple primer melting temp. calculation as employed by MWG Biotech.
 * 
 * @author Sebastian Fršhler
 *
 */
public class MWGSimplePrimerTMCalc implements PrimerMeltingTempCalculation {

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerMeltingTempCalculation#computeTM(java.lang.String)
	 */
	public double computeTM(String primer) {
		
		primer = primer.toUpperCase();
		
		double gcCount = 0;
		for(int i=0; i<primer.length(); i++){
			if(primer.charAt(i) == 'G' || primer.charAt(i) == 'C') gcCount++;
		}
		
		return 69.3 + 41 * (gcCount/primer.length()) - 650/primer.length();
	}
	

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerMeltingTempCalculation#computeMonoAndDivalentCationCorrectedTM(java.lang.String, double, double, double, double)
	 */
	public double computeMonoAndDivalentCationCorrectedTM(String primr, double primerConcentration, double cationConcentration, double divalentCationConcentration, double dNTPConcentration) {
		throw new IllegalStateException("This method does not support concentration-dependent primer TM calculations!");
	}

	public static void main(String[] args){
		String primer = args[0];
		MWGSimplePrimerTMCalc calc = new MWGSimplePrimerTMCalc();
		System.out.println(primer + "\t" + calc.computeTM(primer));
	}
}
