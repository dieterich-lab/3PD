/**
 * 
 */
package primerDesign.Test;

import java.util.HashMap;


/**
 * Implements some primer stability calculations as in SantaLucia: Proc.Natl.Acad.Sci. 1998 p.1460ff.
 * 
 * The implementation uses the unified nearest neighbor parameters!
 * 
 * @author froehler
 *
 */
public class PrimerStabilityCalc {
	private static HashMap<String, Double> deltaG = new HashMap<String, Double>();
	private static boolean isInitialized = false;

	
//	private static double threeBp_hairpin_inc = 5.2;
//	private static double fourBp_hairpin_inc = 4.5;
//	private static double fiveBp_hairpin_inc = 4.4;
//	private static double sixBp_hairpin_inc = 4.3;
//	private static double sevenBp_hairping_inc = 4.1;
//	private static double eightBp_hairpin_inc = 4.1;
	
	public PrimerStabilityCalc(){
		if(!isInitialized) init();
	}
	
	private void init(){
		deltaG.put("AA", -1.0);
		deltaG.put("AT", -0.88);
		deltaG.put("AG", -1.28);
		deltaG.put("AC", -1.44);
		
		deltaG.put("CA", -1.45);
		deltaG.put("CT", -1.28);
		deltaG.put("CG", -2.17);
		deltaG.put("CC", -1.84);
		
		deltaG.put("GA", -1.3);
		deltaG.put("GT", -1.44);
		deltaG.put("GG", -1.84);
		deltaG.put("GC", -2.24);
		
		deltaG.put("TA", -0.58);
		deltaG.put("TT", -1.0);
		deltaG.put("TG", -1.45);
		deltaG.put("TC", -1.3);
		
		deltaG.put("A", 1.03);
		deltaG.put("T", 1.03);
		deltaG.put("G", 0.98);
		deltaG.put("C", 0.98);
	}
	
	/**
	 * Calculates the stability of a oligonucleotide sequence ('delta G').
	 * 
	 * For simplification, the stability is calculated at 37¡C as in SantaLucia (see above)
	 * 
	 * @param sequence
	 * @return
	 */
	public double calcSelfStability(String primer){
		double result = 0;

		// compute delta G according to Table 1 in SantaLucia et.al. (unified parameters)
		for(int i=0; i < primer.length()-1; i++){
			if(i == 0){
				result += deltaG.get(primer.substring(0, 1));
			}
			else if(i == primer.length()-2){
				// add terminal bp penalty
				result += deltaG.get(primer.substring(i+1, i+2));
			}
			result += deltaG.get(primer.substring(i, i+2));
			
		}
		return result;
	}
	
	/**
	 * Calculates the stability of a candidate hairpin in a sequence.
	 * 
	 * Hairpins are supposed to have length three to eight!
	 * 
	 * @param sequence the sequence to calculate the hairpin stability for
	 * @return the max stability of a candidate hairpin in sequence 'sequence'
	 */
//	public static double calcHairpinStability(String sequence){
//		double result = 0;
//		double temp;
//		
//		for(int i=3; i<=8; i++){
//			for(int j=0; j<=sequence.length()-i; j++){
//				temp += calcStability(sequence.substring(j, j+i));
//				if(i == 3) temp += threeBp_hairpin_inc;
//				else if(i == 4) temp += fourBp_hairpin_inc;
//				else if(i == 5) temp += fiveBp_hairpin_inc;
//				else if(i == 6) temp += sixBp_hairpin_inc;
//				else if(i == 7) temp += sevenBp_hairping_inc;
//				else if(i == 8) temp += eightBp_hairpin_inc;
//				else throw new IllegalStateException("Unhandled case!");
//				
//				
//			}
//			result = Math.max(result, temp);
//		}
//		return result;
//	}
	
	public boolean hasMoreStableFivePrimeEnd(String sequence){
		int mean = sequence.length()/2;
		if(calcSelfStability(sequence.substring(0, mean)) < calcSelfStability(sequence.substring(mean + 1))) return true;
		else return false;
	}
	
//	private static String comparePairs(String seq1, String seq2){
//		StringBuffer buffer = new StringBuffer();
//		
//		buffer.append("First:\n");
//		buffer.append("Stability: " + calcStability(seq1));
//		buffer.append(" Hairpin: " + calcHairpinStability(seq1));
//		buffer.append(" End: " + calcStability(seq1.substring(seq1.length()-5, seq1.length())));
//		buffer.append(" More stable 5' end: " + hasMoreStableFivePrimeEnd(seq1) + " (" + calcStability(seq1.substring(0, seq1.length()/2)) + "/" + calcStability(seq1.substring((seq1.length()/2)+1)) + ")\n");
//		buffer.append("Second\n");
//		buffer.append("Stability: " + calcStability(seq2));
//		buffer.append(" Hairpin: " + calcHairpinStability(seq2));
//		buffer.append(" End: " + calcStability(seq2.substring(seq2.length()-5, seq2.length())));
//		buffer.append(" More stable 5' end: " + hasMoreStableFivePrimeEnd(seq2) + " (" + calcStability(seq2.substring(0, seq2.length()/2)) + "/" + calcStability(seq2.substring((seq2.length()/2)+1)) + ")\n");
//		
//		return buffer.toString();
//	}
	
	public static void main(String[] args){
		PrimerStabilityCalc calc = new PrimerStabilityCalc();
//		System.out.println("TATAT: " + calc.calcSelfStability("TATAT"));
//		System.out.println("GCGCG: " + calc.calcSelfStability("GCGCG"));
		
		System.out.println("CATCG: " + calc.calcSelfStability("CATCG"));
		System.out.println("AATAC: " + calc.calcSelfStability("AATAC"));
		System.out.println();
		System.out.println("CGAAT: " + calc.calcSelfStability("CGAAT"));
		System.out.println("AAATT: " + calc.calcSelfStability("AAATT"));
		System.out.println();
		System.out.println("GAATAG: " + calc.calcSelfStability("GAATAG"));
		System.out.println("AATAC: " + calc.calcSelfStability("AATAC"));
	}
}
