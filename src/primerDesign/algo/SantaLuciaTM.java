package primerDesign.algo;

import java.util.HashMap;

import primerDesign.util.Constants;
import primerDesign.util.PrimerSearchParameters;

/**
 * This class implements the melting temperature calculation for a PCR primer according to SantaLucia et.al.
 *
 * Reference: SantaLucia JR. (1998). A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci., 95, 1460-65.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SantaLuciaTM implements PrimerMeltingTempCalculation {

	private static HashMap<String, Double> deltaH = new HashMap<String, Double>();
	private static HashMap<String, Double> deltaS = new HashMap<String, Double>();
	private static boolean isInitialized = false;
	private static final int kilocalories = 1000;
	private static final int temperatureDeltaG = 310;
	private static final PrimerSearchParameters SEARCH_PARAMS = new PrimerSearchParameters();
	
	public SantaLuciaTM(){
		if(!isInitialized){
			this.init();
		}
	}
	
	/**
	 * Initializes deltaH and deltaS parameters.
	 *
	 */
	public void init(){
		// set values for deltaH
		deltaH.put("AA", -7.9);
		deltaH.put("AT", -7.2);
		deltaH.put("AG", -7.8);
		deltaH.put("AC", -8.4);
		
		deltaH.put("TA", -7.2);
		deltaH.put("TT", -7.9);
		deltaH.put("TG", -8.5);
		deltaH.put("TC", -8.2);
		
		deltaH.put("CA", -8.5);
		deltaH.put("CT", -7.8);
		deltaH.put("CG", -10.6);
		deltaH.put("CC", -8.0);
		
		deltaH.put("GA", -8.2);
		deltaH.put("GT", -8.4);
		deltaH.put("GG", -8.0);
		deltaH.put("GC", -9.8);
		

		deltaH.put("A", 2.3);
		deltaH.put("T", 2.3);
		deltaH.put("G", 0.1);
		deltaH.put("C", 0.1);
		
		// set values for deltaS
		deltaS.put("AA", -22.2);
		deltaS.put("AT", -20.4);
		deltaS.put("AG", -21.0);
		deltaS.put("AC", -22.4);
		
		deltaS.put("TA", -21.3);
		deltaS.put("TT", -22.2);
		deltaS.put("TG", -22.7);
		deltaS.put("TC", -22.2);
		
		deltaS.put("CA", -22.7);
		deltaS.put("CT", -21.0);
		deltaS.put("CG", -27.2);
		deltaS.put("CC", -19.9);
		
		deltaS.put("GA", -22.2);
		deltaS.put("GT", -22.4);
		deltaS.put("GG", -19.9);
		deltaS.put("GC", -24.4);
		

		deltaS.put("A", 4.1);
		deltaS.put("T", 4.1);
		deltaS.put("G", -2.8);
		deltaS.put("C", -2.8);
		
		SantaLuciaTM.isInitialized = true;
	}
	
	/**
	 * Computes the melting temperature of a primer.
	 * 
	 * @param primr the primer
	 * @param primerConcentration the primer concentration in mol!
	 * @param cationConcentration the concentration of monovalent cations in mol!
	 * @param divalentCationConcentration the concentration of divalent cations in mol!
	 * @param dNTPConcentration the concentration of dNTPs in mol!
	 * 
	 * @return the melting temperature of primer 'primr'
	 */
	public double computeMonoAndDivalentCationCorrectedTM(String primr, double primerConcentration, double cationConcentration, double divalentCationConcentration, double dNTPConcentration){
		double tm = 0;
		double dH = 0;
		double dS = 0;
		
		if(primr.length() < 2) throw new IllegalArgumentException("Primer length must be >= 2!");
		if(primr.length() > Constants.MAX_PRIMER_LENGTH_FOR_TM_CALC) throw new IllegalArgumentException("Tm calculation for primers is only valid for primer lengths <= " + Constants.MAX_PRIMER_LENGTH_FOR_TM_CALC + "bp");
		if(primerConcentration <= 0 ) throw new IllegalArgumentException("Primer concentration (in mol!) must be > 0!!");
		if(cationConcentration < 0) throw new IllegalArgumentException("The cation concentration must be >= 0!!!");
		if(divalentCationConcentration < 0) throw new IllegalArgumentException("The divalent cation concentration must be >= 0!!!");
		if(dNTPConcentration < 0) throw new IllegalArgumentException("The dNTP concentration must be >= 0!!!");
		
		String primer = primr.toUpperCase();
		if(divalentCationConcentration > 0 && dNTPConcentration > 0){
			cationConcentration += divalentToMonovalent(divalentCationConcentration*1000, dNTPConcentration*1000)/1000;
		}
		
		// symmetry correction as in Primer3 & SantaLucia et.al.
		if(isSymmetric(primr)) dS += -1.4;
		
		// compute delta S and delta H according to Table 2 in SantaLucia et.al.
		for(int i=0; i < primer.length()-1; i++){
			if(i == 0){
				// add initial bp penalty
				dH += deltaH.get(primer.substring(0, 1));
				dS += deltaS.get(primer.substring(0, 1));
			}
			if(i == primer.length()-2){
				// add terminal bp penalty
				dH += deltaH.get(primer.substring(i+1, i+2));
				dS += deltaS.get(primer.substring(i+1, i+2));
			}
			// add dinucleotide penalty
			try{
			dH += deltaH.get(primer.substring(i, i+2));
			dS += deltaS.get(primer.substring(i, i+2));
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		
		if(cationConcentration > 0){	
			dS += 0.368 * (primer.length() - 1) * Math.log(cationConcentration);
		}
			
		boolean symmetry = isSymmetric(primr);
		if(symmetry){
			tm = dH * SantaLuciaTM.kilocalories/(dS + Constants.GAS_CONST_R * Math.log(primerConcentration)) + Constants.ABSOLUTE_ZERO_TEMP;
		}
		else{
			tm = dH * SantaLuciaTM.kilocalories/(dS + Constants.GAS_CONST_R * Math.log(primerConcentration/4)) + Constants.ABSOLUTE_ZERO_TEMP;
		}
		return tm;
	}
	
	/**
	 * Computes the worst-case (optimal, complete binding) binding energy of the oligonucleotide sequence provided (in kcal/mol) at 37¡C.
	 * 
	 * @param primer the oligonucleotide sequence
	 * @param primerConcentration the oligonucleotide concentration
	 * @param cationConcentration the monovalen cation concentration
	 * @param divalentCationConcentration the divalent cation concentration
	 * @param dNTPConcentration the dNTP concentration
	 * 
	 * @return the worst-case (optimal, complete binding) binding energy of the oligonucleotide sequence provided (in kcal/mol) at 37¡C
	 */
	public double computeDeltaG(String primer, double primerConcentration, double cationConcentration, double divalentCationConcentration, double dNTPConcentration){
		
		primer = primer.toUpperCase();
		
		double dH = 0;
		double dS = 0;
		
		if(divalentCationConcentration > 0 && dNTPConcentration > 0){
			cationConcentration += divalentToMonovalent(divalentCationConcentration*1000, dNTPConcentration*1000)/1000;
		}
		
		// symmetry correction as in Primer3 & SantaLucia et.al.
		if(isSymmetric(primer)) dS += -1.4;
		
		// compute delta S and delta H according to Table 2 in SantaLucia et.al.
		for(int i=0; i < primer.length()-1; i++){
			if(i == 0){
				// add initial bp penalty
				dH += deltaH.get(primer.substring(0, 1));
				dS += deltaS.get(primer.substring(0, 1));
			}
			if(i == primer.length()-2){
				// add terminal bp penalty
				dH += deltaH.get(primer.substring(i+1, i+2));
				dS += deltaS.get(primer.substring(i+1, i+2));
			}
			// add dinucleotide penalty
			try{
			dH += deltaH.get(primer.substring(i, i+2));
			dS += deltaS.get(primer.substring(i, i+2));
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		
		if(cationConcentration > 0){	
			dS += 0.368 * (primer.length() - 1) * Math.log(cationConcentration);
		}
		
		return (dH * SantaLuciaTM.kilocalories - SantaLuciaTM.temperatureDeltaG * dS) / SantaLuciaTM.kilocalories;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerMeltingTempCalculation#computeTM(java.lang.String, double)
	 */
	public double computeTM(String primer, double concentration){
		return computeMonoAndDivalentCationCorrectedTM(primer, concentration, SEARCH_PARAMS.getMONOVALENT_CATION_CONCENTRATION(), SEARCH_PARAMS.getDIVALENT_CATION_CONCENTRATION(), SEARCH_PARAMS.getDNTP_CONCENTRATION());
	}
	
	/**
	 * Convenience method to call cation-corrected TM calculation.
	 * 
	 * @param primer the primer
	 * @param concentration the concentration of the primer
	 * @param cationConcentration the monovalent cation concentration
	 * 
	 * @return the melting temperature of primer 'primer'
	 */
	public double computeCationCorrectedTM(String primer, double concentration, double cationConcentration){
		return computeMonoAndDivalentCationCorrectedTM(primer, concentration, cationConcentration, 0, 0);
	}
	
	/**
	 * Computes if a primer sequence is self-symmetric.
	 * 
	 * @param primer the primer for which symmetry should be computed
	 * 
	 * @return true, iff the primer is self-symmetric
	 */
	private boolean isSymmetric(String primer){		
		if(primer.length() % 2 != 0){
			return false;
		}
		else{
			String inversePrimer = (new StringBuilder(primer)).reverse().toString();
			char forward;
			char inverse;
			for(int i=0; i<primer.length(); i++){
				forward = primer.charAt(i);
				inverse = inversePrimer.charAt(i);
				if(inverse != complementBase(forward)){
					return false;
				}
			}
			return true;
		}
	}
	
	/**
	 * Returns the complementary DNA base of a specific base.
	 * @param base the DNA base to return the complementary base for
	 * 
	 * @return the complementary DNA base of the respective base
	 */
	private char complementBase(char base){
		char ret;
		
		switch(base){
			case 'A' : ret = 'T'; break;
			case 'T' : ret = 'A'; break;
			case 'G' : ret = 'C'; break;
			case 'C' : ret = 'G'; break;
			default: throw new IllegalArgumentException("The DNA base " + base + " is no valid base!");
		}
		return ret;
	}
	
	/**
	 * Converts divalent cation concentration to monovalent cation concentration.
	 * 
	 * Implemented as in Primer3
	 * 
	 * @param divalent the divalent cation concentration in mol!
	 * @param dntp the dNTP concentration in mol!
	 * 
	 * @return the monovalent cation concentration in mol!
	 * 
	 * @throws IllegalArgumentException if divalent or dntp is < 0 (Concentration must be positive!)
	 */
	private double divalentToMonovalent(double divalent, double dntp){
		
		if(divalent == 0) dntp = 0;
		if(divalent < 0 || dntp < 0) throw new IllegalArgumentException("Arguments divalent AND dntp must be >= 0!!"); 
		if(divalent < dntp){
			/* According to theory, melting temperature doesn't depend on divalent cations */
			divalent = dntp;
		}
		
		return 120 * (Math.sqrt(divalent-dntp));
	}
		
	/**
	 * Some testing of the class
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		SantaLuciaTM lucia = new SantaLuciaTM();
		String primer = "ATATATATATATATATATATATATATATAT";
		double tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 0, 0, 0);
		System.out.println("T_m (pc,0,0,0) for primer: " + primer + " is: " + tm);
		
		primer = "ATATATATATATATATATATATATATATAT";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 50E-03, 0, 0);
		System.out.println("T_m (pc,mv,0,0) for primer: " + primer + " is: " + tm);
		
		primer = "ATATATATATATATATATATATATATATAT";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 0, 4E-03, 0);
		System.out.println("T_m (pc,0,dv,0) for primer: " + primer + " is: " + tm);
		
		primer = "ATATATATATATATATATATATATATATAT";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 0, 4E-03, 0.2E-03);
		System.out.println("T_m (pc,0,dv,dntp) for primer: " + primer + " is: " + tm);
		
		primer = "ATATATATATATATATATATATATATATAT";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 50E-03, 4E-03, 0.2E-03);
		System.out.println("T_m (pc,mv,dv,dntp) for primer: " + primer + " is: " + tm);
		
		primer = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
		tm = lucia.computeTM(primer, 50E-09);
		System.out.println("T_m (pc,0,0,0) for primer: " + primer + " is: " + tm);
		
		primer = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
		tm = lucia.computeCationCorrectedTM(primer, 50E-09, 50E-03);
		System.out.println("T_m mono(pc,mv,0,0) for primer: " + primer + " is: " + tm);
		
		primer = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 50E-03, 0, 0);
		System.out.println("T_m di(pc,mv,0,0) for primer: " + primer + " is: " + tm);
		
		primer = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 0, 4E-03, 0.2E-03);
		System.out.println("T_m (pc,0,dv,dntp) for primer: " + primer + " is: " + tm);
	
		primer = "ATAT";
		tm = lucia.computeTM(primer, 50E-09);
		System.out.println("T_m for primer: " + primer + " is: " + tm);
		tm = lucia.computeCationCorrectedTM(primer, 50E-09, 50E-03);
		System.out.println("T_m for primer: " + primer + " is: " + tm);
		
		tm = lucia.computeCationCorrectedTM("ATGC", 50E-09, 50E-03);
		System.out.println("T_m for primer: " + "ATGC" + " is: " + tm);
		
		primer = "TAATACGACTCACTATAGGG";
		tm = lucia.computeCationCorrectedTM(primer, 50E-09, 50E-03);
		System.out.println("T_m for primer: " + primer + " is: " + tm);
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, 50E-09, 50E-03, 1.5E-03, 0.2E-03);
		System.out.println("T_m for primer: " + primer + " is: " + tm);
		
		PrimerSearchParameters params = new PrimerSearchParameters();
		primer = "AGAGAGAGAGAGAGAGAGAGA";
		primer = "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
		tm = lucia.computeMonoAndDivalentCationCorrectedTM(primer, params.getPRIMER_CONCENTRATION(), params.getMONOVALENT_CATION_CONCENTRATION(), params.getDIVALENT_CATION_CONCENTRATION(), params.getDNTP_CONCENTRATION());
		System.out.println("T_m for primer: " + primer + " is: " + tm);
	}

	/* (non-Javadoc)
	 * @see primerDesign.algo.PrimerMeltingTempCalculation#computeTM(java.lang.String)
	 */
	public double computeTM(String primer) {
		return this.computeMonoAndDivalentCationCorrectedTM(primer, SEARCH_PARAMS.getPRIMER_CONCENTRATION(), SEARCH_PARAMS.getMONOVALENT_CATION_CONCENTRATION(), SEARCH_PARAMS.getDIVALENT_CATION_CONCENTRATION(), SEARCH_PARAMS.getDNTP_CONCENTRATION());
		//return this.computeCationCorrectedTM(primer, Constants.PRIMER_CONCENTRATION, Constants.MONOVALENT_CATION_CONCENTRATION);
	}
}
