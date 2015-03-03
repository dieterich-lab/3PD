/**
 * 
 */
package primerDesign.algo;

/**
 * Specifies a generic primer melting temperature calculation method.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface PrimerMeltingTempCalculation {

	/**
	 * This method specifies the melting temperature calculation for a dna oligonucleotide.
	 * 
	 * @param primer the primer to calculate the melting temperature for
	 * 
	 * @return the melting temperature of primer 'primer'
	 */
	public abstract double computeTM(String primer);
	
	/**
	 * This method specifies the melting temperature calculation for a dna oligonucleotide using concentration corrections.
	 * 
	 * @param primr the primer to calculate the melting temperature for
	 * @param primerConcentration the concentration of the primer
	 * @param cationConcentration the monovalent cation concentration - usually Na+
	 * @param divalentCationConcentration the divalent cation concentration - usually Mg2+
	 * @param dNTPConcentration the dNTP concentration
	 * 
	 * @return the melting temperature of primer 'primer'
	 */
	public abstract double computeMonoAndDivalentCationCorrectedTM(String primr, double primerConcentration, double cationConcentration, double divalentCationConcentration, double dNTPConcentration);

}