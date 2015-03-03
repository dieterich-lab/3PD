package primerDesign.util;

import primerDesign.algo.PrimerMeltingTempCalculation;
import primerDesign.algo.SantaLuciaTM;

/**
 * Class for constants used in this program
 * 
 * @author Sebastian Fršhler
 *
 */
public class Constants {
	public static final PrimerMeltingTempCalculation PRIMER_TM_CALC_METHOD = new SantaLuciaTM();
	
	// general constants
	public static final int MAX_PRIMER_LENGTH_FOR_TM_CALC = 36;  // max primer length for "reliable" tm calculation
	public static final double ABSOLUTE_ZERO_TEMP = -273.15;  // the absolute point zero temperature
	public static final double GAS_CONST_R = 1.987;  // the universal gas constant
	
	public static final String REPETITIVE_ELEMENT_CHARACTER = "N";
	
	// parallelization parameters
	public static int MAX_NUM_PICKING_THREADS = Runtime.getRuntime().availableProcessors();
	public static final int MAX_NUM_RESTRICTION_MAPPER_THREADS = 1; // There appears to be a bug in the biojava 1.5 restriction mapper when using multiple processors/threads

	// optimizations
	public static boolean doEarlyMMScan = true; // whether to scan for misprimings at primer enumeration instead of at primer pair set evaluation
	public static boolean doSortedDOPTScreen = false; // whether to sort pairs matrix by dist to opt pair prior to screening - screen for homonegity at the cost of alignment values

	// screening variants of matrix (greedy = sorted matrix column, increasing # primer pairs)
	public static boolean doSortedSitesScreening = true;
	public static boolean doGreedySitesScreening = false;
}
