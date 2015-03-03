package primerDesign.dsc;

/**
 * Enumerates some primer types.
 * 
 * @author Sebastian Fršhler
 *
 */
public enum PrimerTypes {
	forwardPrimer, reversePrimer, hybridizationProbe;
	
	public static boolean isCompatible(Enum first, Enum second){
		if(first.equals(second)
				|| (first.equals(forwardPrimer) && second.equals(reversePrimer))
				|| (first.equals(reversePrimer) && second.equals(forwardPrimer))
				|| (first.equals(hybridizationProbe) && second.equals(hybridizationProbe))) return true;
		else return false;
	}
}
