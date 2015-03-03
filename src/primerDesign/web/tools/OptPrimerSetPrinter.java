/**
 * A simple optimal primer set printer class.
 */
package primerDesign.web.tools;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerPairSet;

/**
 * @author Sebastian Fršhler
 *
 */
public class OptPrimerSetPrinter {
	public String printBestPrimerSet(PrimerPairSet bestPrimerPairSet){
		StringBuffer buffer = new StringBuffer();
		buffer.append(String.format("%10s\t%s%s", "Type", Primer.toFormattedStringDescription(), "\n"));
		for(int i=0; i<bestPrimerPairSet.getNumPrimerPairs(); i++){
			buffer.append(((PrimerPair)bestPrimerPairSet.getPrimerPair(i)).toFormattedString() + "\n\n");
		}
		buffer.append("avgdOptPrimerPAir: " + bestPrimerPairSet.getAvgDistOptPrimerPair() + "\n");
		buffer.append("Homogenity score: " + bestPrimerPairSet.getHomogenityScore() + "\n");
		buffer.append("worst primer -  primer PA/PEA score (incl probes): " + bestPrimerPairSet.getMaxPairAlignScore() + " " + bestPrimerPairSet.getMaxPairAlignEndScore() + "\n");
		return buffer.toString();
	}
}
