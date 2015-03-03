/**
 * A simple candidates printer class.
 */
package primerDesign.web.tools;

import primerDesign.dsc.RestrictionSite;

/**
 * @author Sebastian Fršhler
 *
 */
public class PrimerProbesCandidatesPrinter {
	public String printCandidates(RestrictionSite[] optimalSites){
		StringBuffer buffy = new StringBuffer();
		for(int i=0; i< optimalSites.length;i++){
			buffy.append("Valid upstream primers\n");
			for(int j=0; j<optimalSites[i].getValidUpstreamPrimers().length; j++) buffy.append(optimalSites[i	].getUpstreamPrimer(j).toString() + "\n");
			buffy.append("Valid downstream primers\n");
			for(int j=0; j<optimalSites[i].getValidDownstreamPrimers().length; j++) buffy.append(optimalSites[i].getDownstreamPrimer(j).toString() + "\n");
			buffy.append("Valid hybridization probes primers\n");
			for(int j=0; j<optimalSites[i].getValidTaqManProbes().length; j++) buffy.append(optimalSites[i].getTaqManProbe(j).toString() + "\n");
		}
		return buffy.toString();
	}
}
