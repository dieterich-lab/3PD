/**
 * 
 */
package primerDesign.testSuite.algo;

import junit.framework.TestCase;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.algo.GAPrimerPicker;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.RestrictionSite;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;

/**
 * @author froehler
 *
 */
public class GAPrimerPickerTest extends TestCase {

	/**
	 * Test method for {@link primerDesign.algo.GAPrimerPicker#pickBestPrimerSet(primerDesign.dsc.RestrictionSite[])}.
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 */
	public void testPickBestPrimerSet() throws IllegalAlphabetException, IllegalSymbolException {
		
		PrimerSearchParameters params = new PrimerSearchParameters();
		
		Primer bestPrimer = new Primer("GCGCGC", PrimerTypes.forwardPrimer, params);
		
		Primer[] upstreamFirst = new Primer[3];
		upstreamFirst[0] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		upstreamFirst[1] = bestPrimer;
		upstreamFirst[2] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		
		Primer[] downstreamFirst = new Primer[3];
		downstreamFirst[0] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		downstreamFirst[1] = bestPrimer;
		downstreamFirst[2] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		
		Primer[] taqManFirst = new Primer[1];
		taqManFirst[0] = bestPrimer;
		
		RestrictionSite first = new RestrictionSite(params);
		first.setValidUpstreamPrimers(upstreamFirst);
		first.setValidDownstreamPrimers(downstreamFirst);
		first.setValidTaqManProbes(taqManFirst);
		
		Primer[] upstreamSecond = new Primer[3];
		upstreamSecond[0] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		upstreamSecond[1] = bestPrimer;
		upstreamSecond[2] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		
		Primer[] downstreamSecond = new Primer[1];
		downstreamSecond[0] = bestPrimer;
		
		Primer[] taqManSecond = new Primer[2];
		taqManSecond[0] = new Primer(SeqTools.getRandomPrimerSequence(params.getMIN_PRIMER_LENGTH(), params.getMAX_PRIMER_LENGTH()), PrimerTypes.forwardPrimer, params);
		taqManSecond[1] = bestPrimer;
		
		RestrictionSite second = new RestrictionSite(params);
		second.setValidUpstreamPrimers(upstreamSecond);
		second.setValidDownstreamPrimers(downstreamSecond);
		second.setValidTaqManProbes(taqManSecond);
		
		RestrictionSite[] candidates = new RestrictionSite[2];
		candidates[0] = first;
		candidates[1] = second;
		
		GAPrimerPicker picker = new GAPrimerPicker();
		PrimerPair[] bestPairs = picker.pickBestPrimerSet(candidates);
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[0]).getForwardPrimer());
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[0]).getReversePrimer());
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[0]).getHybridizationProbe());
		
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[1]).getForwardPrimer());
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[1]).getReversePrimer());
		assertEquals(bestPrimer, ((PrimerPair)bestPairs[1]).getHybridizationProbe());
	}
}
