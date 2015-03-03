/**
 * 
 */
package primerDesign.dsc.indexStructures.primerMisprimingCheck;

import java.io.Serializable;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.algo.PrimerMispriming;
import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.DNASequenceIndexDeserializer;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a mispriming scan using an enhanced suffix array index structure.
 * 
 * @author Sebastian Fršhler
 *
 */
public class ESA3CPrimerMisprimingScan implements PrimerMisprimingCheck, Serializable {

	private static final long serialVersionUID = 1L;
	private DNASequenceIndex index;
	private PrimerSearchParameters params;
	private RestrictionEnzyme enzyme;
	
	/**
	 * Initializes a new mispriming scan.
	 * 
	 * @param index the enhanced suffix array to be used
	 * @param params the primer search parameters for 3PD
	 */
	public ESA3CPrimerMisprimingScan(DNASequenceIndex index, PrimerSearchParameters params){
		this.index = index;
		this.params = params;
		this.enzyme = this.params.getEnzyme();
	}

	/**
	 * Scans for misprimings of 3C primer's end 'primer' in a given background sequence (in forward and reverse direction).
	 * 
	 * The length of the primer's end is defined in the primer search parameters 'params'.
	 * Misprimings of 'primer' are only considered iff the last x basepairs of primer sequence are included in the mispriming,
	 * iff the mispriming is sufficiently close to a restriction site used in the 3C assay and iff the melting temperature is 
	 * 'sufficiently close' (as defined in 'params') to the primers melting temperature!
	 * The first mispriming which is 'sufficiently close' to a restriction site is considered to be the REAL priming of 'primer'!
	 * 
	 * @param primer the primer to check
	 * 
	 * @return true iff the primer has misprimings
	 */
	@Override
	public boolean hasMisprimings(Primer primer) {
		ObjectArrayList hits = this.index.findHitPositions(primer.getSequence().substring(primer.getLength()  - this.params.getPRIMER_END_MISMATCH_SCAN_LENGTH()));
		if(hits.size() <= 1) return false;
		else{
			if(PrimerMispriming.hasSafeDistanceToNextRSS(primer, this.enzyme, hits, this.params)) return false;
			else return true;
		}
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck#setPrimerSearchParams(primerDesign.util.PrimerSearchParameters)
	 */
	@Override
	public void setPrimerSearchParams(PrimerSearchParameters prams) {
		this.params = prams;
	}	
	
	/**
	 * Sets the enhanced suffix array index to be used.
	 * 
	 * @param index the enhanced suffix array index to be used
	 */
	public void setIndex(DNASequenceIndex index){
		this.index = index;
	}
	
	public static void main(String[] args) throws IllegalAlphabetException, IllegalSymbolException{
		int count = Integer.parseInt(args[1]);
		DNASequenceIndex index = DNASequenceIndexDeserializer.deserialize(args[0]);
		System.out.print("Running benchmark, querying index: " + args[0] + " with: " + count + " random primers.");
		PrimerSearchParameters params = new PrimerSearchParameters();
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		params.setEnzyme(enzyme);
		ESA3CPrimerMisprimingScan scanner = new ESA3CPrimerMisprimingScan(index, params);
		SimpleTimer timer = new SimpleTimer();
		for(int i=0; i<count; i++){
			scanner.hasMisprimings(new Primer(SeqTools.getRandomPrimerSequence(18, 30), PrimerTypes.forwardPrimer, 0, 0, 0, params));
		}
		System.out.println(" - done in " + timer.getTimeString());
		System.out.print("Generating " + count + " random primers itself took");
		timer.getTimeString();
		for(int i=0; i<count; i++){
			new Primer(SeqTools.getRandomPrimerSequence(18, 30), PrimerTypes.forwardPrimer, 0, 0, 0, params);
		}
		System.out.println(" - " + timer.getTimeString());
	}
}