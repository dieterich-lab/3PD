/**
 * 
 */
package primerDesign.dsc.indexStructures.primerMisprimingCheck;

import java.io.IOException;
import java.io.Serializable;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.dsc.Primer;
import primerDesign.dsc.PrimerTypes;
import primerDesign.dsc.indexStructures.blat.BLATQueryClient;
import primerDesign.dsc.indexStructures.blat.MisprimingVerifier;
import primerDesign.dsc.indexStructures.blat.RestrictionSitesIndex;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a mispriming scan for a primer in a BLAT-based index.
 * 
 * @author Sebastian Fršhler
 *
 */
public class Blat3CPrimerMisprimingScan implements PrimerMisprimingCheck, Serializable{
	
	private static final long serialVersionUID = 1L;
	private BLATQueryClient client;
	private RestrictionSitesIndex index;
	PrimerSearchParameters params;
	private int scanLength;

	/**
	 * Initializes a blat primer mispriming scan for 3C primers.
	 * 
	 * @param client the wrapper class to the blat index
	 * @param index the index of restriction sites for all organisms and enzymes of interest
	 * @param params the primer search parameters
	 */
	public Blat3CPrimerMisprimingScan(BLATQueryClient client, RestrictionSitesIndex index, PrimerSearchParameters params){
		this.client = client;
		this.index = index;
		this.params = params;
		this.scanLength = params.getPRIMER_END_MISMATCH_SCAN_LENGTH();
	}
	
	/**
	 * Scans for misprimings of 3C primer 'primer' in a given background sequence (in forward and reverse direction).
	 * 
	 * Misprimings of 'primer' are only considered iff the last x basepairs of primer sequence are included in the mispriming
	 * and iff the mispriming is sufficiently close to a restriction site used in the 3C assay!
	 * The first mispriming which is 'sufficiently close' to a restriction site is considered to be the REAL priming of 'primer'!
	 */
	@Override
	public boolean hasMisprimings(Primer primer){
		boolean result = false;
		try{
			ObjectArrayList hits = this.client.getPrimings("Primer", primer.getSubsequence(Math.max(0, primer.getSequenceLength() - this.scanLength)));
			if(hits.size() > 1){
				result = MisprimingVerifier.hasUnsafeMisprimings(hits, this.index, this.params);
			}
		}catch(IOException e){
			e.printStackTrace();
			//System.exit(1);
		}catch(InterruptedException e){
			e.printStackTrace();
			//System.exit(1);
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck#setPrimerSearchParams(primerDesign.util.PrimerSearchParameters)
	 */
	@Override
	public void setPrimerSearchParams(PrimerSearchParameters prams) {
		this.params = prams;
	}	
	
	/**
	 * Sets the restriction site index.
	 * 
	 * @param index the restriction site index to be set
	 */
	public void setIndex(RestrictionSitesIndex index){
		this.index = index;
	}
	
	/**
	 * Sets the BLAT client.
	 * 
	 * @param client the BLAT client
	 */ 
	public void setClient(BLATQueryClient client){
		this.client = client;
	}
	
	public static void main(String[] args) throws IllegalAlphabetException, IllegalSymbolException{
		int count = Integer.parseInt(args[1]);
		String host = "localhost";
		String port = "10110";
		System.out.print("Running benchmark, querying index: " + args[0] + "(host: " + host + ", port: " + port + ") with: " + count + " random primers.");
		PrimerSearchParameters params = new PrimerSearchParameters();
		Blat3CPrimerMisprimingScan scan = new Blat3CPrimerMisprimingScan(new BLATQueryClient("/Users/froehler/Downloads/BLAT-bin/gfClient", host, port, "/Users/froehler/Downloads/BLAT-bin"), RestrictionSitesIndex.deserialize(args[0]), params);
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		params.setEnzyme(enzyme);
		SimpleTimer timer = new SimpleTimer();
		for(int i=0; i<count; i++){
			scan.hasMisprimings(new Primer(SeqTools.getRandomPrimerSequence(18, 30), PrimerTypes.forwardPrimer, 0, 0, 0, params));
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
