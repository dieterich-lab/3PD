/**
 * 
 */
package primerDesign.dsc.indexStructures.primerMisprimingCheck;

import java.io.IOException;

import primerDesign.algo.PrimerMispriming;
import primerDesign.dsc.Primer;
import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.dsc.indexStructures.blat.BLATQueryClient;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleSlimContig;
import primerDesign.util.SlimGenome;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a mispriming scan for a primer in a BLAT-based index using sequence information instead of a restriction site index.
 * 
 * @author Sebastian Fr�hler
 *
 */
public class Blat3CPrimerMisprimingScanPlusSeq implements PrimerMisprimingCheck{
	private BLATQueryClient client;
	private String slimGenomeFilename;
	private SlimGenome genomeSeq;
	private PrimerSearchParameters params;
	private int scanLength;
	private boolean debug = true;
	
	/**
	 * Initilizes a mispriming scan.
	 * 
	 * @param client the BLAT client ot be used
	 * @param slimGenomeFilename the name of the slim serialized reference genome
	 * @param params the primer search parameters for 3PD
	 */
	public Blat3CPrimerMisprimingScanPlusSeq(BLATQueryClient client, String slimGenomeFilename, PrimerSearchParameters params){
		this.client = client;
		this.slimGenomeFilename = slimGenomeFilename;
		this.params = params;
		this.scanLength = params.getPRIMER_END_MISMATCH_SCAN_LENGTH();
		if(this.debug) System.err.println("DEBUG: init Blat3CPrimerMisprimingScanPlusSeq with params: " + this.client.toString() + " " + this.slimGenomeFilename + " " + params);
	}
	
	/**
	 * Loads a slim serialized genome sequence into memory.
	 * 
	 * @param filename the filename of the slim serialized genome to be loaded
	 */
	public void loadSlimGenome(String filename){
		this.genomeSeq = SlimGenome.deserialize(filename);
		if(this.debug) System.err.println("DEBUG: Deserializing slim genome from file: " + filename);
	}
	
	/**
	 * Scans for misprimings of 3C primer 'primer' in a given background sequence (in forward and reverse direction).
	 * 
	 * Misprimings of 'primer' are only considered iff the last x basepairs of primer sequence are included in the mispriming
	 * and iff the mispriming is sufficiently close to a restriction site used in the 3C assay!
	 * The first mispriming which is 'sufficiently close' to a restriction site is considered to be the REAL priming of 'primer'!
	 * 
	 * @param primer the primer to check
	 * 
	 * @return true iff the primer has misprimings
	 */
	@Override
	public boolean hasMisprimings(Primer primer){
		// load genomic sequence iff not loaded
		if(this.genomeSeq == null) loadSlimGenome(this.slimGenomeFilename);
		if(this.debug) System.err.println("DEBUG: exec 'hasMisprimings', genomeSeq == null?: " + (this.genomeSeq == null));
		
		boolean result = false;
		try{
			// get primings
			ObjectArrayList hits = this.client.getPrimings("Primer", primer.getSubsequence(Math.max(0, primer.getSequenceLength() - this.scanLength)));
			// post-process primings, add ref to contig sequence
			SimpleSlimContig contig;
			for(int i=0; i<hits.size(); i++){
				contig = (SimpleSlimContig)((IndexHitImpl)hits.getQuick(i)).getContig();
				contig.setSlimSequence(((SimpleSlimContig)this.genomeSeq.getContig(contig.getID())).getSlimSequence());
			}
			
			// screen if primings are unsafe
			if(hits.size() <= 1) return false;
			else{
				if(PrimerMispriming.hasSafeDistanceToNextRSS(primer, this.params.getEnzyme(), hits, this.params)) return false;
				else return true;
			}
		}catch(IOException e){
			e.printStackTrace();
			//System.exit(1);
		}catch(InterruptedException e){
			e.printStackTrace();
			//System.exit(1);
		}
		catch(NullPointerException e){
			e.printStackTrace();
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
}
