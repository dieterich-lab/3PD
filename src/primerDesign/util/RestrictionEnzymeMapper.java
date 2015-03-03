/**
 * 
 */
package primerDesign.util;

import java.util.Iterator;
import java.util.TreeSet;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;

/**
 * A restriction enzyme mapper class.
 * 
 * @author Sebastian Fršhler
 *
 */

public class RestrictionEnzymeMapper {
	/**
	 * Maps restriction sites onto a sequence
	 * 
	 * @param sequence the sequence to screen
	 * @param enzyme the restriction enzyme to be screened with
	 * 
	 * @return a list of positions of the enzyme in sequence
	 */
	public int[] mapFeature(String sequence, RestrictionEnzyme enzyme){
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(1, true));
		mapper.addEnzyme(enzyme);
		RestrictionEnzymeManager.register(enzyme, new TreeSet());
		
		int[] result = null;
		
		try{
			Sequence newSeq = mapper.annotate(DNATools.createDNASequence(sequence, ""));
			
			result = new int[newSeq.countFeatures()];
			int position = 0;
			Iterator iter = newSeq.features();
			while(iter.hasNext()){
				result[position++] = ((Feature)iter.next()).getLocation().getMin();
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IllegalAlphabetException, IllegalSymbolException{
		String sequence = "ATTTCAAAGATCAAAATGGTATTTTCAGGATCCTCGAAGGGCGGTGAGAAGCGACAACGAACTGCGTACACTCGGAATCAAGTATTGGAATTGGAGAAGGAATTCCACTTCAATAAATATTTGACGAGAAAGAGACGGATTGAGATATCGCATTCGTTGATGCTCAGTGAGAGACAGGTTAGCTCTTGTTATTATAAATA";
		RestrictionEnzymeMapper mapper = new RestrictionEnzymeMapper();
		int[] results = mapper.mapFeature(sequence, new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0));
		System.out.println(sequence.length());
		System.out.println(results.length);
		System.out.println(results[0]);
	}
}
