/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;
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

import primerDesign.util.SimpleContig;
import primerDesign.util.SlimFastaParser;

/**
 * @author froehler
 *
 */
public class RestrictionSitesEnumerator {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws IOException, NumberFormatException, IllegalAlphabetException, IllegalSymbolException {
		SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
		SimpleContig contig = parser.parseNextContigIgnoreCase();
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme(args[1], DNATools.createDNA(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(1, true));
		mapper.addEnzyme(enzyme);
		RestrictionEnzymeManager.register(enzyme, new TreeSet());
		
		Sequence dnaSeq = DNATools.createDNASequence(new String(contig.getSequence()), "");
		Sequence newSeq = mapper.annotate(dnaSeq);
		
		Iterator iter = newSeq.features();
		while(iter.hasNext()){
			System.out.println(((Feature)iter.next()).getLocation().getMin());
		}

	}

}
