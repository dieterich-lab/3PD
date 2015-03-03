/**
 * 
 */
package primerDesign.util;

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

/**
 * Class for extracting sequence regions for a specific restriction enzyme.
 * 
 * @author Sebastan Fršhler
 *
 */
public class ExtractSequenceRegions {

	public static void main(String[] args) throws IOException, IllegalAlphabetException, IllegalSymbolException {
		int seqs = Integer.MAX_VALUE; //100;
		int position;
		int length = 100;
		
		SlimFastaParser parser = new SlimFastaParser(new File("/Users/froehler/SEQUENCE_INDICES/seq.fa"));
		SimpleContigImpl contig;
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
		
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(1, true));
		mapper.addEnzyme(enzyme);
		RestrictionEnzymeManager.register(enzyme, new TreeSet());
		while(parser.hasNextContig() && seqs > 0){
			contig = parser.parseNextContigIgnoreCase();
			try{
				Sequence dnaSeq = DNATools.createDNASequence(new String(contig.getSequence()), contig.getID());
				Sequence newSeq = mapper.annotate(dnaSeq);
				
				Iterator<Sequence> iter = newSeq.features();
				while(iter.hasNext() && seqs > 0){
					position = ((Feature)iter.next()).getLocation().getMin();
					if(position - length >= 0){
						System.out.println(">" + contig.getID() + "_" + seqs);
						System.out.println(contig.getSubsequence(position-length-5, position+length+5));
						seqs--;
					}
				}
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
	}
}
