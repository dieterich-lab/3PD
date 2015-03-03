/**
 * 
 */
package primerDesign.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;


/**
 * A 3C restriction site screener class.
 * 
 * @author Sebastian Fršhler
 *
 */
public class CCCRSSScreener {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws IOException, NumberFormatException, IllegalAlphabetException, IllegalSymbolException {
		SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
		ArrayList<SimpleContigImpl> contigs = new ArrayList<SimpleContigImpl>();
		ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File(args[1]));
		
		while(parser.hasNextContig()){
			contigs.add(parser.parseNextContigIgnoreCase());
		}
		
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
		
		for(SimpleContig contig : contigs){
			for(int i=0; i<enzymes.size(); i++){
					RestrictionEnzyme enzyme = enzymes.get(i);
				RestrictionEnzymeManager.register(enzyme, new TreeSet());
				mapper.clearEnzymes();
				mapper.addEnzyme(enzyme);
				
				Sequence dnaSeq = null;
				
				try{
					dnaSeq = DNATools.createDNASequence(new String(contig.getSequence()), "");
				}catch(IllegalSymbolException e){
					e.printStackTrace();
				}
				Sequence newSeq = mapper.annotate(dnaSeq);
				Iterator iter = newSeq.features();
				
				if(iter.hasNext()){
					enzymes.remove(i);
				}
			}
		}
		
		System.out.println("The following enzymes have NO recognition sites within the contigs provided:");
		System.out.println("----------------------------------------------------------------------------");
		for(RestrictionEnzyme enzyme : enzymes){
			System.out.println(enzyme.toString());
		}
	}
}
