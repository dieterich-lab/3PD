/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.TreeSet;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;

import primerDesign.util.Constants;
import primerDesign.util.RestrictionEnzymeListParser;
import primerDesign.util.SimpleContig;
import primerDesign.util.SlimFastaParser;

/**
 * @author froehler
 *
 */
public class GenomeRSSScreener {
	private static RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
	private static Sequence dnaSeq;
	private static Sequence newSeq;
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws IOException, NumberFormatException, IllegalAlphabetException, IllegalSymbolException {
		NumberFormat format = NumberFormat.getInstance();
		
		ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File(args[1]));
		
		for(RestrictionEnzyme enzyme : enzymes){
			SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
			
			RestrictionEnzymeManager.register(enzyme, new TreeSet());
			mapper.clearEnzymes();
			mapper.addEnzyme(enzyme);
			
			int count = 0;
			int seqLength = 0;
			SimpleContig contig = null;
			
			while(parser.hasNextContig()){
				try{
					contig = parser.parseNextContigIgnoreCase();
					dnaSeq = DNATools.createDNASequence(new String(contig.getSequence()), "");
				}catch(IllegalSymbolException e){
					e.printStackTrace();
				}
				newSeq = mapper.annotate(dnaSeq);
				
				count += newSeq.countFeatures();
				seqLength += contig.getSequenceLength();
			}
			//System.out.println("There are " + format.format(count) + " RSSs for enzyme " + enzyme.getName() + " in sequence file " + args[0] + " (1/" + format.format(seqLength/count) + " bp on average)");
			System.out.println(args[0] + "\t" + enzyme.getName() + "\t" + format.format(count) + "\t" + "1/" + format.format(seqLength/count));
		}
	}
}
