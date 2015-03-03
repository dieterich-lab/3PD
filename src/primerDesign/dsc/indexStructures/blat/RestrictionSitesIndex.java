/**
 * 
 */
package primerDesign.dsc.indexStructures.blat;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.molbio.RestrictionSite;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.SimpleThreadPool;

import primerDesign.util.Constants;
import primerDesign.util.FileTools;
import primerDesign.util.RestrictionEnzymeListParser;
import primerDesign.util.SimpleContig;
import primerDesign.util.SlimFastaParser;

/**
 * Wraps an index of restriction site positions within a genome.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RestrictionSitesIndex implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3998532871429223030L;
	private HashMap<String, int[]> index;
	private HashSet<String> organisms;
	private HashSet<String> enzymes;
	private HashSet<String> contigs;
	
	private static final String SERIALIZATION_SUFFIX = ".RSSIndex.ser";
	
	public RestrictionSitesIndex(){
		this.index = new HashMap<String, int[]>();
		this.organisms = new HashSet<String>();
		this.enzymes = new HashSet<String>();
		this.contigs = new HashSet<String>();
	}
	
	/**
	 * Adds a list of restriction sites to the index.
	 * 
	 * @param enzyme the restriction enzyme
	 * @param organism the organism
	 * @param chromosome the chromosome of the organism
	 * @param sites the list of restriction sites for enzyme 'enzyme' in organism 'organism' on chromosome 'chromosome'
	 */
	public void addSites(String enzyme, String organism, String chromosome, int[] sites){
		this.index.put(enzyme + organism + chromosome, sites);
		this.enzymes.add(enzyme);
		this.organisms.add(organism);
		this.contigs.add(chromosome);
	}
	
	/**
	 * Returns the list of restriction sites for the given parameters or NULL of no such list can be found.
	 * 
	 * @param enzyme the restriction enzyme
	 * @param Organism the organism
	 * @param chromosome the chromosome of the organism
	 * @return the list of restriction sites for the given parameters or NULL of no such list can be found
	 */
	public int[] getSites(String enzyme, String Organism, String chromosome){
		return this.index.get(enzyme + Organism + chromosome);
	}
	
	/**
	 * Tests whether a genomic position is sufficiently close to a restriction site w.r.t the threshold specified.
	 * 
	 * @param position the position to test
	 * @param threshold the threshold defining what 'sufficiently close' means
	 * @param enzyme the enzyme generating the restriction site
	 * @param organism the organism to scan
	 * @param chromosome the chromosome of the organism
	 * 
	 * @return true iff 'position' has distance <= 'threshold' to the next restriction site of enzyme 'enzyme' in organism 'organism'
	 */
	public boolean isSufficientlyClose(int position, int threshold, String enzyme, String organism, String chromosome){
		chromosome = chromosome.replace('_', ' ').substring(0, chromosome.length() - 1);
		int[] positions = getSites(enzyme, organism, chromosome);
		if(positions == null){
			System.err.println("Contains Contig?: " + this.containsContig(chromosome));
			System.err.println("Contains Enzyme?: " + this.containsEnzyme(enzyme));
			System.err.println("Contains Organism?: " + this.containsOrganism(organism));
			throw new IllegalStateException("No positions were specified for organism: " + organism + " and chromosome: " + chromosome + " and enzyme: " + enzyme);
		}
		
		for(int pos : positions){
			if(Math.abs(position - pos) <= threshold) return true;
		}
		return false;
	}
	
	public boolean containsOrganism(String organism){
		return this.organisms.contains(organism);
	}
	
	public boolean containsEnzyme(String enzyme){
		return this.enzymes.contains(enzyme);
	}
	
	public boolean containsContig(String contig){
		return this.contigs.contains(contig);
	}
	
	public String[] getOrganisms(){
		String[] result = new String[this.organisms.size()];
		Iterator<String> iter = this.organisms.iterator();
		int i=0;
		while(iter.hasNext()){
			result[i++] = iter.next();
		}
		return result;
	}
	
	/**
	 * Serializes the current index to file 'filename'.
	 * 
	 * @param filename the filename to serialize the index to
	 */
	public final void serialize(File filename) {
		try{
			//ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename + SERIALIZATION_SUFFIX));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Deserializes an index of this datatype from file 'filename'.
	 * 
	 * @param filename the filename to deserialize the index from
	 * 
	 * @return an index of this datatype
	 */
	public static RestrictionSitesIndex deserialize(String filename){
		RestrictionSitesIndex result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (RestrictionSitesIndex) in.readObject();
//			System.out.println("### Org: " + result.organisms.contains("Ppacificus"));
//			System.out.println("### Contig: " + result.contigs.contains("SUPERCONTIG:1:CONTIG40:1:1108810:1_SUPERCONTIG_CONTIG40_"));
//			System.out.println("### Enzyme: " + result.enzymes.contains("EcoRI"));
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws NumberFormatException, IllegalAlphabetException, IllegalSymbolException, IOException{
		if(args.length < 2){
			System.out.println(getUsage());
			System.exit(1);
		}
		
		ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File(args[1]));
		if(enzymes.size() == 0) throw new IllegalArgumentException("No restriction enzymes were parsed from file!\n" + getUsage());
		
		RestrictionSitesIndex index = new RestrictionSitesIndex();
		
		SlimFastaParser fastaParser;
		SimpleContig contig;
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
		Sequence dnaSeq = null;
		Sequence newSeq = null;
		String organismName;
		
		// process each sequence
		for(int i=2; i<args.length; i++){
			fastaParser = new SlimFastaParser(new File(args[i]));
			organismName = FileTools.extractFilename(new File(args[i])).replaceAll(".fa", "");
			// process each chromosome
			while(fastaParser.hasNextContig()){
				contig = fastaParser.parseNextContigIgnoreCase();
				// process each
				for(RestrictionEnzyme enzyme : enzymes){
					RestrictionEnzymeManager.register(enzyme, new TreeSet());
					mapper.clearEnzymes();
					mapper.addEnzyme(enzyme);
					
					try{
						dnaSeq = DNATools.createDNASequence(new String(contig.getSequence()), "");
					}catch(IllegalSymbolException e){
						e.printStackTrace();
					}
					newSeq = mapper.annotate(dnaSeq);
					Iterator<RestrictionSite> iter = newSeq.features();
					
					int[] rss = new int[newSeq.countFeatures()];
					int temp = 0;
					while(iter.hasNext()){
						rss[temp++] = iter.next().getPosition();
					}
					index.addSites(enzyme.getName(), organismName, contig.getID(), rss);
				}
			}
		}
		if(index.index.size() > 0) index.serialize(new File(args[0]));
		else throw new IllegalArgumentException(getUsage());
	}
	
	private static String getUsage(){
		return "Usage: RestrictionSitesIndex <outfile> <restriction enzymes list> <Organism1> ...";
	}
}
