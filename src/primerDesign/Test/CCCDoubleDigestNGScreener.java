/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
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

import primerDesign.util.Constants;
import primerDesign.util.RestrictionEnzymeListParser;
import primerDesign.util.SimpleContig;
import primerDesign.util.SlimFastaParser;

/**
 * @author froehler
 *
 */
public class CCCDoubleDigestNGScreener {
	private static int MIN_FRAGMENT_LENGTH;
	private static int MAX_FRAGMENT_LENGTH;
	private static boolean printTraceLog = false;

	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 * @throws NumberFormatException 
	 */
	public static void main(String[] args) throws NumberFormatException, IllegalAlphabetException, IllegalSymbolException, IOException {
		String enzymesList = args[0];
		String sequence = args[1];
		MIN_FRAGMENT_LENGTH = Integer.parseInt(args[2]);
		MAX_FRAGMENT_LENGTH = Integer.parseInt(args[3]);
		
		HashMap<String, Integer> pairs = new HashMap<String, Integer>();
		
		// read in enzymes list'enzymesList'
		ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File(enzymesList));
		
		// init counters for all enzymes combinations
		for(int i=0; i<enzymes.size(); i++){
			for(int j=0; j<enzymes.size(); j++){
				if(i != j) pairs.put(enzymes.get(i) + ":" + enzymes.get(j), 0);
			}
		}
				
		RestrictionMapper mapperA = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
		RestrictionMapper mapperB = new RestrictionMapper(new SimpleThreadPool(Constants.MAX_NUM_RESTRICTION_MAPPER_THREADS, true));
		Sequence dnaSeqA = null;
		Sequence newSeqA;
		Sequence dnaSeqB = null;
		Sequence newSeqB;
		int positionA;
		int positionB;
		Iterator<Feature> iterA;
		Iterator<Feature> iterB;
		String seq;
		
		// read in sequence'sequence'
		SlimFastaParser parser = new SlimFastaParser(new File(sequence));
		//foreach contig of sequence
		SimpleContig currentContig;
		long genomeSize = 0;
		int maxUpstreamDiff;
		int maxDownstreamDiff;
		int closestUpstreamSite = -1;
		int closestDownstreamSite = -1;
		while(parser.hasNextContig()){
			currentContig = parser.parseNextContig();
			genomeSize += currentContig.getSequenceLength();
			if(printTraceLog) System.err.println("Screening Contig: " + currentContig.getID());
			seq = new String(currentContig.getSequence());
			// foreach enzyme 1
			RestrictionEnzyme a;
			RestrictionEnzyme b;
			for(int i=0; i<enzymes.size(); i++){
				a = enzymes.get(i);
				if(printTraceLog) System.err.println("\tEnzyme A: " + a.getName());
				// scan contig sequence for RSSs of enzyme 1 
				RestrictionEnzymeManager.register(a, new TreeSet());
				mapperA.clearEnzymes();
				mapperA.addEnzyme(a);
				try{
					dnaSeqA = DNATools.createDNASequence(seq, "");
				}catch(IllegalSymbolException e){
					e.printStackTrace();
				}
				newSeqA = mapperA.annotate(dnaSeqA);
				
				// foreach RSS of enzyme 1: scan vicinity of RSS for RSSs for enzyme 2
				for(int j=0; j<enzymes.size(); j++){
					if(i == j) continue;
					else{
						b = enzymes.get(j);
						if(printTraceLog) System.err.println("\tEnzyme B: " + b.getName());
						RestrictionEnzymeManager.register(b, new TreeSet());
						mapperB.clearEnzymes();
						mapperB.addEnzyme(b);
						iterA = newSeqA.features();
						// scan each restriction site for enzyme A for associated restriction sites for enzyme B
						while(iterA.hasNext()){
							positionA = iterA.next().getLocation().getMin();
							try{
								if(positionA - MAX_FRAGMENT_LENGTH/2 >= 0 && positionA + MAX_FRAGMENT_LENGTH/2 < currentContig.getSequenceLength()){
									dnaSeqB = DNATools.createDNASequence(seq.substring(positionA - MAX_FRAGMENT_LENGTH/2, positionA + MAX_FRAGMENT_LENGTH/2), "");
								}
								else continue;
							}catch(IllegalSymbolException e){
								e.printStackTrace();
							}
							newSeqB = mapperB.annotate(dnaSeqB);
							iterB = newSeqB.features();
							
							// get closest restriction sites for enzyme B
							maxUpstreamDiff = Integer.MAX_VALUE;
							maxDownstreamDiff = Integer.MAX_VALUE;
							while(iterB.hasNext()){
								// get closests RSS for enzyme 2 - check whether distance constraints are satisfied, iff yes: increment count of RE combination
								positionB = iterB.next().getLocation().getMin();
								if(positionB < positionA && positionA - positionB < maxUpstreamDiff){
									maxUpstreamDiff = positionA - positionB;
									closestUpstreamSite = positionB;
								}
								else if(positionB > positionA && positionB - positionA < maxDownstreamDiff){
									maxDownstreamDiff = positionB - positionA;
									closestDownstreamSite = positionB;
								}
							}
							
							// scan for proper distance constraints of sites for enzyme B
							// if there is at least one RSS for enzyme 2 in acceptable distance interval on both sides: accept site
							if(maxUpstreamDiff >= MIN_FRAGMENT_LENGTH/2 && maxDownstreamDiff >= MIN_FRAGMENT_LENGTH/2 && maxUpstreamDiff <= MAX_FRAGMENT_LENGTH/2 && maxDownstreamDiff <= MAX_FRAGMENT_LENGTH/2){
								pairs.put(a + ":" + b, pairs.get(a + ":" + b) + 1);
								if(printTraceLog) System.err.println(a + " - " + b + " Contig: " + currentContig.getID() + " 2: " + closestUpstreamSite + " 1: " + positionA + " 2': " + closestDownstreamSite);
							}
						}
					}
				}
			}
		}
		
		// output enzyme pair with highest number of hits
		String best= null;
		String candidate;
		int max = 0;
		Set<String> keys = pairs.keySet();
		Iterator<String> iter = keys.iterator();
		while(iter.hasNext()){
			candidate = iter.next();
			if(pairs.get(candidate) > max){
				max = pairs.get(candidate);
				best = candidate;
			}
		}
		NumberFormat format = NumberFormat.getInstance();
		String[] temp = null;
		if(best != null){
			 temp = best.split(":");
			 System.out.println("The best enzyme pair is:");
				System.out.println("\t"+ temp[0] + " - " + temp[1]);
			System.out.println("This pair has: " + format.format(pairs.get(best)) + " detection sites corresponding to an average resolution of: 1 site per " + format.format(genomeSize / pairs.get(best)) + " bp");
		}
		else{
			System.out.println("No best enzyme combination can be found - check your setup!");
		}
	}
	
//	class DoubleDigestStat{
//		private double numSites;
//		private ArrayList<Site> sites;
//		
//		public Site[] getSortedSites(){
//			Site[] result = new Site[this.sites.size()];
//			result = this.sites.toArray(result);
//			Arrays.sort(result);
//			return result;
//		}
//	}
	
	class Site implements Comparable{
		private String contig;
		private int position;

		public int compareTo(Object o) {
			Site other = (Site) o;
			if(this.contig.compareTo(other.contig) != 0) return this.contig.compareTo(other.contig);
			else return (this.position < other.position) ? -1 : (this.position > other.position) ? +1 : 0;
		}
	}
}
