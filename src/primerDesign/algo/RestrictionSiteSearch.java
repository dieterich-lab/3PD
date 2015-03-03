package primerDesign.algo;

import java.util.Iterator;
import java.util.TreeSet;
import java.util.Vector;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.utils.SimpleThreadPool;

import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.SequenceRegion;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.MyExtendedMath;
import primerDesign.util.NoRestrictionSitesFoundException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import weka.core.FastVector;
import cern.colt.list.ObjectArrayList;

/**
 * This class implements the search for a specific number of best 
 * (most homogenuously distributed) restriction sites in a give input sequence.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RestrictionSiteSearch {
	
	private Vector<RestrictionSite> allRestrictionSites = new Vector<RestrictionSite>();

	/**
	 * Computes the position of the most-centred restriction site in an input sequence.
	 * 
	 * @param sequence the sequence to scan for restriction sites
	 * @param enzyme the restriction enzyme 'generating' the restriction site
	 * @param sequenceRegion the sequence region this restriction site is located on
	 * @param offset the offset for computing absolute sequence positions
	 * @param searchParams the 3PD search parameters
	 * 	 
	 * @return the position of the most-centred restriction site
	 */	
	private RestrictionSite findBestRestrictionSite(String sequence, RestrictionEnzyme enzyme, SequenceRegion sequenceRegion, int offset, PrimerSearchParameters searchParams){
			
			RestrictionSite bestRSS = null;
			int intervalMean = MyExtendedMath.round(offset + sequence.length()/2);
			int bestDistance = Integer.MAX_VALUE;
			
			RestrictionSite[] allRestrictionSites = this.findAllRestrictionSites(sequence, enzyme, sequenceRegion, offset, searchParams);

			RestrictionSite currentRSS;
			int currentRSSPos;
			int currentDistance;
			for(int i=0; i< allRestrictionSites.length; i++){
				// remember position of best RSS
				currentRSS = (RestrictionSite) allRestrictionSites[i];
				currentRSSPos = currentRSS.getPosition();
				currentDistance = Math.abs(currentRSSPos - intervalMean);
				if(currentDistance < bestDistance){
					bestRSS = currentRSS; 
					bestDistance = currentDistance;
				}
			}
			return bestRSS;
		}
	
//	public RestrictionSite getNextClosestRestrictionSite(int position){
//		RestrictionSite site;
//	}
	
	/**
	 * Computes the position of all restriction site in a sequence.
	 * 
	 * @param sequence the sequence to scan
	 * @param enzyme the restriction enzyme generating the restriction sites
	 * @param sequenceRegion the sequence region this restriction site is located on
	 * @param offset the offset for computing absolute sequence positions
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return all restriction sites of enzyme 'enzyme' in sequence 'sequence
	 */
	private RestrictionSite[] findAllRestrictionSites(String sequence, RestrictionEnzyme enzyme, SequenceRegion sequenceRegion, int offset, PrimerSearchParameters searchParams){
		RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(1, true));
		mapper.addEnzyme(enzyme);
		RestrictionEnzymeManager.register(enzyme, new TreeSet());
		ObjectArrayList restrictionSites = new ObjectArrayList();
		//Vector<RestrictionSite> restrictionSites = new Vector<RestrictionSite>();
		
		try{
			Sequence dnaSeq = DNATools.createDNASequence(sequence, "");
			Sequence newSeq = mapper.annotate(dnaSeq);
			
			Iterator iter = newSeq.features();
			while(iter.hasNext()){
				RestrictionSite site = new RestrictionSite(((Feature)iter.next()).getLocation().getMin() + offset, enzyme, searchParams);
				site.setSequenceRegion(sequenceRegion);
				site.setDistanceToIntervalMean(Math.abs(sequenceRegion.getMean()-site.getPosition()));
				restrictionSites.add(site);
				sequenceRegion.addRestrictionSite(site);
				this.allRestrictionSites.add(site);
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		restrictionSites = RestrictionFragmentSizeFilter.excludeSitesMinMaxLength(restrictionSites, searchParams.getMIN_RESTRICTION_FRAGMENT_LENGTH(), searchParams.getMAX_RESTRICTION_FRAGMENT_LENGTH());
		
		RestrictionSite[] result = new RestrictionSite[restrictionSites.size()];
		restrictionSites.toArray(result);
		return result;
	}
	
	/**
	 * Computes the 'nbEnzymes' best (most-homogenuously distributed on 'sequence') restriction sites.
	 * 
	 * @param enzyme the restriction enzyme 'generating' the restriction sites
	 * @param nbSites the number of restriction sites to return
	 * @param searchParams the 3PD search parameters
	 * 
	 * @return the 'nbSites' most-homogenuously distributed restriction sites of restriction enzyme 'enzyme' on sequence 'sequence'
	 */
	public FastVector getMostHomogenuousRSSs(RestrictionEnzyme enzyme, int nbSites, PrimerSearchParameters searchParams){
		FastVector bestRSSs = new FastVector();
		if(searchParams.getNumContigs() == 1){
			assert(nbSites >= 1);
			String sequence = new String(searchParams.getContigs()[0].getSequence());
			if(sequence.length() < enzyme.getRecognitionSite().length()) throw new IllegalArgumentException("Sequence length must be >= enzyme length!");
			if(nbSites <= 0) throw new IllegalArgumentException("nbSites must be > 0!");
					
			int length = sequence.length();
			int intervalLength = MyExtendedMath.round(((double)length)/ nbSites);
			if(nbSites == 1){
				SequenceRegion sequenceRegion = new SequenceRegion(searchParams.getContigs()[0], 0, sequence.length());
				RestrictionSite site = findBestRestrictionSite(sequence, enzyme, sequenceRegion, 0, searchParams);
				if(site != null) bestRSSs.addElement(site);
				else throw new NoRestrictionSitesFoundException("No restriction sites for enzyme " + enzyme.getName() + " can be found in the input sequence!");
			}else{
			
				//for(int i=0; i<= length - intervalLength; i += intervalLength){
				//for(int i=0; i< length - intervalLength + 2; i += intervalLength){
				int start;
				int end = -1;
				for(int i=0; i<nbSites; i++){
					start = end + 1;
					end = Math.min(start + intervalLength, length);
	//				String interval = sequence.substring(i, Math.min(i+intervalLength, length));
	//				SequenceRegion sequenceRegion = new SequenceRegion(i,i+intervalLength-1);
	//				RestrictionSite bestSite = findBestRestrictionSite(interval, enzyme, sequenceRegion, i, searchParams);
					String interval = sequence.substring(start, end);
					SequenceRegion sequenceRegion = new SequenceRegion(searchParams.getContigs()[0], start, end);
					RestrictionSite bestSite = findBestRestrictionSite(interval, enzyme, sequenceRegion, start, searchParams);
					if(bestSite != null){
						bestSite.getSequenceRegion().sortRestrictionSites();
						bestRSSs.addElement(bestSite);
					}				
					else throw new NoRestrictionSitesFoundException("No restriction sites for enzyme " + enzyme.getName() + " can be found in the input sequence in interval " + sequenceRegion.getSeqRegionStart() + "-" + sequenceRegion.getSeqRegionEnd());
				}
			}
		}
		else{
			assert(nbSites >= 2);
			
			if(nbSites != searchParams.getNumContigs()) throw new IllegalArgumentException("nbSites must be == #Contigs for targeted search mode!");
			SimpleContig currentContig;
			String sequence;
			int minFragmentLength = searchParams.getMAX_AMPLICON_LENGTH()/2;
			
			for(int i=0; i<searchParams.getNumContigs(); i++){
				
				currentContig = searchParams.getContigs()[i];
				sequence = new String(currentContig.getSequence());
				
				if(currentContig.getSequenceLength() < searchParams.getMAX_AMPLICON_LENGTH()) throw new IllegalArgumentException("The sequence of contig: " + currentContig.getID() + " must have length > MAX_AMPLICON_LENGTH/2 on either side of the restriction site!");
				
				SequenceRegion sequenceRegion = new SequenceRegion(currentContig, 0, currentContig.getSequenceLength());
				RestrictionSite bestSite = findBestRestrictionSite(sequence, enzyme, sequenceRegion, 0, searchParams);
				if(bestSite.getPosition() - minFragmentLength < 0 || bestSite.getPosition() + minFragmentLength > sequence.length()) throw new IllegalArgumentException("The sequence of contig: " + currentContig.getID() + " must have length > MAX_AMPLICON_LENGTH/2 on either side of the restriction site!");
				
				if(bestSite != null){
					bestSite.getSequenceRegion().sortRestrictionSites();
					bestRSSs.addElement(bestSite);
				}				
				else throw new NoRestrictionSitesFoundException("No restriction sites for enzyme " + enzyme.getName() + " can be found in the input sequence of contig: " + currentContig.getID());
			}
		}
		if(bestRSSs.size() == 0) throw new EmptyResultSetException("No restriction sites for enzyme " + enzyme.getName() + " could be found in the input sequence.");
		
		return bestRSSs;
	}
	
	/**
	 * Main method for testing the RSS scan
	 * 
	 * @param args args[0] = input sequence to scan for RSSs
	 */
	public static void main(String[] args){
		
		RestrictionSiteSearch sites = new RestrictionSiteSearch();
		
		String[] sequences = {"100kb.dna", "200kb.dna", "400kb.dna", "800kb.dna"};
		String path = "/export/Sebastian/PrimerDesign/Testsequenzen/";
		final PrimerSearchParameters searchParams = new PrimerSearchParameters();
		
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		for(int i=0; i<sequences.length; i++){
			try{
				System.out.println("Reading sequence " + sequences[i]);
				SlimFastaParser parser = new SlimFastaParser(path + sequences[i]);
				SimpleContig contig = parser.parseNextContigIgnoreCase();
				System.out.println("Reading sequence " + sequences[i] + " took " + timer.getTimeString());
				
				searchParams.setContigs(new SimpleContig[]{contig});
				
				System.out.println("Scanning sequence " + sequences[i] + " for 4 homogenous restriction sites");
				RestrictionEnzyme enzyme = new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0);
				FastVector hits = sites.getMostHomogenuousRSSs(enzyme, 4, searchParams);
				System.out.println("Scanning sequence " + sequences[i] + " for 4 homogenous restriction sites took " + timer.getTimeString());
				System.out.println("The input sequence " + sequences[i] + " of length " + contig.getSequenceLength() + " contains " + hits.size() + " instances of restriction sites " + enzyme.getName());
				
				for(int j=0; j< hits.size(); j++){
					System.out.println("\tat location " + ((RestrictionSite) hits.elementAt(j)).getPosition());
				}								
			}
			catch(Exception e){
				e.printStackTrace();
			}	
			System.out.println();
		}
	}
}
