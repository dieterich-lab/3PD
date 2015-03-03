/**
 * 
 */
package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.TreeSet;

import org.biojava.bio.BioException;
import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionEnzymeManager;
import org.biojava.bio.molbio.RestrictionMapper;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.SimpleThreadPool;

import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class DigestScreener {
	private static final int MAX_NUM_THREADS = 1;
	private static RestrictionMapper mapper = new RestrictionMapper(new SimpleThreadPool(MAX_NUM_THREADS, true));
	
	public static void main(String[] args) throws BioException, IOException{
		
		if(args.length != 3){
			System.out.println("Usage: java ... DigestScreener <ScanSequences> <EnzymeList> <BackgroundSequences>");
			System.exit(1);
		}
		
		// read in sequences to screen for restriction sites
		SequenceDB scanSequences = SeqIOTools.readFasta(new FileInputStream(args[0]), DNATools.getDNA());
		SequenceDB backgroundSequences = SeqIOTools.readFasta(new FileInputStream(args[2]), DNATools.getDNA());
		
		BufferedReader reader = new BufferedReader(new FileReader(args[1]));
		
		String line;
		RestrictionEnzyme currentEnzyme;
		String[] values;
//		FastVector absentEnzymes = new FastVector();
		ObjectArrayList absentEnzymes = new ObjectArrayList();
		HashMap<String, Integer> backgroundMatches = new HashMap<String, Integer>();
		// read in and process restriction enzymes sequentially, scan each sequence for presence of restriction sites, report all non-occurring enzymes
		while((line = reader.readLine()) != null){
			values = line.split("\t");
			if(values.length == 4){
				currentEnzyme = new RestrictionEnzyme(values[0], DNATools.createDNA(values[1]), Integer.parseInt(values[2]), Integer.parseInt(values[3]));
				RestrictionEnzymeManager.register(currentEnzyme, new TreeSet());
				DigestScreener.mapper.clearEnzymes();
				DigestScreener.mapper.addEnzyme(currentEnzyme);
				
				SequenceIterator iter = scanSequences.sequenceIterator();
				int count = 0;
				while(iter.hasNext()){
					if(DigestScreener.mapper.annotate(iter.nextSequence()).countFeatures() > 0)	count++;
				}
//				if(count == 0) absentEnzymes.addElement(currentEnzyme);
				if(count == 0){
					absentEnzymes.add(currentEnzyme);
					
					// scan occurrences of restriction site for current enzyme in background sequence
					SequenceIterator iter2 = backgroundSequences.sequenceIterator();
					while(iter2.hasNext()){
						backgroundMatches.put(currentEnzyme.getName(), DigestScreener.mapper.annotate(iter2.nextSequence()).countFeatures());
					}
				}
			}
		}
		// print all enzymes having no restroction site in all of the sequences
		//absentEnzymes.sort();
		System.out.println("Enzymes having no restriction site in any of the sequences in file " + args[0] + ":");
		System.out.println("Format: Enzyme name, enzyme site, position forward, position backward, hits in background sequence");
		for(int i=0; i<absentEnzymes.size(); i++){
//			System.out.println(absentEnzymes.elementAt(i).toString());
			System.out.println(absentEnzymes.get(i).toString() + "\t" + backgroundMatches.get(((RestrictionEnzyme)absentEnzymes.get(i)).getName()));
		}
	}
	
	class MyRestrictionEnzyme extends RestrictionEnzyme implements Comparable{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * @throws IllegalAlphabetException 
		 * 
		 */
		public MyRestrictionEnzyme(String name, SymbolList site, int fw, int rev) throws IllegalAlphabetException {
			super(name, site, fw, rev);
		}

		public int compareTo(Object other) {
			return this.name.compareTo(((MyRestrictionEnzyme)other).getName());
		}
		
	}
}
