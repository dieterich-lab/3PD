/**
 * A wrapper class for construction and querying an SSAHA index as provided by biojava.
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.text.NumberFormat;

import org.biojava.bio.BioException;
import org.biojava.bio.program.ssaha.DataStore;
import org.biojava.bio.program.ssaha.NIODataStoreFactory;
import org.biojava.bio.program.ssaha.SearchException;
import org.biojava.bio.program.ssaha.SearchListener;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.DNANoAmbPack;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

import primerDesign.dsc.indexStructures.IndexHit;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContig;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class SSAHA_Index implements Serializable {
	private static final String SUFFIX = ".ssaha";
	private static final String OBJECT_SERIAL_SUFFIX = ".ser";
	private static final int WORD_SIZE = 10;
	private static final int IGNORE_THRESHOLD = 10000;
	private DataStore index;
	private HashSequenceDB seqDB;
	
	public SSAHA_Index(File file) throws IOException, IllegalIDException, ChangeVetoException, IllegalSymbolException, BioException{
		// init new dummy sequence db
		this.seqDB = new HashSequenceDB();
		
		// fill dummy db with sequences from fasta file
		SlimFastaParser parser = new SlimFastaParser(file);
		SimpleContig contig;
		while(parser.hasNextContig()){
			contig = parser.parseNextContigIgnoreCase();
			seqDB.addSequence(DNATools.createDNASequence(new String(contig.getSequence()), contig.getID()));
		}
		
		// create ssaha data store
		NIODataStoreFactory factory = new NIODataStoreFactory();
		this.index = factory.buildDataStore(new File(file.getAbsolutePath() + SUFFIX), seqDB, new DNANoAmbPack(DNATools.n()), WORD_SIZE, IGNORE_THRESHOLD);
	}
	
	private ObjectArrayList findHits(String query, boolean isForwardSearch) throws IllegalAlphabetException, IllegalSymbolException, SearchException, org.biojava.bio.program.ssaha.SearchException, IllegalIDException{
		if(query.length() < WORD_SIZE) throw new IllegalArgumentException("The query has to be >= the minimum word length of this index (" + WORD_SIZE + " in this case)!");
		
		ObjectArrayList hits = new ObjectArrayList();
		SearchListener listener = new SSAHA_SearchListener(hits, this, isForwardSearch);		
		//SearchListener listener = new SearchListener.Echo(System.out);
		
		if(isForwardSearch){
			this.index.search("Dummy", DNATools.createDNA(query), listener);
		}
		else{
			query = SeqTools.revcompDNA(query.toCharArray());
			this.index.search("Dummy", DNATools.createDNA(query), listener);
		}
		boolean isRealHit;
		for(int i=0; i<hits.size(); i++){
			//verify hit
			isRealHit = verifyHit((IndexHit) hits.get(i), query);
			// delete hit if not veryfiable
			if(!isRealHit) hits.remove(i);
			if(i > 0) i--;
		}
		
		return hits;
	}
	
	/**
	 * Searches index for hits of query 'query'.
	 * 
	 * @param query the query sequence
	 * @return a list of hits (->IndexHitImpl) of query 'query'
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 * @throws IllegalIDException
	 * @throws SearchException
	 */
	public ObjectArrayList findHits(String query) throws IllegalAlphabetException, IllegalSymbolException, IllegalIDException, SearchException{
		ObjectArrayList result = findHits(query, true);
		ObjectArrayList reverse = findHits(query, false);
		if(reverse.size() > 0) result.addAllOfFromTo(reverse, 0, reverse.size() - 1);
		
		return result;
	}
	
	public final DataStore getDataStore(){
		return this.index;
	}
	
	public final HashSequenceDB getSequenceDB(){
		return this.seqDB;
	}
	
	/**
	 * Returns a subsequence of a sequence of a queries hit.
	 * 
	 * @param hit the hit for which to return the sequence
	 * @param start the start of the sequence (includive)
	 * @param end the end of the sequence (inclusive)
	 * 
	 * @return a subsequence of a sequence of a queries hit
	 * @throws IllegalIDException 
	 * @throws IndexOutOfBoundsException 
	 */
	public String getSequence(IndexHit hit, int start, int end) throws IndexOutOfBoundsException, IllegalIDException{
		return this.seqDB.getSequence(hit.getContigName()).subStr(start, end);
	}
	
	/**
	 * Verifies whether a hit of length WORD_SIZE really is a hit of the complete query.
	 * 
	 * This is a mechanism required because of the way queries are searched by SSAHA.
	 * 
	 * @param hit the hit to verify
	 * @param query the query that generated the hit
	 * @return whether a hit of length WORD_SIZE really is a hit of the complete query
	 * @throws IllegalIDException
	 */
	private boolean verifyHit(IndexHit hit, String query) throws IllegalIDException{
		query = query.toLowerCase();
		String seq;
		int start;
		int end;
		if(hit.isForwardHit()){
			end = hit.getPosition() + query.length() - 1;
			if(end >= this.seqDB.getSequence(hit.getContigName()).length()) return false;
			seq = this.seqDB.getSequence(hit.getContigName()).subStr(hit.getPosition(), end);
		}
		else{
			start = hit.getPosition() + WORD_SIZE - query.length() + 1;
			end = hit.getPosition() + WORD_SIZE;
			if(start < 1 || end >= this.seqDB.getSequence(hit.getContigName()).length()) return false;
			seq = this.seqDB.getSequence(hit.getContigName()).subStr(start, end);
		}
		for(int i=0; i<query.length(); i++){
			if(seq.charAt(i) != query.charAt(i)) return false;
		}
		return true;
	}
	
	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#serialize(java.io.File)
	 */
	public void serialize(File file) {
		try{
			this.index = null;
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	

	/* (non-Javadoc)
	 * @see esaRmi.DNASequenceIndex#deserialize(java.io.File)
	 */
	public static SSAHA_Index deserialize(File file) {
		SSAHA_Index result = null;
		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file));
			result = (SSAHA_Index) in.readObject();
			NIODataStoreFactory factory = new NIODataStoreFactory();
			file = new File(file.getAbsolutePath().replace(OBJECT_SERIAL_SUFFIX, ""));
			result.index = factory.getDataStore(file);
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IllegalIDException, ChangeVetoException, IllegalSymbolException, IOException, BioException, SearchException, org.biojava.bio.program.ssaha.SearchException{
		SSAHA_Index index;
		
		String mode = args[0];
		String file = args[1];
		SimpleTimer timer = new SimpleTimer();
		NumberFormat format = NumberFormat.getInstance();
		
		if(mode.equals("create")){
			System.out.print("Creating index from file: " + file);
			index = new SSAHA_Index(new File(file));
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.print("Serializing index");
			index.serialize(new File(file + SUFFIX + OBJECT_SERIAL_SUFFIX));
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.println("Total runtime: " + timer.getTotalTimestring());
		}else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			
			System.out.print("Deserializing index");
			index = SSAHA_Index.deserialize(new File(file));
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.print("Querying index " + format.format(numQueries) + " times for query " + query);
			ObjectArrayList hits = new ObjectArrayList();
			for(int i=0; i<numQueries; i++){
				hits = index.findHits(query);
			}
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Found " + format.format(hits.size()) + " matches to query " + query);
		}
		else if(mode.equals("randomquery")){
			int numQueries = Integer.parseInt(args[2]);
			
			System.out.print("Deserializing index");
			index = SSAHA_Index.deserialize(new File(file));
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.print("Querying index " + format.format(numQueries) + " times for random string of length " + 18 + "-" + 30);
			ObjectArrayList hits = new ObjectArrayList();
			for(int i=0; i<numQueries; i++){
					hits = index.findHits(SeqTools.getRandomPrimerSequence(18, 30));
			}
			System.out.println(" " + format.format(hits.size()) + " matches - done in " + timer.getTimeString());
			System.out.print("Creating " + format.format(numQueries) + " random primers itself took");
			for(int i=0; i<numQueries; i++){
				SeqTools.getRandomPrimerSequence(18, 30);
			}
			System.out.println(" " + timer.getTimeString());
		}else{
			System.out.println("Unknown option: " + args[0]);
			System.exit(1);
		}
	}
}
