/**
 * 
 */
package primerDesign.dsc.indexStructures;

import java.io.File;
import java.util.HashMap;

import primerDesign.util.SimpleContigImpl;
import cern.colt.list.ObjectArrayList;

/**
 * Specifies an index structure over at least one contig.
 * 
 * An index structure covering more than one contig has to handle each contig separately
 * in order to provide the name of the contig a hit is found on.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface DNASequenceIndex {
	/**
	 * Creates an index structure covering all sequences in file 'file'.
	 * 
	 * An index structure has to provide methods to create and query the index structure.
	 * As a convencience, the index also has to be serializable in order to re-use it.
	 * 
	 */
	public abstract void createIndex();
	
	/**
	 * Returns the name of the contig, this index was created on.
	 * 
	 * The name of a single sequence index is the fasta id,
	 * the name of a multi sequence index is the name of the file containing all sequences.
	 * 
	 * @return the name of the contig, this index was created on
	 */
	public abstract String getName();
	
	/**
	 * Returns the contig this index was created on.
	 * 
	 * @return the contig this index was created on
	 */
	public abstract SimpleContigImpl[] getContig();
	
	/**
	 * Searches for hits of query 'sequence' in an index structure.
	 * 
	 * Each hit consists of: the name of the contig the hit was found on AND the position of the hit on that contig.
	 * 
	 * @param sequence the sequence to query the index structure with
	 * 
	 * @return a list of all hits found in the index structure, each hit consists of: the name of the contig the hit was found on AND the position of the hit on that contig
	 */
	public abstract ObjectArrayList findHitPositions(String sequence);
	
	/**
	 * Searches for the number of hits of query 'sequence' in an index structure.
	 * 
	 * @param sequence the sequence to query the index structure with
	 * 
	 * @return the number of hits of query 'sequence' in the index structure
	 */
	public abstract int findHitCount(String sequence);
	
	//public abstract String getSequence();
	
	/**
	 * Serializes the index structure to file 'file'.
	 * 
	 * @param file the file to serialize the index structure to
	 */
	public abstract void serialize(File file);
	
	/**
	 * Returns some statistics about a DNA sequence index.
	 * 
	 * @return some statistics about a DNA sequence index
	 */
	public abstract HashMap<Character, Integer> getStatistics();
}
