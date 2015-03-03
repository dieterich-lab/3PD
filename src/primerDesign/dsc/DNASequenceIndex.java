package primerDesign.dsc;

import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * Specifies a DNA sequence index.
 * 
 * @author Sebastian Fršhler
 *
 */
public interface DNASequenceIndex {

	/**
	 * Creates a new suffix tree of sequence 'sequence'.
	 * 
	 * @param sequence the sequence to construct the suffix tree from
	 * @param maxWordSize the maximum word size of a suffix in the tree
	 * 
	 * @throws NoSuchElementException 
	 * @throws BioException
	 */
	public abstract void createIndex(String sequence, int maxWordSize, boolean includesScanRegion);

	/**
	 * Queries the suffix tree for number of occurrences of 'searchString'.
	 * 
	 * @param searchString the query string whose number of occurrences should be retrieved
	 * @return the number of occurrences of search string 'searchString' - 0 iff searchString cannot be found!
	 * 
	 * @throws IllegalSymbolException
	 * @throws NoSuchElementException
	 * @throws BioException
	 */
	public abstract int searchNbMatchesInIndex(String searchString);
	
	/**
	 * Queries the suffix tree for match positions of 'searchString'
	 * 
	 * @param searchString the query string whose number of occurrences should be retrieved
	 * 
	 * @return the positions of occurrences of search string 'searchString' - empty FastVector iff searchString cannot be found!
	 */
	public abstract Integer[] searchMatchPositionsInIndex(String searchString);

	/**
	 * Returns whether the region to scan for valid primers is included in the region, this index is constructed from.
	 * 
	 * @return true iff sequence region includes primer scan region, false else
	 */
	public abstract boolean includesScanRegion();
	
	/**
	 * Returns the maximum word size contained in this index.
	 * 
	 * @return the maximum word size contained in this index
	 */
	public int getMaxWordSize();
	
	/**
	 * Returns the sequence this index was created on.
	 * 
	 * @return the sequence this index was created on
	 */
	public abstract String getSequence();
	
//	/**
//	 * Serializes an already computed index structure for later-on reuse.
//	 * 
//	 * @param filename the filename to write the (compressed) serialized index structure to
//	 */
//	public abstract void serialize(String filename);
//	
//	/**
//	 * Deserializes an already computed index structure for reuse.
//	 * 
//	 * @param filename the filename to read the (compressed) serialized index structure from
//	 */
//	public abstract static DNASequenceIndex deserialize(String filename);
	
}