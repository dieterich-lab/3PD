/**
 * 
 */
package primerDesign.dsc.indexStructures;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Set;

import primerDesign.util.SimpleContigImpl;
import cern.colt.list.ObjectArrayList;

/**
 * Defines properties and methods common to all multi-sequence index structures.
 * 
 * @author Sebastian Fršhler
 *
 */
public abstract class MultiSeqIndex implements DNASequenceIndex, Serializable {

	private static final long serialVersionUID = 1L;
	private HashMap<String, DNASequenceIndex> indices;
	private String name;
	protected File fastaFile;
	
	/**
	 * Initializes a multi sequence index structure.
	 * 
	 * @param file
	 */
	public MultiSeqIndex(File file){
		this.indices = new HashMap<String, DNASequenceIndex>();
		this.fastaFile = file;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#createIndex(java.io.File)
	 */
	public abstract void createIndex();
	
	/**
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getName()
	 */
	public String getName(){
		return this.name;
	}
	
	/**
	 * Returns the number of indices stored in this multi sequence index.
	 * 
	 * @return the number of indices stored in this multi sequence index
	 */
	public int getNumIndices(){
		return this.indices.size();
	}
	
	/**
	 * Returns a specific index.
	 * 
	 * @param name the name of the index to return
	 * @return the index with name 'name'
	 */
	public DNASequenceIndex getIndex(String name){
		return this.indices.get(name);
	}
	
	/**
	 * Adds a new sequence index.
	 * 
	 * @param name the name of the index to add
	 * @param index the index to add
	 */
	public void putIndex(String name, DNASequenceIndex index){
		this.indices.put(name, index);
	}
	
	/**
	 * Returns all names of the indices of this multi sequence index.
	 * 
	 * @return all names of the indices of this multi sequence index
	 */
	public Set<String> getIndicesKeyset(){
		return this.indices.keySet();
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#findHitCount(java.lang.String)
	 */
	public int findHitCount(String sequence) {
		int count = 0;
		// count matches in all contigs
		for(String contig : this.indices.keySet()){
			count += this.indices.get(contig).findHitCount(sequence.toUpperCase());
		}
		return count;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#findHitPositions(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence) {
		ObjectArrayList result = new ObjectArrayList();
		ObjectArrayList currentList;
		// search matches in all contigs
		for(String contig : this.indices.keySet()){
			// retrieve match positions in current contig
			currentList = this.indices.get(contig).findHitPositions(sequence.toUpperCase());
			for(int i=0; i<currentList.size(); i++){
				result.add(currentList.get(i));
			}
		}
		return result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getContig()
	 */
	public SimpleContigImpl[] getContig() {
		SimpleContigImpl[] result = new SimpleContigImpl[this.indices.size()];
		int i=0;
		for(String name : this.indices.keySet()){
			result[i++] = getIndex(name).getContig()[0];  // a multi index consists of single indices, each constructed over exactly one! contig
		}
		return result;
	}
	
	/**
	 * Returns the sequence length.
	 * 
	 * @return the sequence length
	 */
	public int getSequenceLength() {
		int result = 0;
		for(String contig : this.indices.keySet()){
			result += this.indices.get(contig).getContig()[0].getSequenceLength();
		}
		return result;
	}
	
	/**
	 * Returns some statistics of single character frequency in the sequence
	 * 
	 * @return a hasmap of: character to integer frequency counts
	 */
	public HashMap<Character, Integer> getStatistics(){
		HashMap<Character, Integer> result = new HashMap<Character, Integer>();
		
		for(String contig : this.indices.keySet()){
			HashMap<Character, Integer> temp = this.indices.get(contig).getStatistics();
			for(Character x : temp.keySet()){
				if(result.containsKey(x)){
					result.put(x, result.get(x) + temp.get(x));
				}
				else{
					result.put(x, temp.get(x));
				}
			}
		}
		
		return result;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#serialize(java.io.File)
	 */
	public void serialize(File file) {
		try{
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(file));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
}
