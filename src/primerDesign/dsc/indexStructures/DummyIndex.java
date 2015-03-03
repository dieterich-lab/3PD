/**
 * 
 */
package primerDesign.dsc.indexStructures;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;

import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * A dummy index always returning NO hits at all - just for debugging purpose.
 * 
 * @author Sebastian Fršhler
 *
 */
public class DummyIndex implements DNASequenceIndex, primerDesign.dsc.DNASequenceIndex, Serializable {

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#createIndex()
	 */
	public void createIndex() {
		// TODO Auto-generated method stub

	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#findHitCount(java.lang.String)
	 */
	public int findHitCount(String sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#findHitPositions(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence) {
		// TODO Auto-generated method stub
		return new ObjectArrayList();
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getContig()
	 */
	public SimpleContigImpl[] getContig() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getName()
	 */
	public String getName() {
		// TODO Auto-generated method stub
		return "Dummy Index!";
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getStatistics()
	 */
	public HashMap<Character, Integer> getStatistics() {
		// TODO Auto-generated method stub
		return null;
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
	
	public static DummyIndex deserialize(File filename){
		DummyIndex result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (DummyIndex) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#createIndex(java.lang.String, int, boolean)
	 */
	public void createIndex(String sequence, int maxWordSize,
			boolean includesScanRegion) {
		// TODO Auto-generated method stub
		
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#getMaxWordSize()
	 */
	public int getMaxWordSize() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#getSequence()
	 */
	public String getSequence() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#includesScanRegion()
	 */
	public boolean includesScanRegion() {
		// TODO Auto-generated method stub
		return true;
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchMatchPositionsInIndex(java.lang.String)
	 */
	public Integer[] searchMatchPositionsInIndex(String searchString) {
		// TODO Auto-generated method stub
		return new Integer[0];
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.DNASequenceIndex#searchNbMatchesInIndex(java.lang.String)
	 */
	public int searchNbMatchesInIndex(String searchString) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String mode = args[0];
		String file = args[1];
		if(mode.equals("create")){
			DummyIndex index = new DummyIndex();
			
			System.out.println("Creating single indices from file " + file);
			index.createIndex();
			System.out.println("Creating single indices from file " + file + " - done in " + timer.getTimeString());

			System.out.print("Serializing index");
			index.serialize(new File(file + ".DummyIndex.esaidx"));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.println("Index construction done in " + timer.getTotalTimestring());
		}
		else if(mode.equals("query")){
			String query = args[2];
			int numQueries = Integer.parseInt(args[3]);
			System.out.print("Deserializing index");
			DNASequenceIndex index = DummyIndex.deserialize(new File(file));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for string " + query);
			int result = 0;
			for(int i=0; i<numQueries; i++){
					result = index.findHitCount(query);
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
		}
		else if(mode.equals("randomquery")){
			int numQueries = Integer.parseInt(args[2]);
			System.out.print("Deserializing index");
			DummyIndex index = DummyIndex.deserialize(new File(file));
			System.out.println(" - done in " + timer.getTimeString());
			System.out.print("Querying index " + numQueries + " times for random string of length " + 18 + "-" + 30);
			int result = 0;
			for(int i=0; i<numQueries; i++){
					result = index.findHitCount(SeqTools.getRandomPrimerSequence(18, 30));
			}
			System.out.println(" " + result + " matches - done in " + timer.getTimeString());
			System.out.print("Creating " + numQueries + " random primers itself took");
			for(int i=0; i<numQueries; i++){
				SeqTools.getRandomPrimerSequence(18, 30);
			}
			System.out.println(" " + timer.getTimeString());
		}
		else{
			System.out.println("Usage: DummyIndex query <indexFile> <queryString> <numQueries>");
			System.out.println("or");
			System.out.println("Usage: DummyIndex create <serializedIndexFile>");
			System.out.println("\tNote: Filenames have to be provided including absolute path to file!");
		}
	}

}
