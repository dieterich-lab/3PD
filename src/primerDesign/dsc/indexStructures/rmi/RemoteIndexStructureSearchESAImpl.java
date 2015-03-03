/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.io.File;
import java.rmi.RemoteException;
import java.rmi.server.UnicastRemoteObject;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import cern.colt.list.ObjectArrayList;

/**
 * Encapsulates an enhanced suffix array as a remote index structure.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RemoteIndexStructureSearchESAImpl extends UnicastRemoteObject implements RemoteIndexStructureSearch {

	private static final long serialVersionUID = 2890559659516611074L;

	private final DNASequenceIndex index;

	/**
	 * Reads in (de-serializes) a pre-existing index structure.
	 */
	public RemoteIndexStructureSearchESAImpl(File file) throws RemoteException{
		super();
		// read in serialized index structure
		this.index = MultiSeqESAFatOpt.deserialize(file);
	}
	
	/* (non-Javadoc)
	 * @see esaRmi.RemoteIndexStructureSearch#findHits(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence) throws RemoteException {
		return this.index.findHitPositions(sequence);
	}

}
