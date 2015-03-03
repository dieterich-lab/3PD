/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.io.File;
import java.rmi.RemoteException;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.dsc.indexStructures.esa.MultiSeqMemoryMappedESAIndex;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a search in a remote index structure.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RemoteIndexStructureSearchESAmemoryMappedImpl implements RemoteIndexStructureSearch {

	public final DNASequenceIndex index;
	
	public RemoteIndexStructureSearchESAmemoryMappedImpl(File file) throws RemoteException{
		this.index = MultiSeqMemoryMappedESAIndex.deserialize(file);
	}
	/* (non-Javadoc)
	 * @see esaRmi.RemoteIndexStructureSearch#findHitPositions(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence) throws RemoteException {
		return this.index.findHitPositions(sequence);
	}
}
