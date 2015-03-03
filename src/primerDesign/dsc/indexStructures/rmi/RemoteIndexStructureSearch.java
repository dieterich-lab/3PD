/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.rmi.Remote;
import java.rmi.RemoteException;

import cern.colt.list.ObjectArrayList;

/**
 * @author Sebastian Fršhler
 *
 */
public interface RemoteIndexStructureSearch extends Remote{
	
	/**
	 * Searches (via rmi) a remote index structure for hits of sequence 'sequence'.
	 * 
	 * @param sequence the sequence to scan the index structure for
	 * 
	 * @return a list of hits of sequence 'sequence' in the queried index structure
	 * 
	 * @throws RemoteException
	 */
	public ObjectArrayList findHitPositions(String sequence) throws RemoteException;
}
