/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.SimpleContigImpl;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a wrapper around an index hit.
 * 
 * @author Sebastian Fršhler
 *
 */
public class DummyRemoteIndexStructure implements RemoteIndexStructureSearch {
	
	/**
	 * Connects to a remote index structure, creates a reference to this remote index structure.
	 * 
	 * @param adress the ip-adress of the computer hosting the remote index structure
	 * @param port the port to use
	 * @param serviceName the name of the remote index structure
	 */
	public DummyRemoteIndexStructure(String adress, int port, String serviceName){
		
	}

	/* (non-Javadoc)
	 * @see esaRmi.RemoteIndexStructureSearch#findHitPositions(java.lang.String)
	 */
	public ObjectArrayList findHitPositions(String sequence){
		ObjectArrayList result = new ObjectArrayList();
		result.add(new IndexHitImpl(new SimpleContigImpl("Contig1","abcd".toCharArray()), 10, true));
		return result;
	}
}
