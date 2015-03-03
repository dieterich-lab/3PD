/**
 * 
 */
package primerDesign.Test;

import org.biojava.bio.program.ssaha.SearchException;
import org.biojava.bio.program.ssaha.SearchListener;

import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.SimpleContigImpl;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class SSAHA_SearchListener implements SearchListener{
	
	private SSAHA_Index index;
	private ObjectArrayList hits;
	private boolean isForwardSearch;
	
	public SSAHA_SearchListener(ObjectArrayList hits, SSAHA_Index index, boolean isForwardSearch){
		this.hits = hits;
		this.index = index;
		this.isForwardSearch = isForwardSearch;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.program.ssaha.SearchListener#endSearch(java.lang.String)
	 */
	@Override
	public void endSearch(String seqID) {
		// TODO Auto-generated method stub
		
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.program.ssaha.SearchListener#hit(int, int, int, int)
	 */
	@Override
	public void hit(int hitID, int queryOffset, int hitOffset, int hitLength) {
		try {
			this.hits.add(new IndexHitImpl(new SimpleContigImpl(this.index.getDataStore().seqNameForID(hitID), null), hitOffset, this.isForwardSearch));
		} catch (IndexOutOfBoundsException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SearchException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.program.ssaha.SearchListener#startSearch(java.lang.String)
	 */
	@Override
	public void startSearch(String seqID) {
		
	}
}
