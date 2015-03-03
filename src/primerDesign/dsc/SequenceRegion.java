package primerDesign.dsc;

import java.util.Collections;
import java.util.Iterator;
import java.util.Vector;

import primerDesign.util.MyExtendedMath;
import primerDesign.util.SimpleContig;

/**
 * Encapsulates a stretch of dna and all of its restriction sites for a specific restriction enzyme.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SequenceRegion {
	private SimpleContig contig;
	private int start;
	private int end;
	private int mean;
	private Vector<RestrictionSite> restrictionSites = new Vector<RestrictionSite>();
	private boolean isSorted;
	private Iterator<RestrictionSite> restSiteIterator = null;
	private int iteratorPosition;
	
	/**
	 * Initializes a new sequence region.
	 * 
	 * @param contig the contig the region is on
	 * @param start the start position of the region
	 * @param end the end of the region
	 */
	public SequenceRegion(SimpleContig contig, int start, int end){
		if(end < start) throw new IllegalArgumentException("End must be >= start!");
		this.contig = contig;
		this.start = start;
		this.end = end;
		this.mean = MyExtendedMath.round(start + (end - start)/2);
		this.isSorted = false;
	}
	
	/**
	 * Adds a restriction site to this interval.
	 * 
	 * @param site a restriction site to be added to interval
	 */
	public void addRestrictionSite(RestrictionSite site){
		this.restrictionSites.add(site);
		//this.sortRestrictionSites();
		this.isSorted = false;
	}
	
	/**
	 * Removes a restriction site from the region.
	 * 
	 * @param site the site to remove
	 */
	public void removeRestrictionSite(RestrictionSite site){
		this.restrictionSites.remove(this.restrictionSites.indexOf(site));
	}
	
	/**
	 * Sorts the restriction sites in this interval in ascending order from interval mean.
	 *
	 */
	public void sortRestrictionSites(){
		if(!this.restrictionSites.isEmpty()){
			Collections.sort(restrictionSites);
			this.isSorted = true;
			this.restSiteIterator = null;
			this.iteratorPosition = 0;
		}
	}
	
	/**
	 * Returns an iterator over all restriction sites within this interval.
	 * 
	 * Restriction sites are ordered ascending w.r.t. distance to sequence region mean.
	 * 
	 * @param removeFirstElement remove first element to obtain iterator over alternative best restriction sites?
	 * 
	 * @return an iterator over all restriction sites within this interval
	 */
	public Iterator<RestrictionSite> getRestrictionSitesIterator(boolean removeFirstElement){
		if(!this.isSorted){
			sortRestrictionSites();
			this.restSiteIterator = null;
			this.iteratorPosition = 0;
		}
		if(this.restSiteIterator == null){
			this.restSiteIterator = restrictionSites.iterator();
			this.iteratorPosition = 0;
		}
		if(removeFirstElement && this.restSiteIterator.hasNext()){
			this.restSiteIterator.next();
			this.iteratorPosition++;
		}
		return this.restSiteIterator;
	}
	
	/**
	 * Returns the current restriction site iterator being at the position of the last restriction site queried.
	 * 
	 * @return the current restriction site iterator being at the position of the last restriction site queried
	 */
	public Iterator<RestrictionSite> getCurrentRestrictionSitesIterator(){
		return this.restSiteIterator;
	}
	
//	public int getIteratorPosition(){
//		return this.iteratorPosition;
//	}
	
	public int getRestricionSiteIndex(RestrictionSite rss){
		return this.restrictionSites.indexOf(rss);
	}
	
	public int getNumRestrictionSites(){
		return this.restrictionSites.size();
	}
	
	/**
	 * Returns the mean of this sequence region.
	 * 
	 * @return the mean of this sequence region
	 */
	public int getMean(){
		return this.mean;
	}
	
	/**
	 * Returns the start of this sequence region.
	 * 
	 * @return the start of this sequence region
	 */
	public int getSeqRegionStart(){
		return start;
	}
	
	/**
	 * Returns the end of this sequence region.
	 * 
	 * @return the end of this sequence region
	 */
	public int getSeqRegionEnd(){
		return end;
	}
	
	public SimpleContig getContig(){
		return this.contig;
	}
	
	public String getContigSequence(){
		return new String(this.contig.getSequence());
	}
	
	public int getContigSequenceLength(){
		return this.contig.getSequenceLength();
	}
	
	public boolean equals(SequenceRegion otherRegion){
		if(this.start == otherRegion.start && this.end == otherRegion.end
				&& this.mean == otherRegion.mean && this.restrictionSites.equals(otherRegion.restrictionSites)) return true;
		else return false;
	}
}
