/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.PrimerPair;
import primerDesign.dsc.PrimerPairSet;
import cern.colt.list.ObjectArrayList;

/**
 * A class for iterating over a list of primer pairs in a sorted manner!
 * 
 * Primer pairs are returned by increasing distance to the optimal distance to the primer pair set provided.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPairIterator {
	private ObjectArrayList pairs;
	private PrimerPairSet pairSet;
	private int low;
	private int high;
	private int current;
	private boolean isInit;
	
	/**
	 * Initializes a primer pair iterator.
	 * 
	 * @param pairs the primer pairs to iterate over
	 * @param pairSet the primer pair to compare and order each element in 'pairs' with
	 */
	public PrimerPairIterator(ObjectArrayList pairs, PrimerPairSet pairSet){
		this.pairs = pairs;
		this.pairSet = pairSet;
		this.isInit = false;
	}
	
	/**
	 * Returns the next primer pair w.r.t distance.
	 * 
	 * @return the next primer pair w.r.t distance
	 */
	public PrimerPair getNext(){
		if(!this.isInit){
			// find pair with d(pairs(x), pair) = argmin
			this.low = 0;
			this.high = pairs.size() -1;
			this.current = (this.high + this.low)/2;
			while(this.low != this.high){
				if(((PrimerPair)pairs.getQuick(current)).getDistanceToOptimalPrimerPair() < pairSet.getAvgDistOptPrimerPair()){
					if(this.low == this.current) this.current++;
					this.low = this.current;
					this.current = Math.round((1.0F + this.high + this.low)/2);
				}
				else{
					if(this.high == this.current) this.current--;
					this.high = this.current;
					this.current = Math.round((1.0F + this.high + this.low)/2);
				}
			}
			this.low = this.current - 1;
			this.high = this.current +1;
			this.isInit = true;
			return ((PrimerPair)this.pairs.getQuick(this.current));
		}
		else if(this.low < 0 && this.high < this.pairs.size()) return ((PrimerPair) this.pairs.getQuick(this.high++));
		else if(this.low >= 0 && this.high >= this.pairs.size()) return ((PrimerPair) this.pairs.getQuick(this.low--));
		else if(this.low < 0 && this.high >= this.pairs.size()) throw new IllegalStateException("No more pair to output - complete list scanned!");
		else if(Math.abs(((PrimerPair) this.pairs.getQuick(this.low)).getDistanceToOptimalPrimerPair() - this.pairSet.getAvgDistOptPrimerPair()) < 		Math.abs(((PrimerPair) this.pairs.getQuick(this.high)).getDistanceToOptimalPrimerPair() - this.pairSet.getAvgDistOptPrimerPair())){
			return ((PrimerPair) this.pairs.getQuick(this.low--));
		}
		else{
			return ((PrimerPair) this.pairs.getQuick(this.high++));
		}
	}
	
	public boolean hasNext(){
		return this.low >= 0 || this.high < this.pairs.size();
	}
}
