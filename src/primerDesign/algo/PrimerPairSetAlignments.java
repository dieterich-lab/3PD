/**
 * 
 */
package primerDesign.algo;

import primerDesign.dsc.AlignmentType;
import primerDesign.dsc.PrimerAlignmentScores;
import primerDesign.dsc.SequenceRegionAlignment;
import primerDesign.util.PrimerSearchParameters;

/**
 * Encapsulates an alignment between a set of primer pairs.
 * 
 * For each primer pair combination, the forward primers and the forward primers and the probes are aligned.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerPairSetAlignments {
	Object[][] forward1Forward2;
	Object[][] forward1Probe2;
	Object[][] forward2Probe1;
	
	Object[] forwardSequences;
	int numFwSeq;
	Object[] probeSequences;
	int numPrSeq;
	
	/**
	 * Initializes a primer pairs alignment.
	 * 
	 * @param elements the number of elements
	 */
	public PrimerPairSetAlignments(int elements){
		this.forward1Forward2 = new Object[elements][elements];
		this.forward1Probe2 = new Object[elements][elements];
		this.forward2Probe1 = new Object[elements][elements];
		
		this.forwardSequences = new Object[elements];
		this.probeSequences = new Object[elements];
		
 		this.numFwSeq = 0;
		this.numPrSeq = 0;
	}
	
	/**
	 * Initializes a primer pairs alignment.
	 * 
	 * @param other the primer pair set to initialize this set with
	 */
	public PrimerPairSetAlignments(PrimerPairSetAlignments other){
		this.forward1Forward2 = other.forward1Forward2.clone();
		this.forward1Probe2 = other.forward1Probe2.clone();
		this.forward2Probe1 = other.forward2Probe1.clone();
		
		this.forwardSequences = other.forwardSequences.clone();
		this.probeSequences = other.probeSequences.clone();
		
		this.numFwSeq = other.numFwSeq;
		this.numPrSeq = other.numPrSeq;
	}
	
	/**
	 * Adds a scan regions to the alignment set.
	 * 
	 * @param scanRegionForward the forward scan region to add
	 * @param scanRegionProbe the probe scan region to add
	 * @param searchParams the 3PD search parameters
	 */
	public void addScanRegion(char[] scanRegionForward, char[] scanRegionProbe, PrimerSearchParameters searchParams){
		char[] current;
		assert(numFwSeq == numPrSeq);
		if(numFwSeq > 0){
			for(int i=0; i<numFwSeq; i++){
				current = (char[])forwardSequences[i];
				this.forward1Forward2[i][numFwSeq] = SequenceRegionAligner.alignSequenceRegions(current, scanRegionForward, searchParams.getA_t_basepair_score(), searchParams.getG_c_basepair_score());
				if(searchParams.isPickTaqManProbe()){
					this.forward1Probe2[i][numFwSeq] = SequenceRegionAligner.alignSequenceRegions(current, scanRegionProbe, searchParams.getA_t_basepair_score(), searchParams.getG_c_basepair_score());
					this.forward2Probe1[i][numFwSeq] = SequenceRegionAligner.alignSequenceRegions(scanRegionForward, (char[])probeSequences[i], searchParams.getA_t_basepair_score(), searchParams.getG_c_basepair_score());
				}
			}
		}
		this.forwardSequences[numFwSeq++] = scanRegionForward;
		this.probeSequences[numPrSeq++] = scanRegionProbe;
	}
	
	public void deleteLastRegion(){
		this.numFwSeq--;
		this.numPrSeq--;
	}
	
	/**
	 * Returns a specific sub-alignment.
	 * 
	 * @param regionIndex1 the index of region 1
	 * @param regionPos1 the start position within region 1
	 * @param regionLength1 the length of the alignment in region 1
	 * @param regionIndex2 the index of region 2
	 * @param regionPos2 the start position within region 2
	 * @param regionLength2 the length of the alignment in region 2
	 * @param type the alignment type
	 * 
	 * @return a specific sub-alignment from a larger alignment
	 */
	public PrimerAlignmentScores getAlignment(int regionIndex1, int regionPos1, int regionLength1, int regionIndex2, int regionPos2, int regionLength2, Enum<AlignmentType> type){
		if(regionIndex2 < regionIndex1) throw new IllegalArgumentException("Index 1 must be < index 2!");
		if(regionIndex1 >= numFwSeq || regionIndex2 >= numFwSeq || regionIndex1 >= numPrSeq || regionIndex2 >= numPrSeq) throw new ArrayIndexOutOfBoundsException();
		if(type.equals(AlignmentType.forward1Forward2)){
			return ((SequenceRegionAlignment)this.forward1Forward2[regionIndex1][regionIndex2]).getGlobalAlignmentValues(regionPos1, regionLength1, regionPos2, regionLength2);
		}
		else if(type.equals(AlignmentType.forward1Probe2)){
			return ((SequenceRegionAlignment)this.forward1Probe2[regionIndex1][regionIndex2]).getGlobalAlignmentValues(regionPos1, regionLength1, regionPos2, regionLength2);
		}
		else if(type.equals(AlignmentType.forward2Probe1)){
			return ((SequenceRegionAlignment)this.forward2Probe1[regionIndex1][regionIndex2]).getGlobalAlignmentValues(regionPos1, regionLength1, regionPos2, regionLength2);
		}
		else throw new IllegalArgumentException("Unsupported alignment type");
	}
}
