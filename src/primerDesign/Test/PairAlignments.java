/**
 * 
 */
package primerDesign.Test;

import primerDesign.dsc.SequenceRegionAlignment;

/**
 * Encapsulates a set of pair alignments.
 * 
 * @author froehler
 *
 */
public class PairAlignments {
	private SequenceRegionAlignment fwRev;
	private SequenceRegionAlignment fwProbe;
	private SequenceRegionAlignment revProbe;
	
//	public PairAlignments(String forwardRegion, String reverseRegion, String probeRegion){
//		this.fwRev = SequenceRegionAligner.alignSequenceRegions(forwardRegion.toCharArray(), reverseRegion.toCharArray());
//		this.fwProbe = SequenceRegionAligner.alignSequenceRegions(forwardRegion.toCharArray(), probeRegion.toCharArray());
//		this.revProbe = SequenceRegionAligner.alignSequenceRegions(reverseRegion.toCharArray(), probeRegion.toCharArray());
//	}
//	
//	public int getForwardReverseAlignment(int posForward, int lengthForward, int posReverse, int lengthReverse){
//		return this.fwRev.getGlobalAlignmentValue(posForward, lengthForward, posReverse, lengthReverse, false);
//	}
//	
//	public int getForwardProbeAlignment(int posForward, int lengthForward, int posProbe, int lengthProbe){
//		return this.fwProbe.getGlobalAlignmentValue(posForward, lengthForward, posProbe, lengthProbe, false);
//	}
//	
//	public int getReverseProbeAlignment(int posReverse, int lengthReverse, int posProbe, int lengthProbe){
//		return this.revProbe.getGlobalAlignmentValue(posReverse, lengthReverse, posProbe, lengthProbe, false);
//	}
//	
//	public int getAlignment(Primer first, Primer second){
//		Enum<PrimerTypes> typeFirst = first.getPrimerType();
//		Enum<PrimerTypes> typeSecond = second.getPrimerType();
//		
//		if((typeFirst.equals(PrimerTypes.forwardPrimer) && typeSecond.equals(PrimerTypes.reversePrimer)) || ((typeSecond.equals(PrimerTypes.forwardPrimer) && typeFirst.equals(PrimerTypes.reversePrimer)))){
//			return getForwardReverseAlignment(posForward, first.getLength(), posReverse, second.getLength());
//		}
//		else if((typeFirst.equals(PrimerTypes.forwardPrimer) && typeSecond.equals(PrimerTypes.hybridizationProbe)) || (typeSecond.equals(PrimerTypes.forwardPrimer) && typeFirst.equals(PrimerTypes.hybridizationProbe))){
//			return getForwardProbeAlignment(posForward, first.getLength(), posProbe, second.getLength());
//		}
//		else if((typeFirst.equals(PrimerTypes.reversePrimer) && typeSecond.equals(PrimerTypes.hybridizationProbe)) || (typeSecond.equals(PrimerTypes.reversePrimer) && typeFirst.equals(PrimerTypes.hybridizationProbe)){
//			return getReverseProbeAlignment(posReverse, first.getLength(), posProbe, second.getLength());
//		}
//	}
}
