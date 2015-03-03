package primerDesign.dsc.indexStructures;

import java.io.File;

import primerDesign.dsc.indexStructures.esa.EnhancedSuffixArrayFatOpt;
import primerDesign.dsc.indexStructures.esa.MultiSeqMemoryMappedESAIndex;
import primerDesign.dsc.indexStructures.rmi.MultiSeqESAFatOpt;

/**
 * Deserializer for several kinds of DNA sequence indices.
 * 
 * @author Sebastian Fršhler
 *
 */
public class DNASequenceIndexDeserializer {
	/**
	 * Deserializes a DNA sequence index from a file.
	 * 
	 * @param filname the file containing the DNA sequence index
	 * 
	 * @return the DNA sequence index stored in a file
	 */
	public static DNASequenceIndex deserialize(String filname){
		DNASequenceIndex index;
		if(filname.endsWith("MultiSeqESAFatOpt.esaidx")){
			// multi sequence ESAFatOpt
			index = MultiSeqESAFatOpt.deserialize(new File(filname));
		}
		else if(filname.endsWith("ESAFatOpt.esaidx")){
			// single sequence ESAFatOpt
			index = EnhancedSuffixArrayFatOpt.deserialize(filname);
		}
		else if(filname.endsWith("MultiSeqMemoryMappedESAIndex.esaidx")){
			// multi sequence MemoryMappedESAFatOpt
			index = MultiSeqMemoryMappedESAIndex.deserialize(new File(filname));
		}
		else if(filname.endsWith("DummyIndex.esaidx")){
			index = DummyIndex.deserialize(new File(filname));
		}
		else{
			throw new IllegalArgumentException("Unsupported index type!: " + filname);
		}
		return index;
	}
}
