/**
 * 
 */
package primerDesign.dsc.indexStructures.primerMisprimingCheck;

import java.io.File;

import primerDesign.dsc.indexStructures.DummyIndex;
import primerDesign.dsc.indexStructures.blat.BLATQueryClient;
import primerDesign.dsc.indexStructures.blat.RestrictionSitesIndex;
import primerDesign.dsc.indexStructures.esa.EnhancedSuffixArrayFatOpt;
import primerDesign.dsc.indexStructures.esa.MultiSeqMemoryMappedESAIndex;
import primerDesign.dsc.indexStructures.rmi.MultiSeqESAFatOpt;
import primerDesign.util.PrimerSearchParameters;

/**
 * Deserializes a primer mispriming check to be used during primer mispriming checks.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerMisprimingCheckDeserializer {
	/**
	 * Deserializes an index structure to be used during primer mispriming checks.
	 * 
	 * @param filename the filename containing the index structure
	 * @param params the primer search parameters to be used by 3PD
	 * 
	 * @return a primer mispriming check to be used during 3PD
	 */
	public static PrimerMisprimingCheck deserialize(String filename, PrimerSearchParameters params){
		PrimerMisprimingCheck checker;
		if(filename.endsWith("MultiSeqESAFatOpt.esaidx")){
			// multi sequence ESAFatOpt
			checker = new ESA3CPrimerMisprimingScan(MultiSeqESAFatOpt.deserialize(new File(filename)), params);
		}
		else if(filename.endsWith("ESAFatOpt.esaidx")){
			// single sequence ESAFatOpt
			checker = new ESA3CPrimerMisprimingScan(EnhancedSuffixArrayFatOpt.deserialize(filename), params);
		}
		else if(filename.endsWith("MultiSeqMemoryMappedESAIndex.esaidx")){
			// multi sequence MemoryMappedESAFatOpt
			checker = new ESA3CPrimerMisprimingScan(MultiSeqMemoryMappedESAIndex.deserialize(new File(filename)), params);
		}
		else if(filename.endsWith("DummyIndex.esaidx")){
			checker = new ESA3CPrimerMisprimingScan(DummyIndex.deserialize(new File(filename)), params);
		}
		else if(filename.endsWith("RSSIndex.ser")){
			//checker = new Blat3CPrimerMisprimingScan(new BLATQueryClient(null, null, null, null), RestrictionSitesIndex.deserialize(filename), params);
			checker = new Blat3CPrimerMisprimingScan(new BLATQueryClient(params.getValueFromConfigFile("BLAT_gfClient"), params.getValueFromConfigFile("BLAT_host"), params.getValueFromConfigFile("BLAT_PORT_" + params.getTargetOrganism()), params.getValueFromConfigFile("BLAT_SequenceDir")), RestrictionSitesIndex.deserialize(filename), params);
		}
		else if(filename.endsWith("SlimGenome.ser")){
			checker = new Blat3CPrimerMisprimingScanPlusSeq(new BLATQueryClient(params.getValueFromConfigFile("BLAT_gfClient"), params.getValueFromConfigFile("BLAT_host"), params.getValueFromConfigFile("BLAT_PORT_" + params.getTargetOrganism()), params.getValueFromConfigFile("BLAT_SequenceDir")), filename, params);
		}
		else{
			throw new IllegalArgumentException("Unsupported index type!: " + filename);
		}
		return checker;
	}
}
