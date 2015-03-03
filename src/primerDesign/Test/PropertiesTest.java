/**
 * 
 */
package primerDesign.Test;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;

/**
 * This class reads in the parameters of the program as provided by an XML file.
 * 
 * @author froehler
 *
 */
public class PropertiesTest {
	public static void main(String[] args) throws FileNotFoundException, IOException{
		Properties configFile = PropertiesTest.getDefaultProperties();
		
		System.out.println("Writing properties");
		configFile.storeToXML(new FileOutputStream("/export/Sebastian/PrimerDesign/config.xml"), "SNP mining parameters");
		//configFile.save(new FileOutsetPropertyStream("/export/Sebastian/config.xml"), "My Comment");
		
//		System.out.println("Reading properties");
//		Properties inFile = new Properties();
//		inFile.loadFromXML(new FileInsetPropertyStream(new File("/export/Sebastian/config.xml")));
//		
//		System.out.println("Read properties:");
//		System.out.println(inFile.get("NumInstances"));
//		System.out.println(inFile.get("parameter"));
	}
	
	private static Properties getDefaultProperties(){
		Properties configFile = new Properties();
		
		configFile.setProperty("PRIMER_CONCENTRATION", Double.toString(500E-09));
		configFile.setProperty("TAQMAN_PROBE_CONCENTRATION", Double.toString(150E-09));
		configFile.setProperty("MONOVALENT_CATION_CONCENTRATION", Double.toString(0));
		configFile.setProperty("DIVALENT_CATION_CONCENTRATION", Double.toString(4E-03));
		configFile.setProperty("DNTP_CONCENTRATION", Double.toString(4E-03));
		
		configFile.setProperty("a_t_basepair_score", Integer.toString(2));
		configFile.setProperty("g_c_basepair_score", Integer.toString(4));
		
		configFile.setProperty("REPETITIVE_ELEMENT_CHARACTER", "N");
		
		configFile.setProperty("MAX_PRIMER_LENTH", Integer.toString(30));
		configFile.setProperty("MIN_PRIMER_LENTH", Integer.toString(18));
		configFile.setProperty("OPT_PRIMER_LENTH", Integer.toString(22));
		
		configFile.setProperty("TAQMAN_MAX_PRIMER_LENTH", Integer.toString(35));
		configFile.setProperty("TAQMAN_MIN_PRIMER_LENTH", Integer.toString(18));
		configFile.setProperty("TAQMAN_OPT_PRIMER_LENTH", Integer.toString(25));
		
		configFile.setProperty("MAX_AMPLICON_LENGTH", Integer.toString(200));
		configFile.setProperty("MIN_AMPLICON_LENGTH", Integer.toString(140));
		configFile.setProperty("OPT_AMPLICON_LENGTH", Integer.toString(150));
		
		configFile.setProperty("SAFE_FALSE_POSITIVE_AMPLICON_LENGTH", Integer.toString(500));
		configFile.setProperty("SAFE_MISPRIMING_SEQUENCE_PRECENT_CUTOFF", Double.toString(0.5));
		configFile.setProperty("MAX_PRIMER_MISPRIMING_CUTOFF", Integer.toString(0));
		
		configFile.setProperty("OPT_DISTANCE_TO_RSS", Integer.toString(75));
		
		configFile.setProperty("MAX_TM", Integer.toString(63));
		configFile.setProperty("MIN_TM", Integer.toString(58));
		configFile.setProperty("OPT_TM", Integer.toString(60));
		
		configFile.setProperty("TAQMAN_MAX_TM", Integer.toString(72));
		configFile.setProperty("TAQMAN_MIN_TM", Integer.toString(67));
		configFile.setProperty("TAQMAN_OPT_TM", Integer.toString(69));
		
		configFile.setProperty("MAX_GC", Double.toString(0.8));
		configFile.setProperty("MIN_GC", Double.toString(0.3));
		configFile.setProperty("OPT_GC", Double.toString(0.5));
		
		configFile.setProperty("TAQMAN_MAX_GC", Double.toString(0.8));
		configFile.setProperty("TAQMAN_MIN_GC", Double.toString(0.3));
		configFile.setProperty("TAQMAN_OPT_GC", Double.toString(0.5));
		
		configFile.setProperty("MAX_PRIMER_SELF_ALIGNMENT_SCORE", Integer.toString(24));
		configFile.setProperty("OPT_PRIMER_SELF_ALIGNMENT_SCORE", Integer.toString(0));
		configFile.setProperty("MAX_PRIMER_SELF_END_ALIGNMENT_SCORE", Integer.toString(8));
		configFile.setProperty("OPT_PRIMER_SELF_END_ALIGNMENT_SCORE", Integer.toString(0));
		
		configFile.setProperty("MAX_PRIMER_PAIR_ALIGNMENT_SCORE", Integer.toString(24));
		configFile.setProperty("OPT_PRIMER_PAIR_ALIGNMENT_SCORE", Integer.toString(0));
		configFile.setProperty("MAX_PRIMER_PAIR_END_ALIGNMENT_SCORE", Integer.toString(8));
		configFile.setProperty("OPT_PRIMER_PAIR_END_ALIGNMENT_SCORE", Integer.toString(0));
		
		configFile.setProperty("MAX_TAQMAN_SELF_ALIGNMENT_SCORE", Integer.toString(36));
		configFile.setProperty("OPT_TAQMAN_SELF_ALIGNMENT_SCORE", Integer.toString(0));
		configFile.setProperty("MAX_TAQMAN_SELF_END_ALIGNMENT_SCORE", Integer.toString(36));
		configFile.setProperty("OPT_TAQMAN_SELF_END_ALIGNMENT_SCORE", Integer.toString(0));
		
		configFile.setProperty("MAX_TAQMAN_PAIR_ALIGNMENT_SCORE", Integer.toString(36));
		configFile.setProperty("OPT_TAQMAN_PAIR_ALIGNMENT_SCORE", Integer.toString(0));
		configFile.setProperty("MAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE", Integer.toString(36));
		configFile.setProperty("OPT_TAQMAN_PAIR_END_ALIGNMENT_SCORE", Integer.toString(0));
		
		configFile.setProperty("PRIMER_DELTA_TM_WEIGHT", Double.toString(1));
		configFile.setProperty("RRIMER_DELTA_GC_WEIGHT", Double.toString(0.5));
		configFile.setProperty("PRIMER_DELTA_LENGTH_WEIGHT", Double.toString(0.5));
		configFile.setProperty("PRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT", Double.toString(1));
		configFile.setProperty("SELF_ALIGNMENT_WEIGHT", Double.toString(0.1));
		configFile.setProperty("SELF_END_ALIGNMENT_WEIGHT", Double.toString(0.2));
		configFile.setProperty("PAIR_ALIGNMENT_WEIGHT", Double.toString(0.1));
		configFile.setProperty("PAIR_ALIGNMENT_WEIGHT", Double.toString(0.2));
		configFile.setProperty("PRIMER_FALSE_POSITIVES_WEIGHT", Double.toString(1));
		
		return configFile;
	}
	
	private static Properties parseConfigFile(InputStream stream) throws InvalidPropertiesFormatException, IOException{
		Properties configFile = new Properties();
		configFile.loadFromXML(stream);
		
		return configFile;
	}
}
