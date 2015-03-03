/**
 * 
 */
package primerDesign.web;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;

/**
 * Class for handling the 3PD servelet config file.
 * 
 * @author Sebastian Fršhler
 *
 */
public class PrimerDesignWebProperites {
	public static void main(String[] args) throws FileNotFoundException, IOException{
		Properties configFile = PrimerDesignWebProperites.getDefaultProperties();
		
		System.out.println("Writing properties");
		configFile.storeToXML(new FileOutputStream("/export/Sebastian/PrimerDesign/deploy/3CPrimerDesign-ServletConfig.xml"), "Properites for the web interface");
	}
	
	/**
	 * Returns the default configuration for 3PD.
	 * 
	 * @return the default configuration for 3PD
	 */
	public static Properties getDefaultProperties(){
		Properties configFile = new Properties();
		
		configFile.setProperty("SEQUENCE_INDEX_PATH", "/var/tmp");
		configFile.setProperty("C.elegans-IndexFile", "C.elegans.dummy");
		configFile.setProperty("P.pacificus-IndexFile", "P.pacificus.dummy");
		configFile.setProperty("M.musculus-IndexFile", "Mmusculus/Mmusculus.fa.MultiSeqMemoryMappedESAIndex.esaidx");
		configFile.setProperty("MAX_NUM_PICKING_THREADS", "all");
		configFile.setProperty("ERROR_REPORT_EMAIL", "froehler@tuebingen.mpg.de");
		configFile.setProperty("THIS_SERVICE_SENDER_EMAIL", "froehler@tuebingen.mpg.de");
		configFile.setProperty("SMTP_SERVER", "mailhost.tuebingen.mpg.de");
		
		return configFile;
	}
	
	/**
	 * Parses a 3PD config file.
	 * 
	 * @param stream the stream to the 3PD config file
	 * 
	 * @return the parameters specified in the config file
	 * 
	 * @throws InvalidPropertiesFormatException
	 * @throws IOException
	 */
	public static Properties parseConfigFile(InputStream stream) throws InvalidPropertiesFormatException, IOException{
		Properties configFile = new Properties();
		configFile.loadFromXML(stream);
		
		return configFile;
	}
}
