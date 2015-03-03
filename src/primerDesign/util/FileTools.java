/**
 * 
 */
package primerDesign.util;

import java.io.File;

/**
 * @author Sebastian Fršhler
 *
 */
public class FileTools {
	public static final String PATH_SEPARATOR = File.separator;
	
	public static String extractPath(File file){
		String[] values = file.getAbsolutePath().split(PATH_SEPARATOR);
		StringBuffer buffy = new StringBuffer();
		for(int i=1; i<values.length - 1; i++){
			buffy.append(PATH_SEPARATOR);
			buffy.append(values[i]);
		}
		return buffy.toString();
	}
	
	public static String extractFilename(File filename){
		String[] values = filename.getAbsolutePath().split(PATH_SEPARATOR);

		return values[values.length - 1];
	}
}
