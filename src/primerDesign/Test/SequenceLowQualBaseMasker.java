/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author froehler
 *
 */
public class SequenceLowQualBaseMasker {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		SlimFastQParser parser = new SlimFastQParser(new File("/Users/froehler/TEMP/BAC-Ppa50-C09/Contig46.fastq"));
		SimpleContigWithQuality contig;
		
		int start = 1197940; // ceh-13
		int end = 1246582; // mab-5
		int padding = 10000; // +/- 10kb
		
		int qualThreshold = 40;
		char maskingChar = 'N';
		
		start -= padding;
		end += padding;
		
		while(parser.hasNextContig()){
			contig = parser.parseNextContigIgnoreCase();
			if(contig.getID().matches(".*Contig46.*")){
				char[] ctg = contig.getSequence();
				for(int i=start; i<=end; i++){
					if(contig.getQuality()[i] < qualThreshold){
						ctg[i] = maskingChar;
					}
				}
				System.out.println(new String(ctg).substring(start, end + 1));
			}
		}
	}

}
