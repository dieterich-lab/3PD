/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import primerDesign.util.SeqTools;

/**
 * @author froehler
 *
 */
public class SequencemaskingScreener {
	public static void main(String[] args) throws IOException{
		SlimFastQParser parser = new SlimFastQParser(new File("/Users/froehler/TEMP/BAC-Ppa50-C09/Contig46.fastq")); // args[0]));
		while(parser.hasNextContig()){
			SimpleContigWithQuality contig = parser.parseNextContigIgnoreCase();

			SimpleContigWithQuality ctg = new SimpleContigWithQuality("Contig46-HOX", Arrays.copyOfRange(contig.getSequence(), 1197657, 1267090), Arrays.copyOfRange(contig.getQuality(), 1197657, 1267090));
			ctg.setSequence(SeqTools.maskSequence(ctg.getSequence(), ctg.getQuality(), (byte) 40));
			
//			BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/froehler/TEMP/BAC-Ppa50-C09/Contig46-HOX.fastq"));
//			writer.write(ctg.toFastaString());
//			writer.close();
//			writer = new BufferedWriter(new FileWriter("/Users/froehler/TEMP/BAC-Ppa50-C09/Contig46-HOX.fa"));
//			writer.write(">" + ctg.getName() + "\n" + new String(ctg.getSequence()));
//			writer.close();
			
			boolean mask = false;
			for(int i=0; i<ctg.getSequenceLength(); i++){
				if(!mask && ctg.getSequence()[i] == 'N'){
					mask = true;
					System.out.print(i + " - ");
				}
				else if(mask && ctg.getSequence()[i] != 'N'){
					mask = false;
					System.out.println(i-1);
				}
			}
		}
	}
}
