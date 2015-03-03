/**
 * 
 */
package primerDesign.util;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;

/**
 * Class returning general information about a sequence.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SeqInfo {

	public static void main(String[] args) throws IOException {
		SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
		SimpleContig contig;
		char[] seq;
		int gc = 0;
		int masked = 0;
		int length = 0;
		int contigs = 0;
		while(parser.hasNextContig()){
			contig = parser.parseNextContigIgnoreCase();
			seq = contig.getSequence();
			contigs++;
			for(int i=0; i<seq.length; i++){
				if(seq[i] == 'N') masked++;
				else if(seq[i] == 'C' || seq[i] == 'G') gc++;
				length++;
			}
		}
		NumberFormat format = NumberFormat.getInstance();
		System.out.println("Sequence length: " + format.format(length) + " bp");
		System.out.println("Contigs: " + format.format(contigs));
		System.out.println("GC content: " + format.format((double)gc/length*100) + "%");
		System.out.println("Masked basepairs: " + format.format(masked) + " bp (" + format.format((double)masked/length*100) + "%)");
	}
}
