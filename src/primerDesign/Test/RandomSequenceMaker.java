/**
 * 
 */
package primerDesign.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import primerDesign.util.SeqTools;

/**
 * @author froehler
 *
 */
public class RandomSequenceMaker {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		int startSize = Integer.parseInt(args[0]);
		int steps = Integer.parseInt(args[1]);
		int item = 0;
		
		for(int j=0; j<steps; j++){
			StringBuffer buffy = new StringBuffer();
			for(int i=0; i<startSize; i++){
				buffy.append(SeqTools.getRandomBase());
			}
			BufferedWriter writer = new BufferedWriter(new FileWriter(startSize/1000 + "kb.dna"));
			writer.write(">Testheader_" + item + "\n");
			writer.write(buffy.toString());
			item++;
			startSize *= 2;
		}
	}
}
