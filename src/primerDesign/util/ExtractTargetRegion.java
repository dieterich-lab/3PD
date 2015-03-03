package primerDesign.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SoftMaskedAlphabet;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

/**
 * Class for extracting target regions for 3PD.
 * 
 * @author Sebastian Fršhler
 *
 */
public class ExtractTargetRegion {

	public static void main(String[] args) throws IOException, NoSuchElementException, BioException {
		
		String file = args[0];
		String contigName = args[1];
		int start = Integer.parseInt(args[2]);
		int stop = Integer.parseInt(args[3]);

		BufferedReader reader = new BufferedReader(new FileReader(file));
		SoftMaskedAlphabet alphabet = SoftMaskedAlphabet.getInstance((FiniteAlphabet) AlphabetManager.alphabetForName("DNA"));
		//SequenceIterator stream = SeqIOTools.readFasta(reader, alphabet.getTokenization("token"));
		RichSequenceIterator iter = RichSequence.IOTools.readFasta(reader, alphabet.getTokenization("token"), null);
		while(iter.hasNext()){
			Sequence seq = iter.nextSequence();
			if(seq.getName().matches(".*:" + contigName + ":.*")){
				if(seq.length() < stop) throw new IllegalArgumentException("Sequence has only length " + seq.length());
				if(start > stop) throw new IllegalArgumentException("Start has be < stop");
				if(start < 0) throw new IllegalArgumentException("Start has to be positive");
				System.out.println(">" + contigName + ":" + start + ":" + stop);
				System.out.println(seq.subStr(start, stop).replaceAll("[a-z]", "N"));
				break;
			}
		}
	}
}
