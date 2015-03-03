package primerDesign.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A parser for a list of restriction enzymes.
 * 
 * @author Sebastian Fršhler
 *
 */
public class RestrictionEnzymeListParser {
	/**
	 * Parses a list of restriction enzyme.
	 * 
	 * @param filename the file with the restrition enzymes
	 * 
	 * @return a list of restriction enzymes from file 'filename'
	 * 
	 * @throws NumberFormatException
	 * @throws IllegalAlphabetException
	 * @throws IllegalSymbolException
	 * @throws IOException
	 */
	public static ArrayList<RestrictionEnzyme> parseEnzymesList(File filename) throws NumberFormatException, IllegalAlphabetException, IllegalSymbolException, IOException{
		ArrayList<RestrictionEnzyme> enzymes = new ArrayList<RestrictionEnzyme>();
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String line;
		String[] values;
		while((line = reader.readLine()) != null){
			if(line.startsWith("#")) continue;
			values = line.split("\t");
			if(values.length >= 4){
				enzymes.add(new RestrictionEnzyme(values[0], DNATools.createDNA(values[1]), Integer.parseInt(values[2]), Integer.parseInt(values[3])));
			}
		}
		return enzymes;
	}		
	
	public static void main(String[] args) throws NumberFormatException, IllegalAlphabetException, IllegalSymbolException, IOException{
		ArrayList<RestrictionEnzyme> enzymes = RestrictionEnzymeListParser.parseEnzymesList(new File("/Users/froehler/Desktop/Project/HOX cluster analysis/3C-Candidate_Enzymes-woDamDcmEcoKIEcoBI.txt"));
		for(RestrictionEnzyme enzyme : enzymes){
			System.out.println(enzyme.toString());
		}
	}
}
