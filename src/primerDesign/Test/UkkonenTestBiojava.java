package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

import org.biojava.bio.symbol.UkkonenSuffixTree;

import primerDesign.util.SimpleTimer;

public class UkkonenTestBiojava extends UkkonenSuffixTree{
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		SimpleTimer timer = new SimpleTimer();
		timer.startTimer();
		
		String[] sequences = {"100kb.dna", "200kb.dna", "400kb.dna", "800kb.dna", "../src/yeast-chrXII.dna", "Ppa/Ppacificus.fa", "Cel/Celegans.fa"};
		String path = "/export/Sebastian/PrimerDesign/Testsequenzen/";
		int maxWordLength = 30;
		for(int i=0; i<sequences.length; i++){
			System.out.println("Reading sequence " + sequences[i]);
			BufferedReader reader = new BufferedReader(new FileReader(path + sequences[i]));
			StringBuffer sequence = new StringBuffer();
			String line;
			Pattern fastaStart = Pattern.compile(">.*");
			while((line = reader.readLine()) != null){
				if(!fastaStart.matcher(line).matches()){
					sequence.append(line.trim());
				}
			}
			System.out.println("Reading sequence " + sequences[i] + " (length: " + sequence.length() + ") took " + timer.getTimeString());
			
			for(int j=10; j<= maxWordLength; j += 10){
				System.out.println("Creating Index for sequence " + sequences[i] + " maxWordLength " + j);
				
				System.out.println("Constructing Ukkonen suffix tree");
				UkkonenSuffixTree tree = new UkkonenSuffixTree(sequence.toString());
				System.out.println("Suffix tree construction complete - " + timer.getTimeString());
			}
		}
	}
}
