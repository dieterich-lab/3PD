package primerDesign.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

import junit.framework.TestCase;

public class EnhancedSuffixArrayIntTest extends TestCase {

	public void testLCPTable() throws IOException{
		EnhancedSuffixArrayInt array = new EnhancedSuffixArrayInt(readSequence());
		array.createIndex(Integer.MAX_VALUE, true);
		for(int i=0; i<array.lcptab_compara.length; i++){
			assertEquals(array.lcptab_compara[i], array.getLcpTab(i));
		}
	}
	
	public void testRecompChildTabUp() throws IOException {
		EnhancedSuffixArrayInt array = new EnhancedSuffixArrayInt(readSequence());
		array.createIndex(Integer.MAX_VALUE, true);
		for(int i=0; i<array.childtab_compara_up.length; i++){
			assertEquals(array.childtab_compara_up[i], array.recompChildTabUp(i));
		}
	}

	public void testRecompChildTabDown() throws IOException {
		EnhancedSuffixArrayInt array = new EnhancedSuffixArrayInt(readSequence());
		array.createIndex(Integer.MAX_VALUE, true);
		for(int i=0; i<array.childtab_compara_down.length; i++){
			assertEquals(array.childtab_compara_down[i], array.recompChildTabDown(i));
		}
	}

	public void testRecompChildTabNext() throws IOException {
		EnhancedSuffixArrayInt array = new EnhancedSuffixArrayInt(readSequence());
		array.createIndex(Integer.MAX_VALUE, true);
		for(int i=0; i<array.childtab_compara_next.length; i++){
			assertEquals(array.childtab_compara_next[i], array.recompChildTabNext(i));
		}
	}
	
	public void testChildTable() throws IOException {
		EnhancedSuffixArrayInt array = new EnhancedSuffixArrayInt(readSequence());
		array.createIndex(Integer.MAX_VALUE, true);
		for(int i=0; i<array.childtab_compara_next.length; i++){
			if(i>0 && array.containsUpIndex(i-1)) assertEquals(array.childtab_compara_up[i], array.getChildTabUP(i));
			if(array.containsDownIndex(i)) assertEquals(array.childtab_compara_down[i], array.getChildTabDown(i));
			if(array.containsNextIndex(i)) assertEquals(array.childtab_compara_next[i], array.getChildTabNext(i));
		}
	}
	
	private String readSequence() throws IOException{
		Pattern fastaStart = Pattern.compile(">.*");
				String sequence;
				{
					BufferedReader reader = new BufferedReader(new FileReader("/export/Sebastian/PrimerDesign/Testsequenzen/../src/yeast-chrXII.dna"));
					//BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(path + sequences[i]), "US-ASCII"));
					StringBuffer seq = new StringBuffer();
					String line;
					while((line = reader.readLine()) != null){
						if(!fastaStart.matcher(line).matches()){
							seq.append(line.trim());
						}
					}
					sequence = seq.toString().toUpperCase();
					seq = null;  // clear string buffer -> save 2n space!
				}
		return sequence;
		//return "ACAAACATAT";
	}
}
