package primerDesign.Test;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import cern.colt.list.ObjectArrayList;

public class StringIndex implements DNASequenceIndex {
	private int maxWordSize;
	private boolean includesScanRegion;
	private String sequence;
	
	public void createIndex(String sequence, int maxWordSize, boolean includesScanRegion) {
		this.sequence = sequence;
		this.maxWordSize = maxWordSize;
		this.includesScanRegion = includesScanRegion;
	}

	public int getMaxWordSize() {
		return this.maxWordSize;
	}

	public String getSequence() {
		return this.sequence;
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	public Integer[] searchMatchPositionsInIndex(String searchString) {
		ObjectArrayList list = new ObjectArrayList();
		searchString = searchString.toUpperCase();
		int position = this.sequence.indexOf(searchString);
		while(position != -1){
			list.add(position);
			position = this.sequence.indexOf(searchString, position+1);
		}
		String reverseString = SeqTools.revcompDNA(searchString.toCharArray());
		position = this.sequence.indexOf(reverseString);
		while(position != -1){
			list.add(position);
			position = this.sequence.indexOf(reverseString, position+1);
		}
		Integer[] result = new Integer[list.size()];
		for(int i=0; i<list.size(); i++){
			result[i] = (Integer) list.get(i);
		}
		return result;
	}

	public int searchNbMatchesInIndex(String searchString) {
		return searchMatchPositionsInIndex(searchString).length;
	}

	public static void main(String[] args) throws IOException{
		SimpleTimer timer = new SimpleTimer();
		
		String file = args[0];
		String query = args[1];
		int numQueries = Integer.parseInt(args[2]);
		
		//Pattern fastaStart = Pattern.compile(">.*");
		System.out.println("Reading sequence " + file);
		String sequence;
		{
			{
				//BufferedReader reader = new BufferedReader(new FileReader(file)); //path + sequences[i]));
				//BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(path + sequences[i]), "US-ASCII"));
				StringBuffer seq = new StringBuffer();
				//String line;
				SlimFastaParser parser = new SlimFastaParser(new File(file));
				while(parser.hasNextContig()){
					seq.append(new String(parser.parseNextContigIgnoreCase().getSequence()));
					if(parser.hasNextContig()) seq.append("N");
				}
//				while((line = reader.readLine()) != null){
//					if(!fastaStart.matcher(line).matches()){
//						seq.append(line.trim());
//					}
//					else seq.append("N");
//				}
				sequence = seq.toString().toUpperCase();
				seq = null;  // clear string buffer -> save 2n space!
			}
			System.gc();
			System.out.println("Reading sequence " + file + " (length: " + sequence.length() + ") took " + timer.getTimeString());
		}
		
		StringIndex index = new StringIndex();
		index.createIndex(sequence, sequence.length(), true);
		
		System.out.print("Querying index " + numQueries + " times for string " + query);
		Integer[] result = {};
		for(int i=0; i<numQueries; i++){
			result = index.searchMatchPositionsInIndex(query);
		}
		System.out.println(" " + result.length + " matches - done in " + timer.getTimeString());
	}
}
