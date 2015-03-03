package primerDesign.Test;

import java.io.IOException;

import primerDesign.util.SimpleTimer;

public class PerformanceComparison {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		SimpleTimer timer = new SimpleTimer();
		EnhancedSuffixArrayIntOld array = EnhancedSuffixArrayIntOld.deserialize(args[1]);
		System.out.println("Deserializing " + array.getClass().toString() + " took: " + timer.getTimeString());
//		String[] patterns = {"A", "AT", "ATG", "ATGC", "ATGCA", "ATGCAT", "ATGCATG", "ATGCATGC", "ATGCATGCA", "ATGCATGCAT", "ATGCATGCATG", "ATGCATGCATGC"};
//		for(int i=0; i<patterns.length; i++){
//			System.out.println("Querying ESAintOpt for: " + patterns[i]);
//			System.out.println(array.findMatchPositions(patterns[i]).length + " matches - took " + timer.getTimeString());
//		}
		
		int queries = Integer.parseInt(args[0]);
		String query = "ATGCATGCATGC";
		System.out.print("Performing " + queries + " queries of length " + query.length());
		timer.resetTimer();
		for(int i=0; i<queries; i++){
			array.findMatchPositions(query);
		}
		System.out.println(" - took " + timer.getTimeString());
	}

}
