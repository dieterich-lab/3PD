package primerDesign.Test;

import java.text.NumberFormat;

import primerDesign.util.SeqTools;
import primerDesign.util.SimpleTimer;



/**
 * Implements a simple dynamic programming (w.o. subst. matrix) between two strings represented as char[]s.
 * 
 * The substrings to be aligned can be specified by the respective start end end positons.
 * 
 * @author froehler
 *
 */
public class DP {
	private static int[][] matrix;
	private static final int insertCost = 1;
	private static final int delCost = 1;
	private static final int matchCost = 0;
	private static final int mismatchCost = 1;
	
	/**
	 * Computes the dynamic programming between a pattern and a text of positions patStart to patEnd of the pattern with
	 * positions textStart to textEnd of the text.
	 * 
	 * Pattern and string is accessed in a read-only fashion!
	 * 
	 * @param pattern the pattern
	 * @param patStart the start position in the pattern
	 * @param patEnd the end position in the pattern
	 * @param text the text
	 * @param textStart the start position in the text
	 * @param textEnd the end position in the text
	 * 
	 * @return the DP value of a global alignment of pattern(patStart...patEnd) with text(textStart...textEnd)
	 */
	public static int getDP(final char[] pattern, int patStart, int patEnd, final char[] text, int textStart, int textEnd){
		assert(patStart>=0 && patStart < pattern.length && patEnd >= 0 && patEnd < pattern.length && patStart <= patEnd);
		assert(textStart>=0 && textStart < text.length && textEnd >= 0 && textEnd < text.length && textStart <= textEnd);
		int patternEntries = patEnd - patStart + 2;
		int textEntries = textEnd - textStart + 2;
		
		// set up DP matrix
		DP.matrix = new int[patternEntries][textEntries];
		
		// init first row and column of DP matrix
		for(int i=0; i<patternEntries; i++){
			DP.matrix[i][0] = i;
		}
		for(int i=0; i<textEntries; i++){
			DP.matrix[0][i] = i;
		}
		
		for(int i=1; i<patternEntries; i++){
			for(int j=1; j<textEntries; j++){
				DP.matrix[i][j] = DP.min(DP.matrix[i-1][j] + DP.delCost, DP.matrix[i][j-1] + DP.insertCost, DP.matrix[i-1][j-1] + compare(pattern[patStart+i-1], text[textStart+j-1]));
			}
		}
		
		return DP.matrix[patternEntries - 1][textEntries - 1];
	}
	
	public static int getRevcompDP(final char[] pattern, int patStart, int patEnd, final char[] text, int textStart, int textEnd){
		assert(patStart>=0 && patStart < pattern.length && patEnd >= 0 && patEnd < pattern.length && patStart <= patEnd);
		assert(textStart>=0 && textStart < text.length && textEnd >= 0 && textEnd < text.length && textStart <= textEnd);
		int patternEntries = patEnd - patStart + 2;
		int textEntries = textEnd - textStart + 2;
		
		// set up DP matrix
		DP.matrix = new int[patternEntries][textEntries];
		
		// init first row and column of DP matrix
		for(int i=0; i<patternEntries; i++){
			DP.matrix[i][0] = i;
		}
		for(int i=0; i<textEntries; i++){
			DP.matrix[0][i] = i;
		}
		
		for(int i=1; i<patternEntries; i++){
			for(int j=1; j<textEntries; j++){
				DP.matrix[i][j] = DP.min(DP.matrix[i-1][j] + DP.delCost, DP.matrix[i][j-1] + DP.insertCost, DP.matrix[i-1][j-1] + compareRevComp(pattern[patStart+i-1], text[textEnd-(j-1)]));
			}
		}
		
		return DP.matrix[patternEntries - 1][textEntries - 1];
	}
	
	public static int getCompDP(final char[] pattern, int patStart, int patEnd, final char[] text, int textStart, int textEnd){
		assert(patStart>=0 && patStart < pattern.length && patEnd >= 0 && patEnd < pattern.length && patStart <= patEnd);
		assert(textStart>=0 && textStart < text.length && textEnd >= 0 && textEnd < text.length && textStart <= textEnd);
		int patternEntries = patEnd - patStart + 2;
		int textEntries = textEnd - textStart + 2;
		
		// set up DP matrix
		DP.matrix = new int[patternEntries][textEntries];
		
		// init first row and column of DP matrix
		for(int i=0; i<patternEntries; i++){
			DP.matrix[i][0] = i;
		}
		for(int i=0; i<textEntries; i++){
			DP.matrix[0][i] = i;
		}
		
		for(int i=1; i<patternEntries; i++){
			for(int j=1; j<textEntries; j++){
				DP.matrix[i][j] = DP.min(DP.matrix[i-1][j] + DP.delCost, DP.matrix[i][j-1] + DP.insertCost, DP.matrix[i-1][j-1] + compareComp(pattern[patStart+i-1], text[textStart+j-1]));
			}
		}
		
		return DP.matrix[patternEntries - 1][textEntries - 1];
	}
	
	private static int min(int a, int b, int c){
		return Math.min(a, Math.min(b, c));
	}
	
	private static int compare(char a, char b){
		if(a == b) return matchCost;
		else return mismatchCost;
	}
	
	private static int compareRevComp(char a, char b){
		if(a == SeqTools.revcompChar(b)) return matchCost;
		else return mismatchCost;
	}
	
	private static int compareComp(char a, char b){
		if(a == SeqTools.complementChar(b)) return matchCost;
		else return mismatchCost;
	}
	
	private static void printMatrix(){
		for(int i=0; i<DP.matrix.length; i++){
			for(int j=0; j<DP.matrix[0].length; j++){
				System.out.print(DP.matrix[i][j]);
			}
			System.out.println();
		}
	}
	
	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
		NumberFormat format = NumberFormat.getInstance();
		
		int iter = 1000000;
		char[] pattern = "atgcatgcatgcatgcatgcatgcatgcatgcatgc".toCharArray();
		char[] text = "atgcatgcatgcatgcatgcatgcatgcatgcatgc".toCharArray();
		
		for(int j=1; j<10; j++){
			System.out.print("Computing " + format.format(iter) + " forward-forward DP matrices of size " + pattern.length/j + "x" + text.length/j);
			for(int i=0; i<iter; i++){
				DP.getDP(pattern, 0, 35/j, text, 0, 35/j);
				//System.out.println(DP.getDP(pattern, 0, 35, "atgcatgcatgcatgcatgcatgcatgcatgcatgc".toCharArray(), 0, 35));
				//DP.printMatrix();
			}
			System.out.println(" - forward done in " + timer.getTimeString());
			
			System.out.print("Computing " + format.format(iter) + " forward-reverse DP matrices of size " + pattern.length/j + "x" + text.length/j);
			for(int i=0; i<iter; i++){
				DP.getRevcompDP(pattern, 0, 35/j, text, 0, 35/j);
				//System.out.println(DP.getDP(pattern, 0, 35, "atgcatgcatgcatgcatgcatgcatgcatgcatgc".toCharArray(), 0, 35));
				//DP.printMatrix();
			}
			System.out.println(" - reverse done in " + timer.getTimeString());
		}
//		System.out.println(DP.getDP("aca".toCharArray(), 0, 2, "acaacaca".toCharArray(), 1, 3));
//		printMatrix();
//		System.out.println(DP.getRevcompDP("atgc".toCharArray(), 0, 3, "gcat".toCharArray(), 0, 3));
//		printMatrix();
	}
}
