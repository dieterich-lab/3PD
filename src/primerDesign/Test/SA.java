/**
 * This class implements a simple suffix array using binary search as proposed by Manber & Myers (simple search!).
 * 
 * Reference: Manber, Myers: Sufﬁx arrays: A new method for on-line string searches; Proceedings of the first annual ACM-SIAM symposium on Discrete algorithms 1990:319-327 
 * 
 */
package primerDesign.Test;

import primerDesign.util.MyExtendedMath;

/**
 * @author froehler
 *
 */
public class SA {
	private int[] suffix_array;
    private char[] text;
    private int textLength;
    
    public SA(char[] text, int[] sa){
    	this.text = text;
    	this.textLength = text.length;
    	this.suffix_array = sa;
    }
    
    public int binarySearch(char[] query){
    	int left = 0;
    	int right = this.suffix_array.length - 1;
    	int mid;
    	
    	if(query.length == 0 || compare(query, left, query.length) < 0 || compare(query, right, query.length) > 0) return -1;
    	else{
    		while(right - left > 1){
    			mid = MyExtendedMath.round(((double)(right + left)) / 2);
    			if(compare(query, mid, query.length) <= 0){
    				right = mid;
    			}
    			else{
    				left = mid;
    			}
    		}
    	}
    	return this.suffix_array[right];
    }
	
	private int compare(char[] query, int position, int length){
		int result = 0;
		if(this.suffix_array[position] + length > textLength){
			result = +1;
			length = textLength - this.suffix_array[position];
		}
		for(int i=0; i<length; i++){
			if(this.text[this.suffix_array[position] + i] < query[i]) return +1;
			else if(this.text[this.suffix_array[position] + i] > query[i]) return -1;
		}
		return result;
	}
	
	private String get(int index){
		if(this.suffix_array[index] < this.text.length) return "\t" + this.suffix_array[index] + "\t" + new String(this.text).substring(this.suffix_array[index]);
		else return "$";
	}
	
  
	public static void main(String[] args){
	  char[] text = "abracadabra".toCharArray();
	  KoAluruSuffixSort sorter = new KoAluruSuffixSort();
	  int[] sa = sorter.getSuffixArray(text);
	  
	  SA array = new SA(text, sa);
	  
	  System.out.println("SA-idx\tTxt-idx\tString");
	  for(int i=0; i<=text.length; i++){
		  System.out.println(i + " " + array.get(i));
	  }
	  System.out.println();
	  
	  System.out.println("0: " + array.binarySearch("abracadabra".toCharArray()));
	  System.out.println("1: " + array.binarySearch("bracadabra".toCharArray()));
	  System.out.println("2: " + array.binarySearch("racadabra".toCharArray()));
	  System.out.println("3: " + array.binarySearch("acadabra".toCharArray()));
	  System.out.println("4: " + array.binarySearch("cadabra".toCharArray()));
	  System.out.println("5: " + array.binarySearch("adabra".toCharArray()));
	  System.out.println("6: " + array.binarySearch("dabra".toCharArray()));
	  System.out.println("7:" + array.binarySearch("abra".toCharArray()));
	  System.out.println("8: " + array.binarySearch("bra".toCharArray()));
	  System.out.println("9: " + array.binarySearch("ra".toCharArray()));
	  System.out.println("10: " + array.binarySearch("a".toCharArray()));
	  System.out.println("NO: " + array.binarySearch("".toCharArray()));
    }
}
