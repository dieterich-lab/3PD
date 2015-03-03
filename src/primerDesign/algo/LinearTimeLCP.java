/**
 * 
 */
package primerDesign.algo;


/**
 * Implements a linear-time longest common prefix calculation.
 * 
 * @author Sebastian Fršhler
 *
 */
public class LinearTimeLCP {
	/**
	 * Comoutes the longest common prefix within 'text'.
	 * 
	 * @param text the sequence
	 * @param suffixArray the suffix arrax of 'text'
	 * 
	 * @return the longest common prefix within 'text'
	 */
	public static int[] getLCP(char[] text, int[] suffixArray){
		int length = suffixArray.length;
		
		int[] rank = new int[length];
		int[] lcp =new int[length];
		
		for(int i=0; i<length; i++){
			rank[suffixArray[i]] = i;
		}
		
		int h=0;
		int j;
		int k;
		lcp[0] = 0;
		for(int i=0; i<length; i++){
			k = rank[i];
			if(k == 1){
				lcp[k] = 0;
			}
			else if(k > 1){				
				j = suffixArray[k - 1];
				if(i+h < text.length && j+h < text.length){
					while(i+h < length - 1 && j+h < length - 1 && text[i+h] == text[j+h]){
						h++;
					}
				}
				lcp[k] = h;
				if(h > 0){
					h--;
				}
			}
		}
		return lcp;
	}
	
	public static void main(String[] args){
		char[] text = "acaaacatat".toCharArray();
		KoAluruSuffixSort sort = new KoAluruSuffixSort();
		int[] suffixArray = sort.getSuffixArray(text);
		char[] temp = text;
//		char[] temp = new char[text.length + 1];
//		System.arraycopy(text, 0, temp, 0, text.length);
//		temp[text.length] = '$';
		//int[] suffixArray = new int[]{0, 11, 3,4,1,5,7,9,2,6,8,10};
		//int[] suffixArray = new int[]{10, 2,3,0,4,6,8,1,5,7,9};
		int[] lcp = LinearTimeLCP.getLCP(temp, suffixArray);
		
		System.out.println(printArray(lcp, temp, suffixArray));
		System.out.println(lcpToString(lcp));
		System.out.println("0 0 2 1 3 1 2 0 2 0 1");
	}
	
	private static String printArray(int[] array, char[] text, int[] suffixArray){
		StringBuffer buffy = new StringBuffer();
		for(int i=0; i<array.length; i++){
			buffy.append(suffixArray[i] + " " + new String(text).substring(suffixArray[i]) + " " + array[i] + "\n");
		}
		return buffy.toString();
	}
	
	private static String lcpToString(int[] lcp){
		StringBuffer buffy = new StringBuffer();
		for(int i=0; i<lcp.length; i++){
			buffy.append(lcp[i]);
			if(i<lcp.length-1) buffy.append(" ");
		}
		return buffy.toString();
	}
}
