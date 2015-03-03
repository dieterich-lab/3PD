/**
 * 
 */
package primerDesign.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Arrays;

import primerDesign.util.SimpleTimer;


/**
 * This class implements the suffix sort as proposed by Ko and Aluru 2003.
 * 
 * Reference: P. Ko, A. Aluru:  Space Efficient Linear Time Construction of Suffix Arrays, 
 * Proceedings of the 14th Annual Symposium on Combinatorial Pattern Matching (CPM'03), 
 * Lecture Notes in Computer Science, Vol. 2676, 200-210, 2003
 * 
 * Code is adapted from an efficient implementation thereof by:
 * Sunglim Lee and Kunsoo Park, Efficient Implementations of Suffix Array Construction Algorithms, Proceedings of the 15th Australasian Workshop on Combinatorial Algorithms, pp. 64-72, 2004.
 * 
 * Memory consumption is claimed to be <= 13n
 * 
 * @author froehler
 *
 */
public class KoAluruSuffixSort {
	private static final int CHAR_SIZE = Character.SIZE / 8;
	private static final int INT_SIZE = Integer.SIZE / 8;
	private static final int TERMIN =  0;					// used as unique termination character '$'
	private static final int MAX_INPUT_ALPHABET = (1 << (CHAR_SIZE * 8 - 1)) - 1;		// including the termination character '$'
	private static final boolean ADD_TERMIN = true;			// whether to add the termination character '$'
	private static final boolean ALWAYS_ENCODE_S = false;
	private static final boolean ALWAYS_ENCODE_L = false;
	
	// debug and tuning options
	private static final boolean PRINT_MEMORY_USAGE = false;
	private static final boolean DEBUG = false;
	
	public int[] getSuffixArray(char[] text){
		// pre-process suffix array, shift values by one to account for termination symbol	
		// compute suffix array
		int[] suffixArray;
		if(ADD_TERMIN){
			suffixArray = new int[text.length + 1];
		}
		else{
			suffixArray = new int[text.length];
		}
		int[] temp = sort(shiftSuffixes(text), suffixArray.length, MAX_INPUT_ALPHABET, suffixArray);
		System.arraycopy(temp, 0, suffixArray, 0, suffixArray.length);
		return suffixArray;
	}
	
	private int[] recursivesort(int[] sufenc, int textlength, int iTypeSL, int iTdashMaxAlphabet){
		int i, j;

		final int[] tDash = new int[iTypeSL+1];
		// go through sufenc to type S(L) substrings
		j = 0;
		for (i = 0; i < textlength; i++) 
		{
			if (sufenc[i] != -1) // is a typeSL suffix
			{
				assert(j < iTypeSL);
				assert(sufenc[i] <= Integer.MAX_VALUE && sufenc[i] >= Integer.MIN_VALUE);
				tDash[j++] = sufenc[i];
			}
		}
		tDash[j] = TERMIN;		// assign termination character '$'

		// Tdash generated
		// </sufenc>

		// <SAdash>
		int[] saDash = sufenc;

		// Sort Tdash recursivly (KA p7)
		// Acceleration: skip recursion if iBucketNumberMaxAlphabet >= TdashSA
		if (iTdashMaxAlphabet+1 < iTypeSL+1)	
		{
			// recursion - saDash is implicitely updated, NO return value needs to be catched!
			saDash = sort(tDash, iTypeSL+1, iTdashMaxAlphabet, saDash);
		} 
		else	// generate suffix array directly
		{
			for (i = 0;  i < iTypeSL+1;  i++) saDash[tDash[i]] = i;
		}
		assert(saDash[0] == iTypeSL);		// the first entry should point the end of Tdash: '$'

		return saDash;
	}

	private int[] sort(final int[] text, final int textLength, int iTdashMaxAlphabet, int[] suffixArray){
		// init variables
		int x = Integer.MIN_VALUE, y = Integer.MIN_VALUE;
		int i = 0, j = 0, k = 0;
		
		assert((char)(1 << (INT_SIZE * 8 - 1)) <= Integer.MAX_VALUE && (1 << (INT_SIZE * 8 - 1)) >= Integer.MIN_VALUE);
		final int signBit = 1 << (INT_SIZE * 8 - 1);
		final int alphabetSize = iTdashMaxAlphabet + 1;
		
		if(suffixArray.length < textLength * INT_SIZE){
			int[] temp = new int[textLength * INT_SIZE];
			System.arraycopy(suffixArray, 0, temp, 0, suffixArray.length);
			suffixArray = temp;
		}
		//suffixArray = new int[textLength * INT_SIZE];
		final int[] sufenc;
		
		assert(textLength > 0);					// should have at least one element
		assert(text[textLength-1] == TERMIN);	 // end with '$'
		
		// check termination condition
		if (textLength <= 3)
		{
			switch (textLength)
			{
			case 1:
				suffixArray[0] = 0;
				break;
			case 2:
				suffixArray[0] = 1;
				suffixArray[1] = 0;
				break;
			case 3:
				suffixArray[0] = 2;
				if (text[0] >= text[1])
				{
					suffixArray[1] = 0;
					suffixArray[2] = 1;
				}
				else
				{
					suffixArray[1] = 1;
					suffixArray[2] = 0;
				}
				break;
			}

			return suffixArray;
		}
		
		// type suffixes into S and L (KA p3), use sign bit of input T
		// set S, L flag: S-negative, L-positive
		// set params: iTypeSL
		// iTypeSL first contains iTypeS
		int iTypeSL = 1;	// number of Type S and Type L Suffixes, start with 1 for last $
		j = 0;
		for (i = 0; i < textLength-1; i++)
		{
			if (text[i] > text[i+1])	// type L
			{
				j = i + 1;
			}
			else if (text[i] < text[i+1])	// type S
			{
				iTypeSL += i + 1 - j;
				for (; j < i + 1; j++)
				{
					text[j] |= signBit;
				}
			}
		}
		
		assert(j == textLength-1);	// last $ is both S and L
		
		// Encode and recursively sort type S or L substrings
		int signBitSorL;
		if(ALWAYS_ENCODE_S){
			signBitSorL = signBit;
			text[textLength - 1] |= signBit;
		} 
		else if(ALWAYS_ENCODE_L){
			signBitSorL = 0;
			iTypeSL = textLength + 1 - iTypeSL;
		}
		else{
			if (iTypeSL <= textLength + 1 - iTypeSL)
			{
				signBitSorL = signBit;	// indicate we will encode type S suffixes
				text[textLength-1] |= signBit;
			}
			else
			{
				signBitSorL = 0;	// indicate we will encode type L suffixes
				iTypeSL = textLength + 1 - iTypeSL;	// iTypeSL is iTypeL
			}
		}
		
		{
			// allocate workspace 1 - free memory before calling recursive function
			// total(work space 1 & 2) maximum size is less than 3*input_length word (eg. 12n)
			// n: input size, @: alphabet size
			// size = max(@ + m, n - beta + m)
			final int workSpace1Size = textLength > (iTypeSL * 2) ? textLength : (iTypeSL * 2);
			final int[] workSpace1 = new int[workSpace1Size + iTypeSL];
			
			// name blocks (partitions) of working memory according to their usage
			final int[] dist = workSpace1;		// record S(L) distances
			final int[] sub = workSpace1;		// record S(L) list
			final int[] subLink = workSpace1;	// ????
			final int subLinkOffset = iTypeSL;
			
			// Record S(L)-distances into Dist (KA p6 1.)
			// length of mbucket_count is unknown: reset maximum length
			final int[] mbucket_count = suffixArray;	// record sizes of buckets of m-list
			//Arrays.fill(mbucket_count, 0, Math.min(textLength * INT_SIZE, suffixArray.length), 0);
			Arrays.fill(mbucket_count, 0, textLength * INT_SIZE, 0);
			
			int beta = 0;
			int m = 0;
			dist[0] = 0;
			for(i = 0; sign(i, text, signBit) != (signBitSorL); i++){
				dist[i] = 0;
				beta++;
			}
			dist[i] = 0;
			beta++;
			i++;
			int currentDist = 0;
			
			for(; i<textLength; i++){
				currentDist++;
				dist[i] = currentDist;
				if (currentDist > m) m = currentDist;	// m = max(Dist[i])
				mbucket_count[currentDist-1]++;
				if (sign(i, text, signBit) == (signBitSorL))	// i-1 is type S(L)
				{
					currentDist = 0;
				}
			}
			
			// allocate workspace 2 - free memory before calling recursive function
			final int workSpace2Size = (alphabetSize + m) > (textLength - beta + m) ? (alphabetSize + m) : (textLength - beta + m);
			// allocate working memory
			final int[] workSpace2 = new int[Math.max(workSpace2Size + m, alphabetSize * INT_SIZE)];
			// name blocks of working memory according to their usage
			final int[] mbucket = workSpace2;	// record start positions of buckets of the m-list, size m
			final int[] bucket1 = workSpace2;	// content of the buckets of the m-list, overlap with mlist, size alpha
			final int bucket1Offset = m;
			final int[] mlist = workSpace2;		// content of the buckets of the m-list, size iT - beta
			final int mlistOffset = m;
			
			// build mbucket from mbucket_count
			j = 0; 
			for (i = 0; i < m; i++)
			{
				mbucket[i] = j;
				j += mbucket_count[i];
			}
			
			// bucket all suffixes of T according to first character (KA p4 1.) into array 'a'
			final int[] a = suffixArray;
			//Arrays.fill(bucket1, bucket1Offset, Math.min(alphabetSize * INT_SIZE, bucket1.length), 0);
			Arrays.fill(bucket1, bucket1Offset, Math.min(bucket1Offset + alphabetSize * INT_SIZE, bucket1.length), 0);
			
			// Count occurance of each alphabet
			for (i = 0; i < textLength; i++) bucket1[value(i, text, signBit) + bucket1Offset]++;
			
			// Convert Bucket: alphabet count > index of starting point of each bucket
			// Correctly assign index for alphabet with count 0 (points to the starting index of next nonzero count alphabet)
			j = 0; 
			
			int temp;
			for (i = 0; i < alphabetSize; i++)
			{
				temp = bucket1[i + bucket1Offset];
				bucket1[i + bucket1Offset] = j;
				j += temp;
			}
			
			// bucketing: fill A
			for (i = 0; i < textLength; i++)
			{
				a[bucket1[value(i, text, signBit) + bucket1Offset]] = i;
				bucket1[value(i, text, signBit) + bucket1Offset]++;
			}
			
			// create m-list (KA p6 2.)
			// go through A left-right	
			for (i = 0; i < textLength; i++)	
			{
				final int dst = dist[a[i]];
				if (dst > 0)
				{
					mlist[mbucket[dst-1]++ + mlistOffset] = a[i];

				}
			}
			assert(m == 0 || mbucket[m-1] == textLength - beta);
			
			// sort the type S(L) substrings (KA p6 3.)
			// init sub and sublink
			// initialize, bucketing base on the first character
			// go through A left-right
			// the other way round in case of typeL
			if (signBit == signBitSorL)	// typeS
			{
				j = 0;
				for (i = 0; i < textLength; i++) 
				{
					if (sign(a[i], text, signBit) == signBitSorL) 
					{
						sub[j] = a[i];
						assert(value(a[i], text, signBit) >= Integer.MIN_VALUE && value(a[i], text, signBit) <= Integer.MAX_VALUE);
						y = (char) value(a[i], text, signBit);
						if (j == 0)	// first in sub
						{
							k = 0;
							x = y;
							subLink[0 + subLinkOffset] = 0;
						}
						else
						{
							if (x != y)		// a new bucket front
							{
								k = j;
								x = y;
							}
							subLink[j + subLinkOffset] = k;
						}
						j++;
					}
				}
				assert(j == iTypeSL);
			}
			else						// typeL
			{
				j = iTypeSL-1;
				for (i = textLength-1; i != -1; i--) 
				{
					if (sign(a[i], text, signBit) == signBitSorL) 
					{
						sub[j] = a[i];
						assert(value(a[i], text, signBit) <= Integer.MAX_VALUE && value(a[i], text, signBit) >= Integer.MIN_VALUE);
						y = (char) value(a[i], text, signBit);
						if (j == iTypeSL-1)	// first in sub
						{
							k = iTypeSL-1;
							x = y;
							subLink[iTypeSL-1 + subLinkOffset] = iTypeSL-1;
						}
						else
						{
							if (x != y)		// a new bucket front
							{
								k = j;
								x = y;
							}
							subLink[j + subLinkOffset] = k;
						}
						j--;
					}
				}
				assert(j == -1);
			}
			
			final int[] subR = suffixArray;  // record reverse S(L) list (position of element x in S(L)-list)
			// init subR, which is inverse of sub
			// -1 indicates no corresponding suffix in sub
			//Arrays.fill(subR, 0, Math.min(textLength * INT_SIZE, suffixArray.length), -1);
			Arrays.fill(subR, 0, textLength * INT_SIZE, -1);
			
			for (i = 0; i < iTypeSL; i++)
			{
				subR[sub[i]] = i;
			}
			
			// go through m-list left-right, repeated bucketing on sub
			// go througn m-list right-left in case of typeL to put proper-prefix-substrings to front of bucket
			if (signBit == signBitSorL)	// typeS 
			{
				i = 0;
				for (j = 0; j < m; j++)	// go through mbucket
				{
					final int dst = j + 1;
					for (; i < mbucket[j]; i++)
					{
						final int suffix = mlist[i + mlistOffset] - dst;		// suffix to be moved to its bucket front
						final int current = subR[suffix];	// current position in sub
						int link = subLink[current + subLinkOffset];			

						if (link >= current)	// if current is the front, also this should be the first time this bucket is accessed in this round
						{
							//link = current + 1;	// set correct next position
							subLink[current + subLinkOffset] = current + 1;
						}
						else	// current is not its bucket front
						{
							int bucketCurrent = subLink[link + subLinkOffset];

							// reset bucketcurrent if needed this is the first access in round (also means first access in this round)
							if (bucketCurrent > current || bucketCurrent == link)	// this is the first access of this bucket in this round
							{
								// swap is always necessary because bucketcurrent is bucketfirst and current is not front
								{
									// swap sub value
									int tmp = sub[current];
									sub[current] = sub[link];
									sub[link] = tmp;

									// update subR
									subR[sub[current]] = current;
									subR[sub[link]] = link;
								}

								//bucketCurrent = link + 1;	// bucketcurrent is set correctly
								subLink[link + subLinkOffset] = link + 1;
							}
							else	// this is not first access in this round, thus bucketcurrent is not at bucketfirst
							{
								// swap if necessary
								if (bucketCurrent != current)	
								{
									// swap sub value
									int tmp = sub[current];
									sub[current] = sub[bucketCurrent];
									sub[bucketCurrent] = tmp;

									// update subR
									subR[sub[current]] = current;
									subR[sub[bucketCurrent]] = bucketCurrent;
								}

								// check if bucketcurrent makes new bucketfront for the next round
								// and update sublink
								if (value(mlist[i + mlistOffset], text, signBit) != value(sub[bucketCurrent-1]+dst, text, signBit))	// is new front
								{
									subLink[bucketCurrent + subLinkOffset] = bucketCurrent;
								}
								else	// not new front
								{
									// special case of second place
									if (link == bucketCurrent - 1)
									{
										subLink[bucketCurrent + subLinkOffset] = link;
									}
									else	// general case
									{
										subLink[bucketCurrent + subLinkOffset] = subLink[bucketCurrent-1 + subLinkOffset];
									}
								}

								// increment bucketcurrent
								//bucketCurrent++;
								subLink[link + subLinkOffset]++;
							}
						}
					}
				}
			}
			else	// typeL
			{
				for (j = 0; j < m; j++)	// go through mbucket in reverse order
				{
					final int dst = j + 1;
					final int mbucketfront = (j > 0) ? mbucket[j-1] : 0;
					for (i = mbucket[j] - 1; i >= mbucketfront; i--)
					{
						final int suffix = mlist[i + mlistOffset] - dst;		// suffix to be moved to its bucket front
						final int current = subR[suffix];	// current position in sub
						int link = subLink[current + subLinkOffset];

						if (link <= current)	// if current is the front, also this should be the first time this bucket is accessed in this round
						{
							//! link = current - 1;	// set correct next position
							subLink[current + subLinkOffset] = current - 1;
						}
						else	// current is not its bucket front
						{
							int bucketCurrent = subLink[link + subLinkOffset];

							// reset bucketcurrent if needed this is the first access in round (also means first access in this round)
							if (bucketCurrent < current || bucketCurrent == link)	// this is the first access of this bucket in this round
							{
								// swap is always necessary because bucketcurrent is bucketfirst and current is not front
								{
									// swap sub value
									int tmp = sub[current];
									sub[current] = sub[link];
									sub[link] = tmp;

									// update subR
									subR[sub[current]] = current;
									subR[sub[link]] = link;
								}

								//! bucketCurrent = link - 1;	// bucketcurrent is set correctly
								subLink[link + subLinkOffset] = link - 1;
							}
							else	// this is not first access in this round, thus bucketcurrent is not at bucketfirst
							{
								// swap if necessary
								if (bucketCurrent != current)	
								{
									// swap sub value
									int tmp = sub[current];
									sub[current] = sub[bucketCurrent];
									sub[bucketCurrent] = tmp;

									// update subR
									subR[sub[current]] = current;
									subR[sub[bucketCurrent]] = bucketCurrent;
								}

								// check if bucketcurrent makes new bucketfront for the next round
								// and update sublink
								if (value(mlist[i + mlistOffset], text, signBit) != value(sub[bucketCurrent+1]+dst, text, signBit))	// is new front
								{
									subLink[bucketCurrent + subLinkOffset] = bucketCurrent;
								}
								else	// not new front
								{
									// special case of second place
									if (link == bucketCurrent + 1)
									{
										subLink[bucketCurrent + subLinkOffset] = link;
									}
									else	// general case
									{
										subLink[bucketCurrent + subLinkOffset] = subLink[bucketCurrent+1 + subLinkOffset];
									}
								}

								// increment bucketcurrent
								//! bucketCurrent--;
								subLink[link + subLinkOffset]--;
							}
						}
					}
				}
			}
			// ^-!!!!
			
			// Convert subR to suffix-encoding table
			// subR: suffix to encoding 1 ~ less than (iTypeSL), because some substrings should be assigned identical number
			// (0 is reserved for TERMIN)
			// param iTdashMaxAlphabet
			iTdashMaxAlphabet = 0;
			if (signBit == signBitSorL)	// typeS
			{
				for (i = 0; i < iTypeSL; i++)
				{
					if (subLink[i + subLinkOffset] >= i)	// this substring is a bucket front
					{
						iTdashMaxAlphabet++;
					}
					subR[sub[i]] = iTdashMaxAlphabet;
				}
			}
			else						// typeL
			{

				iTdashMaxAlphabet = 1;
				for (i = 0; i < iTypeSL; i++)
				{
					subR[sub[i]] = iTdashMaxAlphabet;
					if (subLink[i + subLinkOffset] <= i)	// this substring is a bucket front
					{
						iTdashMaxAlphabet++;
					}
				}
				iTdashMaxAlphabet--;
			}
			
			// subR is converted to sufenc
			sufenc = suffixArray;
			
			// free work space 1 & 2 before recursive call by leaving this block
			
			if(PRINT_MEMORY_USAGE){
				NumberFormat format = NumberFormat.getInstance();
				System.gc();
				System.gc();
				System.gc();
				Runtime runtime = Runtime.getRuntime();
				System.out.println("Memory after workspace1&2: " + format.format(runtime.totalMemory() - runtime.freeMemory()));
			}
		}
		
		// generate Tdash (KA p7)
		// determine length and alphabet size of Tdash
		// cf. Tdash length = iTypeSL+1, Tdash alphabet = 0 ~ iTdashMaxAlphabet
		// determine alphabet size
		if (iTdashMaxAlphabet <= Character.MAX_VALUE){
			suffixArray = recursivesort(sufenc, textLength, iTypeSL, iTdashMaxAlphabet);
			if(DEBUG){
				System.out.println("Rec.sort using char alphabet");
			}
		} 
		else if (iTdashMaxAlphabet <= Short.MAX_VALUE){
			suffixArray = recursivesort(sufenc, textLength, iTypeSL, iTdashMaxAlphabet);
			if(DEBUG) System.out.println("Rec.sort using short alphabet");
		}
		else{
			suffixArray = recursivesort(sufenc, textLength, iTypeSL, iTdashMaxAlphabet);
			if(DEBUG) System.out.println("Rec.sort using int alphabet");
		}
		
		final int[] saDash = suffixArray;
		
		{
			// allocate Work Space 3 for post-recursion processing - freed at the end of function - size = max(typeSL * 2, alpha+typeSL)
			final int workSpace3Size = (iTypeSL * 2) > (iTypeSL + alphabetSize) ? (iTypeSL * 2) : (iTypeSL + alphabetSize);
			final int[] workSpace3 = new int[Math.max(workSpace3Size + iTypeSL, alphabetSize * INT_SIZE)];
			final int[] b = workSpace3;
			final int[] saDash2SA = workSpace3;
			final int saDash2SAOffset = iTypeSL;
			final int[] bucket2 = workSpace3;
			final int bucket2Offset = iTypeSL;
			
			// prepare B, B is array of type S(L) suffixes sorted in lexicographic order (KA p4)
			j = 0;
			for (i = 0; i < textLength; i++) 
			{
				if (sign(i, text, signBit) == signBitSorL) 
				{
					saDash2SA[saDash2SAOffset + j++] = i;
				}
			}
			// do mapping, create B
			for (i = 0; i < iTypeSL; i++)		
			{
				b[i] = saDash2SA[saDash[i+1] + saDash2SAOffset];	// SAdash[0] contains termination '$' which we added
			}
			
			// now we have B, whis is S(L) suffixes sorted
			if (signBit == signBitSorL)	// for typeS
			{
				// Move type S suffixes to end of its bucket in A (KA p4 2.)
				// make Bucket array initialized to end of bucket
				// Count occurance of each alphabet
				//Arrays.fill(bucket2, bucket2Offset, Math.min(alphabetSize * INT_SIZE, bucket2.length), 0);
				Arrays.fill(bucket2, bucket2Offset, Math.min(bucket2Offset + alphabetSize * INT_SIZE, bucket2.length), 0);
				
				for (i = 0; i < textLength; i++) bucket2[value(i, text, signBit) + bucket2Offset]++;
				// Convert Bucket: alphabet count -> index of end point of each bucket
				// Correctly assign index for alphabet with count 0 (points to the starting index of next nonzero count alphabet)
				j = 0; 
				for (i = 0; i < alphabetSize; i++)
				{
					int tmp = bucket2[i + bucket2Offset];
					bucket2[i + bucket2Offset] = j + tmp - 1;
					assert(bucket2[i + bucket2Offset] < textLength);
					j += tmp;
				}

				// scan B from right to left, sort type S suffixes inside its buckets 
				for (i = iTypeSL - 1; i != -1; i--)
				{
					assert(bucket2[value(b[i], text, signBit) + bucket2Offset] >= 0 && b[i] >= 0);
					suffixArray[bucket2[value(b[i], text, signBit) + bucket2Offset]--] = b[i];
				}

				// Sort type L suffixes (KA p4 3.)
				// initialize BucketCounts to point the front of each bucket
				// make Bucket array initialized to front of bucket
				// Count occurance of each alphabet
				//Arrays.fill(bucket2, bucket2Offset, Math.min(alphabetSize * INT_SIZE, bucket2.length), 0);
				Arrays.fill(bucket2, bucket2Offset, Math.min(bucket2Offset + alphabetSize * INT_SIZE, bucket2.length), 0);
				
				for (i = 0; i < textLength; i++) bucket2[value(i, text, signBit) + bucket2Offset]++;
				// Convert Bucket: alphabet count -> index of starting point of each bucket
				// Correctly assign index for alphabet with count 0 (points to the starting index of next nonzero count alphabet)
				j = 0; 
				for (i = 0; i < alphabetSize; i++)
				{
					int tmp = bucket2[i + bucket2Offset];
					bucket2[i + bucket2Offset] = j;
					j += tmp;
				}

				// scan SA from left to right
				for (i = 0; i < textLength; i++)
				{
					final int SuffixToPut = suffixArray[i] - 1;
					if (SuffixToPut == -1) continue;	// skip for the first suffix

					if (sign(SuffixToPut, text, signBit) != signBitSorL)	// if it is type L suffix
					{
						assert(bucket2[value(SuffixToPut, text, signBit) + bucket2Offset] >= 0 && SuffixToPut >= 0);
						suffixArray[bucket2[value(SuffixToPut, text, signBit) + bucket2Offset]++] = SuffixToPut;
					}
				}
			}
			else	// Type L
			{
				// Move type L suffixes to FRONT of its bucket in A 
				// make Bucket array initialized to front of bucket
				// Count occurance of each alphabet
				//Arrays.fill(bucket2, bucket2Offset, Math.min(alphabetSize * INT_SIZE, bucket2.length), 0);
				Arrays.fill(bucket2, bucket2Offset, Math.min(bucket2Offset + alphabetSize * INT_SIZE, bucket2.length), 0);
				
				for (i = 0; i < textLength; i++) bucket2[value(i, text, signBit) + bucket2Offset]++;
				// Convert Bucket: alphabet count -> index of starting point of each bucket
				// Correctly assign index for alphabet with count 0 (points to the starting index of next nonzero count alphabet)
				j = 0; 
				for (i = 0; i < alphabetSize; i++)
				{
					int temp = bucket2[i + bucket2Offset];
					bucket2[i + bucket2Offset] = j;
					j += temp;
				}

				// scan B from LEFT to RIGHT, sort type L suffixes inside its buckets 
				for (i = 0; i < iTypeSL; i++)
				{
					assert(bucket2[value(b[i], text, signBit) + bucket2Offset] >= 0 && b[i] >= 0);
					suffixArray[bucket2[value(b[i], text, signBit) + bucket2Offset]++] = b[i];
				}

				// Sort type S suffixes
				// initialize BucketCounts to point the END of each bucket
				// make Bucket array initialized to end of bucket
				// Count occurance of each alphabet
				//Arrays.fill(bucket2, bucket2Offset, Math.min(alphabetSize * INT_SIZE, bucket2.length), 0);
				Arrays.fill(bucket2, bucket2Offset, Math.min(bucket2Offset + alphabetSize * INT_SIZE, bucket2.length), 0);
				
				for (i = 0; i < textLength; i++) bucket2[value(i, text, signBit) + bucket2Offset]++;
				// Convert Bucket: alphabet count -> index of end point of each bucket
				// Correctly assign index for alphabet with count 0 (points to the starting index of next nonzero count alphabet)
				j = 0; 
				for (i = 0; i < alphabetSize; i++)
				{
					int tmp = bucket2[i + bucket2Offset];
					bucket2[i + bucket2Offset] = j + tmp - 1;
					assert(bucket2[i + bucket2Offset] < textLength);
					j += tmp;
				}

				// scan A from RIGHT to LEFT
				for (i = textLength-1; i != -1; i--)
				{
					final int SuffixToPut = suffixArray[i] - 1;
					if (SuffixToPut == -1) continue;	// skip for the first suffix

					if (sign(SuffixToPut, text, signBit) != signBitSorL)	// if it is type S suffix
					{
						assert(bucket2[value(SuffixToPut, text, signBit) + bucket2Offset] >= 0 && SuffixToPut >= 0);
						suffixArray[bucket2[value(SuffixToPut, text, signBit) + bucket2Offset]--] = SuffixToPut;
					}
				}
			}
			
			// free work space 3
			if(PRINT_MEMORY_USAGE){
				NumberFormat format = NumberFormat.getInstance();
				System.gc();
				System.gc();
				System.gc();
				Runtime runtime = Runtime.getRuntime();
				System.out.println("Memory after workspace3: " + format.format(runtime.totalMemory() - runtime.freeMemory()));
			}
		}
		return suffixArray;
	}
	
	/**
	 * Shift suffixes by one position for later insertion of termination character before all other characters.
	 * 
	 * @param suffixes the list of suffixes to be shifted
	 * @return the shifted list of suffixes
	 */
	private int[] shiftSuffixes(char[] text){
		if(KoAluruSuffixSort.ADD_TERMIN){
			int[] result = new int[text.length + 1];
			for(int i=0; i<text.length; i++){
				assert((text[i] + 1) >= Integer.MIN_VALUE && (text[i] + 1) <= Integer.MAX_VALUE);
				result[i] = (char) (text[i] + 1);
			}
			result[text.length] = KoAluruSuffixSort.TERMIN;
			
			return result;
		}
		else{
			int[] result = new int[text.length];
			for(int i=0; i<text.length; i++){
				assert((text[i] + 1) >= Integer.MIN_VALUE && (text[i] + 1) <= Integer.MAX_VALUE);
				result[i] = text[i] + 1;
			}
			return result;
		}
	}
	
	private final int value(int x, int[] text, int signBit){
		assert(x >= 0 && x < text.length);
		return text[x] & ~signBit;
	}
	
	private final int sign(int x, int[] text, int signBit){
		assert(x >= 0 && x < text.length);
		return text[x] & signBit;
	}
	
	public static void main(String[] args) throws IOException{
		
//		char[] text;
//		{
//			SlimFastaParser parser = new SlimFastaParser(new File(args[0]));
//			SimpleContig contig;
//			StringBuffer buffy = new StringBuffer();
//			while(parser.hasNextContig()){
//				contig = parser.parseNextContig();
//				buffy.append(contig.getSequence());
//				if(parser.hasNextContig()) buffy.append('N');
//			}
//			
//			text = buffy.toString().toCharArray();
//		}
//		NumberFormat format = NumberFormat.getInstance();
//		System.out.println("Textlength: " + format.format(text.length));
		
		//char[] text = new char[]{'d','c','b','a'};
		
		//char[] text = new char[]{'a','b','c','d'};
		//char[] text = new char[]{'a','d','c','b'};
		//char[] text = new char[]{'a','c','d','b'};
		
		//char[] text = "acaaacatat".toCharArray();
		//char[] text = "ACAAACATAT".toCharArray();
		char[] text = "MISSISSIPPI".toCharArray();
		
		SimpleTimer timer = new SimpleTimer();
		KoAluruSuffixSort sort = new KoAluruSuffixSort();
		int[] result = sort.getSuffixArray(text);
		for(int i=0; i < result.length; i++){
			if(result[i]<text.length) System.out.println(i + ": " + result[i] + " " + new String(text).substring(result[i]) + "$");
			else System.out.println(i + ": " + result[i] + " $");
		}
		System.out.println("Sorted in: " + timer.getTimeString());
		
		BufferedWriter writer = new BufferedWriter(new FileWriter("ka-java.out"));
		for(int i=0; i<result.length; i++){
			writer.write(result[i] + "\n");
		}
		writer.close();
//		System.gc();
//		System.gc();
//		System.gc();
//		Runtime runtime = Runtime.getRuntime();
//		System.out.println("Memory: " + format.format(runtime.totalMemory() - runtime.freeMemory()));
	}
}
