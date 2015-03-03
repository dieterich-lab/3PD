package primerDesign.Test;

//	 Copyright(c) 2006, Shinichi Morishita. All Rights Reserved.

public class SuffArr {

	// Larsson-Sadakane Algorithm
	public static void suffixArray_LS(int[] S, int[] SA){
		// Assume that all elements in S are positive except that the last is zero.
		int n = S.length;
		int[] ISA = new int[n];
		// Basic step.
		radixSort_LS(S, SA, ISA);
		//Iterative step.
		boolean iteration; int h=1;
		do{		
			iteration = false;
			for(int i=0; i<n; i = ISA[SA[i]]+1)
				if(i < ISA[SA[i]]){ // If some blocks are divided.
					ternary_split_quicksort_LS(SA, i, ISA[SA[i]], ISA, h);
					// Repeat the while-loop. 
					iteration = true;
				}
			h = h*2;
		}while(iteration);
	}
	public static void radixSort_LS(int[] S, int[] SA, int[] ISA){	
		// Treat each element as one digit, execute radix sort, and update ISA.
		int n = S.length;
		int m=0;
		for(int i=0; i<n; i++) 	// Calculate the maximum.
			if(m < S[i]) m = S[i];
		m++; // Total number of elements
		int[] count = new int[m];   
		for(int i=0; i<m; i++){ count[i]=0; } 
		for(int i = 0; i < n; i++){ count[S[i]]++; } 
		for(int i=1; i<m; i++){ count[i] += count[i-1]; } 
		for(int i=n-1; 0<=i; i--){ ISA[i] = count[S[i]]-1; }  // update ISA
		for(int i = n-1; 0 <= i; i--){ SA[count[S[i]]-1] = i; count[S[i]]--; }
	}	
	public static void ternary_split_quicksort_LS
					(int[] SA, int l, int r, int[] ISA, int h) {
		if(l > r) return;	 // Exit when the empty range is input.
		if(l == r){ // Set ISA[SA[l]] to l when the input is singleton.
			ISA[SA[l]] = l; return; }
		// Choose a pivot from ISA[SA[]+h].		
		int v = ISA[SA[l+(int)Math.floor(Math.random()*(r-l+1))]+h]; 
		int i = l; int mi = l; int j = r; int mj = r; int tmp;
		for (;;) {  // Compare values according to ISA[SA[]+h]
			for(; i <= j && ISA[SA[i]+h] <= v; i++) 
				if(ISA[SA[i]+h] == v){  
					tmp = SA[i]; SA[i] = SA[mi]; SA[mi] = tmp; mi++; }
			for(; i <= j && v <= ISA[SA[j]+h]; j--)
				if(ISA[SA[j]+h] == v){ 
					tmp = SA[j]; SA[j] = SA[mj]; SA[mj] = tmp; mj--; }			
			if (i > j)  break; 
			tmp = SA[i];  SA[i] = SA[j];  SA[j] = tmp; i++; j--;
		}
		for(mi--, i--; l <= mi; mi--, i--){
			tmp = SA[i];  SA[i] = SA[mi];  SA[mi] = tmp; }
		for(mj++, j++; mj <= r; mj++, j++){
			tmp = SA[j];  SA[j] = SA[mj];  SA[mj] = tmp; }
		ternary_split_quicksort_LS(SA, l, i, ISA, h);
		for(int k = i+1; k < j; k++) ISA[SA[k]] = j-1;  // update ISA
		ternary_split_quicksort_LS(SA, j, r, ISA, h);
	}
	public static int digit(int num, int j, int m){
		int quotient = num;
		int remainder = 0;
		for(int i=0; i<j; i++){
			remainder = quotient%m;
			quotient = quotient/m;
		}
		return remainder;
	}		
	
	// Karkkainen-Sanders Algorithm
	public static void suffixArray_KS(int[] S, int n, int[] SA, int numRanks){
		if( !(1 < n && 0 < numRanks) ) return;
		// Elements in S must be positive integers except the last three zeros.
		// n = S.length - 3. 0 < n. numRanks = number of symbols in S.
		int n0=(n+2)/3; int n1=(n+1)/3; int n2=n/3;  
		int n02 = n0+n2;
		// If n0=n1+1, an extra element of S12 at position n1 denotes
		// the triplet "000" of S at position n+1. 
		int[] S12 = new int[n02+3]; S12[n02]=S12[n02+1]=S12[n02+2]=0;
		int[] rank12 = new int[n02];				
		for(int j=0; j < n02; j++) S12[j] = f(j, n0);
		int[] count = new int[numRanks]; 	
		// Sort triplets by radix sort.
		for(int d=2; 0<=d; d--){
			for(int i=0; i < numRanks; i++){ count[i]=0; } 
			for(int j=0; j < n02; j++) count[S[S12[j]+d]]++;
			for(int i=1; i < numRanks; i++) count[i] += count[i-1];	
			for(int j=n02-1; 0 <= j; j--) rank12[--count[S[S12[j]+d]]] = S12[j];
			for(int j=0; j < n02; j++) S12[j] = rank12[j];
		}
		// Represent triplets by their ranks.
		int rank = 1;
		for(int j=0; j < n02; j++){ 
			//	Increment rank if neighboring triplets are not equal.
			if(0<j){
				int x = rank12[j-1]; int y = rank12[j];
				if(compare(S[x],S[x+1],S[x+2],S[y],S[y+1],S[y+2])!=0) rank++;
			}
			if(rank12[j]%3==1) S12[rank12[j]/3] = rank;		// i = 1 mod 3
			else 			   S12[rank12[j]/3+n0] = rank;	// i = 2 mod 3
		}	
		// Recursive call to generate the suffix array and its inverse for S12.		
		int[] SA12 = rank12; 	// Reuse rank12 as it is not needed at this point.			
		if(rank == n02) // No duplicate ranks. Sort suffixes by their first integers. 
			for(int j=0; j<n02; j++) SA12[S12[j]-1]=j; // Adjustment to 0-origin. 
		else suffixArray_KS(S12, n02, SA12, rank+1);	
		int[] ISA12 = S12;    	// Reuse S12 array as it is not needed.
		for(int j=0; j<n02; j++) ISA12[SA12[j]] = j+1;  // ISA12[j] = ISA12[j]+1
		
		// Sort elements in S0.		
		int[] S0  = new int[n0];
		// ISA12[j]+1 for j=0,...,n0-1 have been sorted implicitly in SA12.
		for(int j=0, k=0; j < n02; j++)  
			if(SA12[j]<n0) S0[k++] = 3*SA12[j]; // i=f(j) and i=1 mod 3. 

		// Sort the first elements by radix sort.
		int[] SA0 = new int[n0];  // Sorted indexes of S0
		for(int i=0; i < numRanks; i++) count[i]=0; 
		for(int k=0; k < n0; k++)  count[ S[S0[k]] ]++;
		for(int i=1; i < numRanks; i++) count[i] += count[i-1];	
		for(int k=n0-1; 0 <= k; k--) SA0[--count[S[S0[k]]]] = S0[k];		
		// Merge two sorted lists.
		int j, k, l;
		for(j=n0-n1, k=0, l=0; j<n02 && k<n0; l++){
			int comparison;
			int x = f(SA12[j],n0); 
			int y = SA0[k];			
			if(SA12[j]<n0){ // x = 1 mod 3.
				comparison = compare(S[x], ISA12[SA12[j]+n0], S[y], ISA12[y/3]);
				// If y is the last position of S, y/3=n0-1 and n0=n1+1.
				// Since S12[y/3] is arranged to denote "000", ISA[y/3] must be 0. 
			}else{ // x = 2 mod 3.
				comparison = compare(S[x], S[x+1], ISA12[SA12[j]-n0+1], S[y], S[y+1], ISA12[y/3+n0]);	
				// If x is the last position of S, S[x+1] must be 0, 
				// indicating that the former triplet must be lower than the latter.			
			}
			// S[x]... < S[y]... if and only if comparison == -1
			if(comparison == -1){ SA[l]=x; j++; }
			else{ SA[l]=y; k++; }
		}	
		for(; j<n02; l++, j++) SA[l]=f(SA12[j],n0);
		for(; k<n0;  l++, k++) SA[l]=SA0[k];
	}
	public static int f(int j, int n0){
		if(j < n0) return 3*j+1; else return 3*(j-n0)+2;
	}
	public static int compare(int i0, int i1, int j0, int j1){
		// return -1, 0, 1 if i < j, i == j, and i > j, respectively.
		if(i0 < j0) return -1;
		else if(i0 == j0)
			if(i1 <  j1) return -1; else if(i1 == j1) return 0; else return 1;
		else // i0 > j0
			return 1;
	}	
	public static int compare(int i0, int i1, int i2, int j0, int j1, int j2){
		// return -1, 0, 1 if i < j, i == j, and i > j, respectively.
		if(i0 <  j0)      return -1;
		else if(i0 == j0) return compare(i1, i2, j1, j2);
		else              return 1;	
	}

	// Simple but naive binary search algorithm for the query.
	public static int searchLeftmost(int[] S, int[] SA, int[] query){
		// Exit if the query, which is greater than "$", is outside the scope.
		if(compare(S,SA[SA.length-1], query) == -1) return -1;
		int left, right;
		for(left = 0, right = SA.length-1; left+1 < right; ){
			int middle = (left + right) / 2;
			if(compare(S, SA[middle], query) == -1) left = middle;
			else right = middle;
		}
		if( compare(S, SA[right], query) == 0 ) return right;
		else return -1; // No occurrences of the query.
	}
	
	public static int compare(int S[], int start, int[] query){
		// Return 0 if the query occurs as a prefix of the suffix S[start]. 
		// Otherwise, return -1 and 1 if the suffix from "start" is lexicographically  
		// lower or higher than query.
		for(int i=0; start + i < S.length && i < query.length; i++)
			if(S[start+i] < query[i]) return -1;
			else if(S[start+i] > query[i]) return 1;
		return 0;
	}	
	
	// Advanced binary search algorithm for the query with lcp information.
	public static int searchLeftmost(int[] S, int[] SA, int[] query, 
            int[] LCP, int[] LCP_AUX){
		// Compute the lcp values between the query and both end suffixes.
		int left  = 0;            int lcp_left_query  = 0; 
		int right = SA.length-1;  int lcp_right_query = lcp0(0, query, S, SA[right]);
		// Exit if the query, which is greater than "$", is outside the scope.  
		if(!less_eq(lcp_right_query, query, S, SA[right])) return -1;
		// Binary search.
		for(int middle=(left+right)/2; left+1 < right; middle=(left+right)/2){ 
			if(lcp1(left, middle, LCP, LCP_AUX) >= lcp1(middle, right, LCP, LCP_AUX)){
				if(lcp_left_query < lcp1(left, middle, LCP, LCP_AUX))
					left = middle; // lcp_left_query remains unchanged.
				else if(lcp_left_query > lcp1(left, middle, LCP, LCP_AUX)){
					right = middle; lcp_right_query = lcp1(left, middle, LCP, LCP_AUX);
				}else{ // Set the lcp of the query and the middle suffix to i.
					int i = lcp0( lcp1(left, middle, LCP, LCP_AUX), query, S, SA[middle]);
					if(less_eq(i, query, S, SA[middle])){  	
						right = middle; lcp_right_query = i;   
						}else{ left  = middle; lcp_left_query  = i;}				
					}
			}else{
				if(lcp_right_query < lcp1(middle, right, LCP, LCP_AUX)){
					right = middle; // lcp_right_query remains unchanged.
				}else if(lcp_right_query > lcp1(middle, right, LCP, LCP_AUX)){
					left = middle; lcp_left_query = lcp1(middle, right, LCP, LCP_AUX);
				}else{ // Set the lcp of the query and the middle suffix to i.
					int i = lcp0( lcp1(middle, right, LCP, LCP_AUX), query, S, SA[middle]);
					if(less_eq(i, query, S, SA[middle])){  	
						right = middle; lcp_right_query = i;   
					}else{ left  = middle; lcp_left_query  = i;}
				}
			}
		}
		if(lcp_right_query == query.length) return right; else return -1;
	}
	public static boolean less_eq(int lcp, int[] query, int[] S, int start){
		// Return true if the query is lower than or equal to the suffix.
		if(lcp == query.length) return true;
		if(start+lcp < S.length && query[lcp] < S[start+lcp]) return true;
		return false;
	}
	public static int lcp0(int offset, int[] query, int[] S, int start){ 
		// Compute lcp.
		int i = offset;
		for(; i < query.length && start+i < S.length && 
				query[i] == S[start+i];  i++);
		return i;
	}
	public static int lcp1(int i, int j, int[] LCP, int[] LCP_AUX){ 
		// Lookup lcp tables.
		if(i+1 == j) return LCP[j]; else return LCP_AUX[(i+j)/2]; 
	}	
	
	// Build longest common prefixes.
	public static void buildLcp(int[] S, int[] SA, int[] LCP, int[] LCP_AUX){
		// Assume that the last symbol of S is "$", and SA is the suffix array of S.
		// int[] LCP = new int[S.length]; int[] LCP_AUX = new int[S.length-1];
		int n = SA.length;
		int[] ISA = new int[n];
		for(int i=0; i<n; i++) ISA[SA[i]]=i; // Build the inverse suffix array.
		for(int lcp=0, i=0; i<n; i++) {
			int k = ISA[i];		    // Suffixes are processed from longer ones. 
			if(k==0) LCP[k] = -1;   // The last suffix "$".
			else {				
				int j = SA[k-1];	// Note that j=SA[k-1] and i=SA[k] since k=ISA[i].
				while(j+lcp<=n && i+lcp<=n && S[j+lcp]==S[i+lcp]) lcp++;
				LCP[k] = lcp;
			}
			if(lcp > 0) lcp--; 
		}
		// Compute lcp values for all possible intervals and put them into LCP_AUX
		// according to LCP_AUX[(l+r)/2] = lcp(l,r) if l+1 < r.
		buildLcp1(0, SA.length-1, LCP, LCP_AUX);  
	}
	public static int buildLcp1(int l, int r, int[] LCP, int[] LCP_AUX){
		if(l+1 == r) return LCP[r];
		else{ int v = Math.min(buildLcp1(l, (l+r)/2, LCP, LCP_AUX), 
					           buildLcp1((l+r)/2, r, LCP, LCP_AUX));
			LCP_AUX[(l+r)/2] = v;
			return v;
		}
	}
	
	// Calculate occurrence frequencies of k-mers
	public static int[] occurrence_frequency(int[] SA, int[] LCP, int k){
		int[] freq = new int[LCP.length];
		int runLength = 1;
		for(int i = 0; i < LCP.length; i++){
			// Checks if k > lcp(i,i+1) or not.
			if((i < LCP.length-1 && k > LCP[i+1]) || i == LCP.length-1){ 
				// The suffix at SA[i+1] is new, or the end of LCP is hit.
				// Store the occurrence frequency of the current k-mer. 
				for(int j = i+1-runLength; j <= i; j++) 
					freq[SA[j]] = runLength;
				// Initialize for the new k-mer.
				runLength = 1;
			}else // Continue to search the current k-mer.
				runLength++;
		}
		return freq;
	}
	
	// Calculate the longest common factor
	public static int lcf(int[] ISA, int[] LCP, int left, int right){ 
		// The query ranges from left to right in S.
		if(!(0 <= left && left <= right && right < ISA.length)) return 0;
		int answer = 0;
		for(int k = left; k <= right; k++){
			int lcp;
			if(ISA[k]+1 >= ISA.length) lcp = LCP[ISA[k]];
			else	lcp = Math.max(LCP[ISA[k]], LCP[ISA[k]+1]);			
			answer = Math.max(answer, Math.min(right-k+1, lcp));
		}
		return answer;
		
	}		
	
	
	public static int[] nucleotides2IntArray(String targetString){
		int targetLen = targetString.length();
		int[] targetIntArray = new int[targetLen];
		for(int i = 0; i < targetLen; i++)
			targetIntArray[i] = encode_char(targetString.charAt(i));
		return targetIntArray;
	}
	
	public static int encode_char(char c){
		// modified so that $ is the minimum value 0.
		switch (c){
		case '$':	return 0;
		case 'A':	return 1;
		case 'C':	return 2;
		case 'G':	return 3;
		case 'T':	return 4;
		default:	return 5;
		}
	}	

	public static int[] generateRandomNumbers(int targetLen){
		// Generate an array in which elements are selected from {1,2,3,4} at random.
		int[] target = new int[targetLen];
		// modified so that $ is the minimum value 0.
		for(int i=0; i<targetLen; i++) target[i] = (int)Math.floor(Math.random()*4) + 1;
		return target;
	}
	

	
	public static void main(String[] args) {
		int[] S_LS, S_KS, query;
		
		String s = "ATAATACGATAATAA";
		S_LS = nucleotides2IntArray(s+"$");
		S_KS = nucleotides2IntArray(s+"$$$");
		query = nucleotides2IntArray("ATA");

		// Generate random inputs.
/*		
		int targetLen=500000;			
		S_LS = generateRandomNumbers(targetLen+1); S_LS[targetLen]=0;
		S_KS = new int[targetLen+3];
		for(int i=0; i<targetLen; i++) 
			S_KS[i] = S_LS[i];
		S_KS[targetLen] = S_KS[targetLen+1] = S_KS[targetLen+2] = 0;
		query = generateRandomNumbers(8);
*/
		
		// Suffix array construction by Larsson-Sadakane algorithm.
		int[] SA_LS = new int[S_LS.length];
		suffixArray_LS(S_LS, SA_LS);

		// Suffix array construction by Karkkainen-Sanders algorithm.	
		int[] SA_KS = new int[S_KS.length];			
		suffixArray_KS(S_KS, S_KS.length-3, SA_KS, 5);

		// Check if the two algorithms output the same suffix array.
		for(int i=0; i<SA_LS.length-1; i++)
			if(SA_LS[i+1] != SA_KS[i]){
				System.out.println("The suffix arrays built by the two algorithms disagree.");
				break;
			}
		System.out.println("The suffix arrays built by the two algorithms agree.");

		// Use the suffix array built by Larsson-Sadakane algorithm.
		int[] S = S_LS;
		int[] SA = SA_LS;
		
		// Simple but naive binary search algorithm for the query.
		int qPos = searchLeftmost(S, SA, query);
		if(qPos > -1)
			System.out.println("The naive algorithm find the leftmost boundary of the block of the query occurrences at "+qPos+" in SA.");
		else
			System.out.println("The naive algorithm does not find the query.");		
		
		// Build lcp array. 
		int[] LCP = new int[S.length];
		int[] LCP_AUX = new int[S.length-1];
		for(int i=0; i<LCP_AUX.length; i++) LCP_AUX[i]=-1;
		buildLcp(S, SA, LCP, LCP_AUX);
	
		// Advanced binary search algorithm for the query with lcp information.
		qPos = searchLeftmost(S, SA, query, LCP, LCP_AUX);
		if(qPos > -1)
			System.out.println("The advanced algorithm find the leftmost boundary of the block of the query occurrences at "+qPos+" in SA.");
		else
			System.out.println("The advanced algorithm does not find the query.");
		
		// Build the list of occurrence frequencies of all k-mers.
		int k_mer = query.length;
		int[] freq = occurrence_frequency(SA, LCP, k_mer);
/*		
		System.out.println("Occurrence frequency of all "+k_mer+"-mers:");
		for(int i=0; i<freq.length; i++)
			System.out.print(freq[i]+" ");
		System.out.println();
*/
		
		// Compute the longest common factor.
		int[] ISA = new int[SA.length];
		for(int i=0; i<ISA.length; i++) ISA[SA[i]]=i; // Build the inverse suffix array.		
		System.out.println("The longest common factor of the query is "+lcf(ISA, LCP, 0, query.length-1)+".");
	}	
}

