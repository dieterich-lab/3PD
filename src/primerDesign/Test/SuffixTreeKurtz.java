package primerDesign.Test;

import primerDesign.dsc.DNASequenceIndex;


/** 
 * Implements a suffix tree ('SLLI') as proposed by S. Kurtz et.al. (1999).
 * 
 * Ref: S. Kurtz et.al: Reducing the Space Requirement of Suffix Trees. Software - Practice and Experience 29(13), 1149-1171 (1999)
 * 
 * @author froehler
 *
 */
public class SuffixTreeKurtz implements DNASequenceIndex{
	private int[] t_leaf;
	private int leafCount;
	private int[][] t_branch;  // format: <branching node number> <parameter of branching node>, e.g.: t_branch[i][FIRST_CHILD] returns the index of the first child of branching node i
	private int branchCount;
	private String sequence;
	private boolean includesScanRegion;
	private int maxWordLength;
	private int sequenceLength;
	
	// specify allocation of inner array of t_branch[][]
	private static int FIRST_CHILD = 0;
	private static int BRANCH_BROTHER = 1;	// the right branch brother of current node
	private static int DEPTH = 2;			// the depth of current node
	private static int HEAD_POSITION = 3;	// the head position of current node
	private static int SUFFIX_LINK = 4;		// link to the next suffix
	
	// specify the undefined value in t_branch and t_leaf
	private static int UNDEFINED = 0;
	private static int ROOT = 0;
	private static int HEAD_POS_ROOT = 1;

	public void createIndex(String sequence, int maxWordLength, boolean includesScanRegion) {
		this.sequence = sequence;
		this.sequenceLength = sequence.length();
		this.t_leaf = new int[sequenceLength + 1];
		this.leafCount = 0;
		this.t_branch = new int[sequenceLength][5];
		this.branchCount = 0;
		this.maxWordLength = maxWordLength;
		this.includesScanRegion = includesScanRegion;
		
		// init the leaf inserting the root
		setFirstChild(ROOT, UNDEFINED);
		setBranchBrother(ROOT, UNDEFINED);
		setDepth(ROOT, 0);
		setHeadPosition(ROOT, 1);
		setSuffixLink(ROOT, UNDEFINED);
		this.branchCount++;
		
		// extend the leaf by each suffix
		int current_branching_node;
		int current_leaf_node = leafCount++; // offset, leaf counts are valid in interval [1:sequenceLength+1]
		int current_node;
		for(int i=0; i<sequenceLength; i++){
			current_branching_node = ROOT;
			// insert new leaf
			current_leaf_node = leafCount++;
			setTleaf(current_leaf_node, UNDEFINED);
			// update t_branch if necessary
			if(getFirstChild(ROOT) == UNDEFINED){
				setFirstChild(ROOT, -current_leaf_node);
			}
			else{
				boolean inserted = false;
				current_node = getFirstChild(ROOT);
				while(!inserted){
					int headPos;
					if(current_node > UNDEFINED){
						// if current node is branching node
						int[] indices = getBranchingNodeLabelIndeces(current_node, ROOT);
						headPos = computeHeadPosition(indices[0], current_leaf_node, indices[1]);
						if(headPos > HEAD_POS_ROOT && headPos < getHeadPosition(current_node)){
							// insert new branching node before current_node, update references
							int newBranchNode = branchCount++;
							setFirstChild(newBranchNode, current_node);
							setBranchBrother(newBranchNode, UNDEFINED);
							setDepth(newBranchNode, headPos-1);
							setHeadPosition(newBranchNode, UNDEFINED);
							setSuffixLink(newBranchNode, UNDEFINED);
							
							inserted = true;
						}
						else if(headPos > HEAD_POS_ROOT && headPos == getHeadPosition(current_node)){
							// update references of current node OR check whether to insert new branching node AFTER current_node
						}
						else current_node = getBranchBrother(current_node);
					}else{
						// if current node is leaf node
						headPos = computeHeadPosition(current_node, current_leaf_node, sequenceLength);
						if(headPos == HEAD_POS_ROOT && getTleaf(current_node) == UNDEFINED) inserted = true;  // simply add another leaf to the root
						else if(headPos > HEAD_POS_ROOT){
							// insert new branching node
							int newBranchNode = branchCount++;
							setFirstChild(newBranchNode, -current_node);
							setBranchBrother(newBranchNode, UNDEFINED);
							setDepth(newBranchNode, headPos-1);
							setHeadPosition(newBranchNode, UNDEFINED);
							setSuffixLink(newBranchNode, UNDEFINED);
							inserted = true;
						}
						else if(headPos == HEAD_POS_ROOT && getTleaf(current_node) != UNDEFINED) current_node = getTleaf(current_node); // check next sibling of current node (= next child leaf of root)
						else throw new IllegalStateException("Unhandled case!"); // should NOT happen
					}
				}
			}
		}
	}
	
	// some convenience methods and security checks
	private int getTleaf(int leaf){
		if(leaf <1 || leaf>sequenceLength) throw new IllegalArgumentException("Leaf numbers must be in interval [1;sequenceLength]");
		return t_leaf[Math.abs(leaf)];
	}
	
	private void setTleaf(int leaf, int tLeaf){
		if(leaf <1 || tLeaf<1 || leaf>sequenceLength || tLeaf>sequenceLength) throw new IllegalArgumentException("Leaf numbers and t_leaf must be in interval [1;sequenceLength]");
		t_leaf[leaf] = tLeaf;
	}
	
	private int getFirstChild(int branchingNode){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		return t_branch[branchingNode][FIRST_CHILD];
	}
	
	private void setFirstChild(int branchingNode, int firstChild){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		t_branch[branchingNode][FIRST_CHILD] = firstChild;
	}
	
	private int getBranchBrother(int branchingNode){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		return t_branch[branchingNode][BRANCH_BROTHER];
	}
	
	private void setBranchBrother(int branchingNode, int branchBrother){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		t_branch[branchingNode][BRANCH_BROTHER] = branchBrother;
	}
	
	private int getDepth(int branchingNode){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		return t_branch[branchingNode][DEPTH];
	}
	
	private void setDepth(int branchingNode, int depth){
		if(branchingNode<0 || depth<0 || branchingNode>=sequenceLength || depth>sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1] and depth must be in interval [1;sequenceLength]!");
		t_branch[branchingNode][DEPTH] = depth;
	}
	
	private int getHeadPosition(int branchingNode){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		return t_branch[branchingNode][HEAD_POSITION];
	}
	
	private void setHeadPosition(int branchingNode, int headPosition){
		if(branchingNode<0 || headPosition<1 || branchingNode>=sequenceLength || headPosition>sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1], head position must be in interval [1;sequenceLength]!");
		t_branch[branchingNode][HEAD_POSITION] = headPosition;
	}
	
	private int getSuffixLink(int branchingNode){
		if(branchingNode<0 || branchingNode>=sequenceLength) throw new IllegalArgumentException("branchingNode must be in interval [0;sequenceLength-1]!");
		return t_branch[branchingNode][SUFFIX_LINK];
	}
	
	private void setSuffixLink(int branchingNode, int suffixLink){
		if(branchingNode<0 || suffixLink<0 || branchingNode>=sequenceLength || suffixLink >=sequenceLength) throw new IllegalArgumentException("branchingNode and suffix link must be in interval [0;sequenceLength-1]!");
		t_branch[branchingNode][SUFFIX_LINK] = suffixLink;
	}
	
	/**
	 * Retrieves the label of a branching node in a suffix tree.
	 * 
	 * @param nodeNumber the number of the node for which the label is to be computed
	 * @param predecessor the predecessor of node 'nodeNumber'
	 * @return the label of the branching node 'nodeNumber'
	 */
	private String getBranchingNodeLabel(int nodeNumber, int predecessor){
		if(nodeNumber<0 || predecessor<0 || nodeNumber>=sequenceLength || predecessor>=sequenceLength) throw new IllegalArgumentException("nodeNumber and predecessor must be in interval [0;sequenceLength-1]!");
		int[] indices = getBranchingNodeLabelIndeces(nodeNumber, predecessor);
		return sequence.substring(indices[0], indices[0]+indices[1]);
	}
	
	/**
	 * Retrieves the indices of the label of the branching node, format: start, stop+1
	 * @param nodeNumber the branching node number to compute the label for
	 * @param predecessor the predecessor branching node if this node
	 * @return an int[] containing label start, label stop+1
	 */
	private int[] getBranchingNodeLabelIndeces(int nodeNumber, int predecessor){
		if(nodeNumber<0 || predecessor<0 || nodeNumber>=sequenceLength || predecessor>=sequenceLength) throw new IllegalArgumentException("nodeNumber and predecessor must be in interval [0;sequenceLength-1]!");
		int[] indeces = new int[2];
		indeces[0] = t_branch[nodeNumber][HEAD_POSITION] + t_branch[predecessor][DEPTH];
		indeces[1] = t_branch[nodeNumber][DEPTH] - t_branch[predecessor][DEPTH];
		return indeces;
	}

	/**
	 * Retrieves the label of a leaf node in a suffix tree.
	 * 
	 * @param leafNumber the number of the leaf node for which the label is to be computed
	 * @param predecessor the predecessor of 'leafNumber'
	 * @return the label of the leaf node 'leafNumber'
	 */
	private String getLeafNodeLabel(int leafNumber, int predecessor){
		if(leafNumber<1 || predecessor<0 || leafNumber>sequenceLength || predecessor>=sequenceLength) throw new IllegalArgumentException("leafNumber must be in interval [1;sequenceLength], predecessor must be in interval [0;sequenceLength-1]!");
		int start = getLeafNodeLabelStartIndex(leafNumber, predecessor);
		return sequence.substring(start, this.sequenceLength);
	}
	
	/**
	 * Retrieves the start position index of the label of a leaf node.
	 * 
	 * @param leafNumber the leaf number to compute the label for
	 * @param predecessor the predecessor of this leaf
	 * @return the start position index of the label of a leaf node
	 */
	private int getLeafNodeLabelStartIndex(int leafNumber, int predecessor){
		if(leafNumber<1 || predecessor<0 || leafNumber>sequenceLength || predecessor>=sequenceLength) throw new IllegalArgumentException("leafNumber must be in interval [1;sequenceLength], predecessor must be in interval [0;sequenceLength-1]!");
		int i = leafNumber + t_branch[predecessor][DEPTH];
		return i;
	}
	
	/**
	 * Computes the number of common characters between first and second.
	 * 
	 * @param first the start index of the first label
	 * @param second the start index of the second label
	 * @param stop compare characters up to (but not including) character 'stop'
	 * 
	 * @return the number of common characters between first and second
	 */
	private int computeHeadPosition(int first, int second, int stop){
		if(first<0 || second<0 || first>=sequenceLength || second>=sequenceLength || stop<1 || stop>sequenceLength) throw new IllegalArgumentException("Parameter first or second or stop out of range!");
		int headPosition = 1;
		// offset -1 since sequence characters have values [0:length-1]
		int firstChar = Math.abs(first)-1;
		int secondChar = Math.abs(second)-1;
		int end = Math.min(stop, sequenceLength-Math.max(firstChar, secondChar));
		for(int i=0; i<end; i++){
			if(this.sequence.charAt(firstChar+i) == this.sequence.charAt(secondChar+i)) headPosition++;
			else break;
		}
		return headPosition;
	}
	
	public Integer[] searchMatchPositionsInIndex(String searchString) {
		searchString = searchString.toUpperCase();
		// TODO Auto-generated method stub
		return null;
	}

	public int searchNbMatchesInIndex(String searchString) {
		return searchMatchPositionsInIndex(searchString).length;
	}

	public String getSequence() {
		return this.sequence;
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}

	public int getMaxWordSize() {
		return this.maxWordLength;
	}
	
	/**
	 * Returns a textual representation of the tree for debugging purposes.
	 *
	 */
	private void printTree(){
		System.out.println("The Kurtz et.al. suffix tree for the sequence");
		System.out.println();
		System.out.println("\t" + sequence);
		System.out.println();
		System.out.println("T_leaf:");
		System.out.println("#\tt_leaf(#)");
		for(int i=1; i<t_leaf.length; i++){
			System.out.println(i + "\t" + t_leaf[i]);
		}
		System.out.println();
		System.out.println("T_branch");
		System.out.println("#\tfirstchild\tbranchbrother\tdepth\theadpos\tsuffixlink");
		for(int i=0; i<t_branch.length; i++){
			System.out.println(i + "\t" + t_branch[i][FIRST_CHILD] + "\t" + t_branch[i][BRANCH_BROTHER] + "\t" + t_branch[i][DEPTH] + "\t" + t_branch[i][HEAD_POSITION] + "\t" + t_branch[i][SUFFIX_LINK]);
		}
	}
	
	public static void main(String[] args){
		SuffixTreeKurtz tree = new SuffixTreeKurtz();
		tree.createIndex("ABAB", 1, true);
		tree.printTree();
	}
}
