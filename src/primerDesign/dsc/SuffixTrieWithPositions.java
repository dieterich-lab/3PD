package primerDesign.dsc;

import java.util.ArrayList;
import java.util.HashSet;

import org.biojava.bio.symbol.IllegalSymbolException;

import cern.colt.list.ObjectArrayList;

public class SuffixTrieWithPositions extends SuffixTrie {
	private ObjectArrayList alphabet;
	private ArrayList<Integer> counts;
	private SuffixNode root;
	
	public SuffixTrieWithPositions(ObjectArrayList alphabet) {
		super(alphabet);
		this.alphabet = alphabet;
		this.counts = new ArrayList<Integer>();
		this.root = new SimpleNodeWithPositions(alphabet.size(), 0);
	}

	/**
     * Add the match positions for all motifs with length of up to <code>window</code>
     * to this trie.
     *
     * @param sList a SymbolList whose motifs should be added to the trie.
     * @param window The maximum motif length to count.
     */

  public void addSymbols(String sList, int window) throws IllegalSymbolException {	  
    SuffixNode [] buf = new SuffixNode[window];
    int [] counts = new int[window];
    for(int i = 0; i < window; i++) {
      buf[i] = getRoot();
    }
    
    for(int p = 0; p < sList.length(); p++) {
      char s = sList.charAt(p);
      buf[p % window] = getRoot();
      for(int i = 0; i < window; i++) {
        int pi = (p + i) % window;
        if(buf[pi] != null) {
        	int depth = buf[pi].getDepth();
          buf[pi] = getChild(buf[pi], s);
          if(buf[pi] != null) {
            counts[i]++;
            buf[pi].addPosition(p-depth);
          }
        }
      }
    }
    
    for(int i = 0; i < window; i++) {
      incCounts(i+1, counts[i]);
    } 
  }
  
  public SuffixNode getRoot() {
	    return this.root;
	  }
  
  public SuffixNode getChild(SuffixNode node, int i) {
	    if(!node.hasChild(i) && i < this.alphabet.size()) {
	      node.addChild(this, i, new SimpleNodeWithPositions(this.alphabet.size(), node.getDepth()+1));
	    }else if(i >= this.alphabet.size() || i < 0) throw new ArrayIndexOutOfBoundsException("Index " + i + " is out of bounds.");
	    return node.getChild(i);
	  }
  
  private static class SimpleNodeWithPositions extends SuffixNode {
	    /**
		 * 
		 */
		private static final long serialVersionUID = 7382918145922755541L;
		private SuffixNode [] child;
		private HashSet<Integer> matchPositions;
		private int depth;
	    
	    private SuffixNode [] childArray(SuffixTrie tree) {
	      if(child == null)
	        child = new SuffixNode[tree.getAlphabet().size()];
	      return child;
	    }
	    
	    public boolean isTerminal() {
	      return child == null;
	    }
	    
	    public boolean hasChild(int i) {
	      return child != null && child[i] != null;
	    }
	    
	    SuffixNode getChild(int i) {
	      if(hasChild(i))
	        return child[i];
	      else throw new IllegalArgumentException("Current node has no child " + i);
	    }
	    
	    void addChild(SuffixTrie tree, int i, SuffixNode n) {
	      childArray(tree)[i] = n;
	    }
	    
	    SimpleNodeWithPositions(int c, int depth) {
	      this.child = new SuffixNode[c];
	      this.matchPositions = new HashSet<Integer>();
	      this.depth = depth;
	    }
	    
	    public boolean equals(Object otherNode){
	    	SuffixNode other = (SuffixNode) otherNode;
	    	return this.matchPositions.equals(other.getPositions());
	    }

		public float getNumberOfMatches() {
			return this.matchPositions.size();
		}

		public void incCount() {
			throw new NoSuchMethodError("Class SimpleNodeWithPositions does NOT support position-less match counts!");
		}
		
		public void addPosition(int position){
			this.matchPositions.add(position);
		}
		
		public HashSet getPositions(){
			return this.matchPositions;
		}

		int getDepth() {
			return this.depth;
		}
  }
}
