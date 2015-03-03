package primerDesign.dsc;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;

import org.biojava.bio.symbol.IllegalSymbolException;

import cern.colt.list.ObjectArrayList;

/**
 * This code is an adaption of org.biojava.bio.symbol.SuffixTree and now contains positions of matches
 * in the sequence this tree was constructed from instead of just counts of motifs.
 *
 * @author Matthew Pocock
 * @author Thomas Down 
 * @author Sebastian Fršhler
 */

public class SuffixTrie {
	private ObjectArrayList alphabet;
	private SuffixNode root;
	private ArrayList<Integer> counts;

	public ObjectArrayList getAlphabet(){
		return this.alphabet;
	}
	
    /**
     * Return the node object which is the root of this suffix tree.
     * This represents the set of all motifs found in this tree.
     */
  public SuffixNode getRoot() {
    return this.root;
  }
  
    /**
     * Get a child of a SuffixTrie.SuffixNode, constructing a new
     * one if need be.  This method is here due to memory optimisations.
     */

  public SuffixNode getChild(SuffixNode node, char s){
    int index = this.alphabet.indexOf(s, true);
    if(index == -1) throw new IllegalArgumentException("Character " + s + " is no member of the alphabet over which this suffix tree was constructed.");
    return getChild(node, index);
  }
  
    /**
     * Get the n'th child of a node.
     */

  public SuffixNode getChild(SuffixNode node, int i) {
    if(!node.hasChild(i) && i < this.alphabet.size()) {
      node.addChild(this, i, new SimpleNode(this.alphabet.size()));
    }else if(i >= this.alphabet.size() || i < 0) throw new ArrayIndexOutOfBoundsException("Index " + i + " is out of bounds.");
    return node.getChild(i);
  }
  
    /**
     * Add a count for all motifs with length of up to <code>window</code>
     * to this trie.
     *
     * @param sList a SymbolList whose motifs should be added to the trie.
     * @param window The maximum motif length to count.
     */

  public void addSymbols(String sList, int window)
  throws IllegalSymbolException {	  
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
          buf[pi] = getChild(buf[pi], s);
          if(buf[pi] != null) {
            counts[i]++;
            buf[pi].incCount();
          }
        }
      }
    }
    
    for(int i = 0; i < window; i++) {
      incCounts(i+1, counts[i]);
    } 
  }
  
  protected void incCounts(int i, int c) {
	if(i<0) throw new IllegalArgumentException("Index i must be >= 0");
	if(c<0) throw new IllegalArgumentException("Count c must be >= 0");
    if(i < counts.size()) {
      Integer oldC = (Integer) counts.get(i-1);
      Integer newC = new Integer(oldC.intValue() + c);
      counts.set(i-1, newC);
    } else {
      counts.add(new Integer(c));
    }
  }
  
    /**
     * Return the length of the longest motif represented in this
     * SuffixTrie
     */

  public int maxLength() {
    return counts.size();
  }
  
    /**
     * Return the number of motifs of a given length encoded
     * in this SuffixTrie.
     */

  public int frequency(int length) {
    return ((Integer) counts.get(length - 1)).intValue();
  }
  
    /**
     * Construct a new SuffixTrie to contain motifs over the
     * specified alphabet.
     *
     * @param alphabet The alphabet of this SuffixTrie (must be finite).
     */

  public SuffixTrie(ObjectArrayList alphabet) {
	 this.alphabet = alphabet;
    this.counts = new ArrayList<Integer>();
    this.root = new SimpleNode(alphabet.size());
  }
  
  public boolean equals(Object otherTree){
	  SuffixTrie other = (SuffixTrie) otherTree;
	  return this.root.equals(other.root) && this.alphabet.equals(other.alphabet) && this.counts.equals(other.counts);
  }
  
  /**
   * A node in the suffix tree.
   * <p>
   * This class is realy stupid & delegates most work off to a SuffixTrie so
   * that it is as small (in memory-per-object terms) as possible.
   *
   * @author Matthew Pocock
   */
  public static abstract class SuffixNode implements Serializable {
      /**
       * Determine is this node is terminal (has no children).
       *
       * @return <code>true</code> if and only if this node has no children.
       */

    abstract public boolean isTerminal();
      
      /**
       * Determine if this node has a child corresponding to a given
       * alphabet number.
       */

    abstract public boolean hasChild(int i);

      /**
       * Return a number (usually, but not always, a motif count)
       * associated with this node of the tree.
       */

    abstract public float getNumberOfMatches();
    
    abstract public void incCount();
    
    abstract public void addPosition(int position);
    
    abstract public HashSet getPositions();

    abstract SuffixNode getChild(int i);
    abstract void addChild(SuffixTrie tree, int i, SuffixNode n);
    
    abstract int getDepth();
  }

  private static class SimpleNode extends SuffixNode {
    /**
	 * 
	 */
	private static final long serialVersionUID = 7382918145922755541L;
	private SuffixNode [] child;
	private int matchCount;
    
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
    
    SimpleNode(int c) {
      this.child = new SuffixNode[c];
      this.matchCount = 0;
    }
    
    public boolean equals(Object otherNode){
    	SuffixNode other = (SuffixNode) otherNode;
    	return this.matchCount == other.getNumberOfMatches();
    }

	public float getNumberOfMatches() {
		return this.matchCount;
	}

	public void incCount() {
		this.matchCount++;
		
	}

	public void addPosition(int position) {
		throw new NoSuchMethodError("Class SimpleNode does NOT support match positions!");
	}

	public HashSet getPositions() {
		throw new NoSuchMethodError("Class SimpleNode does NOT support match positions!");
	}

	int getDepth() {
		throw new NoSuchMethodError("Class SimpleNode does NOT support node depths!");
	}
  }
}
