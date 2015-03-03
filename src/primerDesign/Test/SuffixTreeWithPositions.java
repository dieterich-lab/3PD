package primerDesign.Test;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * This code is an adaption of org.biojava.bio.symbol.SuffixTree and now contains positions of matches
 * in the sequence this tree was constructed from instead of just counts of motifs.
 *
 * @author Matthew Pocock
 * @author Thomas Down 
 * @author Sebastian Fršhler
 */

public class SuffixTreeWithPositions implements Serializable {
	private FiniteAlphabet alphabet;
	private SuffixNode root;
	private AlphabetIndex indexer;
	private List counts;
  
    /**
     * Return the Alphabet containing all Symbols which might be found in
     * this SuffixTrie.
     */

  public FiniteAlphabet getAlphabet() {
    return alphabet;
  }

    /**
     * Return the node object which is the root of this suffix tree.
     * This represents the set of all motifs found in this tree.
     */

  public SuffixNode getRoot() {
    return root;
  }
  
    /**
     * Get a child of a SuffixTrie.SuffixNode, constructing a new
     * one if need be.  This method is here due to memory optimisations.
     */

  public SuffixNode getChild(SuffixNode node, Symbol s)
  throws IllegalSymbolException {
    if(!getAlphabet().contains(s)) {
    	throw new IllegalArgumentException("Alphabet " + this.getAlphabet().getName() + " does not contain symbol " + s.getName());
      //return null;
    }
    int index = indexer.indexForSymbol(s);
    return getChild(node, index);
  }
  
    /**
     * Get the n'th child of a node.
     */

  public SuffixNode getChild(SuffixNode node, int i) {
    if(!node.hasChild(i)) {
      node.addChild(this, i, new SimpleNode(alphabet.size()));
    }
    return node.getChild(i);
  }
  
    /**
     * Add a count for all motifs with length of up to <code>window</code>
     * to this tree.
     *
     * @param sList a SymbolList whose motifs should be added to the
     *              tree.
     * @param window The maximum motif length to count.
     */

  public void addSymbols(SymbolList sList, int window)
  throws IllegalSymbolException {
    SuffixNode [] buf = new SuffixNode[window];
    int [] counts = new int[window];
    for(int i = 0; i < window; i++) {
      buf[i] = getRoot();
    }
    
    for(int p = 1; p <= sList.length(); p++) {
      Symbol s = sList.symbolAt(p);
      buf[p % window] = getRoot();
      for(int i = 0; i < window; i++) {
        int pi = (p + i) % window;
        if(buf[pi] != null) {
          buf[pi] = getChild(buf[pi], s);
          if(buf[pi] != null) {
            counts[i]++;
            buf[pi].addMatchPosition(p);
          }
        }
      }
    }
    
    for(int i = 0; i < window; i++) {
      incCounts(i+1, counts[i]);
    } 
  }
  
  protected void incCounts(int i, int c) {
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
     * @param alphabet The alphabet of this SuffixTrie (must be
     *                 finite).
     */

  public SuffixTreeWithPositions(FiniteAlphabet alphabet) {
    this.alphabet = alphabet;
    this.indexer = AlphabetManager.getAlphabetIndex(alphabet);
    this.counts = new ArrayList();
    this.root = new SimpleNode(alphabet.size());
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
       * index number.
       */

    abstract public boolean hasChild(int i);

      /**
       * Return a number (usually, but not always, a motif count)
       * associated with this node of the tree.
       */

    abstract public float getNumberOfMatches();
    
    abstract public void addMatchPosition(int position);
    
    abstract public Integer[] getAllMatchPositions();

    abstract SuffixNode getChild(int i);
    abstract void addChild(SuffixTreeWithPositions tree, int i, SuffixNode n);
  }

  private static class SimpleNode extends SuffixNode {
    /**
	 * 
	 */
	private static final long serialVersionUID = -6795489470718766429L;
	private SuffixNode [] child;
    //private FastVector matchPositions = new FastVector();
    private HashSet<Integer> matchPos;
    
    private SuffixNode [] childArray(SuffixTreeWithPositions tree) {
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
    
    public float getNumberOfMatches() {
      return matchPos.size();
    }
    
    public void addMatchPosition(int position){
    	this.matchPos.add(position);
    }
    
    public Integer[] getAllMatchPositions(){
    	return this.matchPos.toArray(new Integer[this.matchPos.size()]);
    }
    
    SuffixNode getChild(int i) {
      if(hasChild(i))
        return child[i];
      return null;
    }
    
    void addChild(SuffixTreeWithPositions tree, int i, SuffixNode n) {
      childArray(tree)[i] = n;
    }
    
    SimpleNode(int c) {
      child = new SuffixNode[c];
      matchPos = new HashSet<Integer>();
    }
  }
}
