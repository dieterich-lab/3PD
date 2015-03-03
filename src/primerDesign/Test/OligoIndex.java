/**
 * 
 */
package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;

import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import primerDesign.util.SlimFastaParser;
import primerDesign.util.SlimGenome;
import cern.colt.list.ObjectArrayList;

/**
 * @author froehler
 *
 */
public class OligoIndex implements Serializable{

	private static final long serialVersionUID = -1523326123131907541L;
	private HashMap<Integer, ObjectArrayList> positions;
	private SlimGenome genome;
	private static final String storeSuffix = ".OligoIndex.ser";
	
	public static void main(String[] args) throws NumberFormatException, IOException{
		if(args[0].equals("create")){
			SimpleTimer timer = new SimpleTimer();
			
			System.out.println("Creating index for sequence(s): " + args[1]);
			OligoIndex index = new OligoIndex(args[1], Integer.parseInt(args[2]));
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.println("Serializing index");
			index.serialize(new File(args[1]));
			System.out.println(" - done in " + timer.getTimeString());
			
		}
		else if(args[0].equals("query")){
			SimpleTimer timer = new SimpleTimer();
			
			System.out.print("Deserializing index");
			OligoIndex index = OligoIndex.deserialize(args[1]);
			System.out.println(" - done in " + timer.getTimeString());
			
			System.out.println("Found: " + index.getHits(args[2]).size() + " hits for query " + args[2] + " in " + timer.getTimeString());
		}
	}
	
	public OligoIndex(String filename, int xMerSize) throws IOException{
		this.positions = new HashMap<Integer, ObjectArrayList>();
		this.genome = new SlimGenome();
		
		SlimFastaParser parser = new SlimFastaParser(new File(filename));
		SimpleContigImpl contig;
		int hash;
		while(parser.hasNextContig()){
			contig = parser.parseNextContigIgnoreCase();
			this.genome.addContig(contig);
			ObjectArrayList temp;
			for(int i=0; i<contig.getSequenceLength() - xMerSize; i++){
				// add forward hits
				hash = contig.getSubsequence(i, i + xMerSize).hashCode();
				if(this.positions.containsKey(hash)){
					temp = this.positions.get(hash);
					temp.add(new IndexHitImpl(new SimpleContigImpl(contig.getID()), i, true));
					this.positions.put(hash, temp);
				}
				else{
					ObjectArrayList current = new ObjectArrayList();
					current.add(new IndexHitImpl(new SimpleContigImpl(contig.getID()), i, true));
					this.positions.put(hash, current);
				}
				// add reverse hits
				hash = SeqTools.revcompDNA(contig.getSubsequence(i, i + xMerSize).toCharArray()).hashCode();
				if(this.positions.containsKey(hash)){
					temp = this.positions.get(hash);
					temp.add(new IndexHitImpl(new SimpleContigImpl(contig.getID()), i, false));
					this.positions.put(hash, temp);
				}
				else{
					ObjectArrayList current = new ObjectArrayList();
					current.add(new IndexHitImpl(new SimpleContigImpl(contig.getID()), i, false));
					this.positions.put(hash, current);
				}
			}
		}
	}
	
	public ObjectArrayList getHits(String query){
		return this.positions.get(query.hashCode());
	}
	
	public void serialize(File filename){
		try{
			//ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename + storeSuffix));
			out.writeObject(this);
			out.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static OligoIndex deserialize(String filename){
		OligoIndex result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (OligoIndex) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	class Xmer implements Comparable{
		private int hash;
		private ObjectArrayList hits;
		
		public Xmer(int hash, ObjectArrayList hits){
			this.hash = hash;
			this.hits = hits;
		}
		
		public void addHit(IndexHitImpl hit){
			if(this.hits == null) throw new IllegalStateException("Hit list not initialized!");
			this.hits.add(hit);
		}
		
		public ObjectArrayList getHits(){
			return this.hits;
		}

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Object o) {
			Xmer other = (Xmer) o;
			if(this.hash < other.hash) return -1;
			else if(this.hash > other.hash) return +1;
			else return 0;
		}		
	}
}
