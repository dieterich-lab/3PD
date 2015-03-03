/**
 * 
 */
package primerDesign.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.HashMap;

/**
 * A wrapper class encapsulating a genome consisting of contigs with name and sequence each.
 * 
 * @author Sebastian Fr�hler
 *
 */
public class SlimGenome implements Serializable{
	private static final long serialVersionUID = -1L;

	private static final String storeSuffix = ".SlimGenome.ser";

	private HashMap<String, SimpleSlimContig> contigs;
	
	/** 
	 * Initializes a slim genome.
	 */
	public SlimGenome(){
		this.contigs = new HashMap<String, SimpleSlimContig>();
	}
	
	/**
	 * Initializes a slim genome.
	 * 
	 * @param contigs the contigs of the slim genome
	 */
	public SlimGenome(SimpleSlimContig[] contigs){
		this.contigs = new HashMap<String, SimpleSlimContig>();
		for(SimpleSlimContig contig : contigs){
			this.contigs.put(contig.getID(), contig);
		}
	}
	
	/**
	 * Adds a contig to a slim genome.
	 * 
	 * @param contig the contig to be added
	 */
	public void addContig(SimpleSlimContig contig){
		this.contigs.put(contig.getID(), contig);
	}
	
	/**
	 * Returns a contig of a slim genome.
	 * 
	 * @param name the name of the contig to return
	 * 
	 * @return a contig of a slim genome
	 */
	public SimpleContig getContig(String name){
		return this.contigs.get(name);
	}
	
	public void serialize(File filename){
		try{
			//ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename + storeSuffix));
			out.writeObject(this);
			out.close();
			System.out.println("finished.");

		}catch(Exception e){
			e.printStackTrace();
		}
	}

	public static SlimGenome deserialize(String filename){
		SlimGenome result = null;
		try{
			//ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)));
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
			result = (SlimGenome) in.readObject();
		}catch(Exception e){
			e.printStackTrace();
		}
		return result;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length < 2){
			System.out.println("Usage: SlimGenome create <FASTA FILE>");
			System.exit(1);
		}
		if(args[0].equals("create")){
			SlimGenome genome = new SlimGenome();
			SlimFastaParser parser = new SlimFastaParser(new File(args[1]));
			SimpleContigImpl current;
			while(parser.hasNextContig()){
				current = parser.parseNextContigIgnoreCase();
				genome.addContig(new SimpleSlimContig(current.getID(), current.getSequence()));
			}
			genome.serialize(new File(args[1]));
		}
		else{
			System.out.println("Usage: SlimGenome create <sequence file(s)>");
		}
	}
}
