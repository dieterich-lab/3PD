/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.io.File;
import java.rmi.Naming;
import java.rmi.RemoteException;
import java.util.HashMap;

import primerDesign.dsc.indexStructures.DNASequenceIndex;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Encapsulates a client to perform index queries on a remote index.
 * 
 * @author Sebastian Fršhler
 *
 */
public class IndexRemoteQueryClient implements DNASequenceIndex{
	private final String host;
	private final int port;
	private final String name;
	private String url;
	private RemoteIndexStructureSearch search;
	private boolean debug = true;
	
	/**
	 * Initializes a new remote index query and connects to a default remote index.
	 */
	public IndexRemoteQueryClient() {
		this.host = "abt4-cd-xeon"; // "abt4-cd-2ghz";
		this.port = 1099;
		this.name = "RemoteIndex";
		this.url = "rmi://" + host + ":" + port + "/" + name;
		this.connect();
	}
	
	/**
	 * Initializes a new remote index query and connects to a specified remote index.
	 * 
	 * @param host the host hosting the remote index
	 * @param port the port on the host
	 * @param name the name of the remote index to perform queries on
	 */
	public IndexRemoteQueryClient(String host, int port, String name){
		this.host = host;
		this.port = port;
		this.name = name;
		this.connect();
	}

	/**
	 * Connects to a  remote index.
	 *
	 */
	private void connect(){
		
		//System.setSecurityManager(new RMISecurityManager()); 
		
		this.url = "rmi://" + host + ":" + port + "/" + name;
		
		if(debug) System.out.println("Looking-up RemoteIndex " + this.url);
		
		try{
			this.search = (RemoteIndexStructureSearch) Naming.lookup(url);
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * Queries a remote index for hits of query 'query'.
	 * 
	 * @param query the query to the index
	 * 
	 * @return a list of index hits in the remote index: @see esaRmi.IndexHit
	 * 
	 * @throws RemoteException
	 */
	public ObjectArrayList findHitPositions(String query){
		ObjectArrayList result = null;
		try{
			result = this.search.findHitPositions(query);
		}
		catch(RemoteException e){
			e.printStackTrace();
		}
		return result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#findHitCount(java.lang.String)
	 */
	public int findHitCount(String sequence) {
		int result = -1;
		try{
			result = this.search.findHitPositions(sequence).size();
		}
		catch(RemoteException e){
			e.printStackTrace();
		}
		return result;
	}
	
	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#createIndex()
	 */
	public void createIndex() {
		throw new IllegalStateException("A query client can perform only query tasks!");
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getContig()
	 */
	public SimpleContigImpl[] getContig() {
		throw new IllegalStateException("Unsupported method!");
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getName()
	 */
	public String getName() {
		throw new IllegalStateException("Unsupported method!");
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#getStatistics()
	 */
	public HashMap<Character, Integer> getStatistics() {
		throw new IllegalStateException("Unsupported method!");
	}

	/* (non-Javadoc)
	 * @see primerDesign.dsc.indexStructures.DNASequenceIndex#serialize(java.io.File)
	 */
	public void serialize(File file) {
		throw new IllegalStateException("Unsupported method!");
	}
	
	public static void main(String[] args) throws RemoteException{
		String query = args[0];
		int numQueries = Integer.parseInt(args[1]);
		SimpleTimer timer = new SimpleTimer();
		
		System.out.print("Fetching reference to remote object");
		IndexRemoteQueryClient client = new IndexRemoteQueryClient();
		System.out.println(" - took " + timer.getTimeString());
		ObjectArrayList matches = null;
		
		System.out.println("Benchmarking match position lookup:");
		for(int i=1; i<=numQueries; i*=10){
			System.out.print("\tPerforming " + i + " queries");
			for(int j=0; j<i; j++){
				matches = client.findHitPositions(query);
				for(int k=0; k<matches.size(); k++){
				}
			}
			System.out.print(" - " + matches.size() + " matches");
			System.out.println(" - done in " + timer.getTimeString());
		}
		
		System.out.println("Benchmarking match count lookup:");
		int hits = 0;
		for(int i=1; i<=numQueries; i*=10){
			System.out.print("\tPerforming " + i + " queries");
			for(int j=0; j<i; j++){
				hits = client.findHitPositions(query).size();
			}
			System.out.print(" - " + hits + " matches");
			System.out.println(" - done in " + timer.getTimeString());
		}
	}
}
