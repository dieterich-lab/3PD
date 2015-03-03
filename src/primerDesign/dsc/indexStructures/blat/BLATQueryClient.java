package primerDesign.dsc.indexStructures.blat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.Random;
import java.util.Arrays;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.dsc.indexStructures.IndexHitImpl;
import primerDesign.util.SimpleSlimContig;
import primerDesign.util.SimpleTimer;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a JAVA client to query a BLAT index.
 * 
 * @author Sebastian Fr�hler
 *
 */
public class BLATQueryClient {
	private String CLIENT = "/bin/gfClient";
	//private String OPTIONS = "-minScore=10";
	private String HOST = "localhost";
	private String PORT = "10110";
	private String SEQDIR = "/Data/";
//	private static final Runtime RUNTIME = Runtime.getRuntime();
//	private String[] cmdarray;
//	private Process proc;
//	private OutputStream stream;
//	private BufferedWriter writer;
	private int sessionID;
	private static Random RANDOM = new Random();
	
	private static final int PRIMER_END_LENGTH = 5;

	public static void main(String[] args) throws IllegalAlphabetException, IllegalSymbolException, IOException, InterruptedException {	
		SimpleTimer timer = new SimpleTimer();
		NumberFormat format = NumberFormat.getInstance();
		
		BLATQueryClient book = new BLATQueryClient(args[0], args[1], args[2], args[3]);
		int numQueries = Integer.parseInt(args[4]);
		System.out.print("Querying " + format.format(numQueries) + " queries on port " + book.PORT);
		int hits = 03;
		
		System.out.println("Server up and running?: " + book.isAvailable());
		
//		PrimerSearchParameters params = new PrimerSearchParameters();
//		int min = params.getMIN_PRIMER_LENGTH();
//		int max = params.getMAX_PRIMER_LENGTH();
		for(int i=0; i<numQueries; i++){
			hits = book.getPrimings("Test", "TATCAATTTTA").size(); //SeqTools.getRandomPrimerSequence(11, 11));
		}
		System.out.println(" - done in " + timer.getTimeString());
		System.out.println("Found: " + hits + " hits!");
	}
	
	/**
	 * Sets the properties of a BLAT query client.
	 * 
	 * After iniializing this query client, one should assert that the server this client should connect is up and running by using the 'isAvailable' method!
	 * 
	 * @param client the absolute path to the 'gfClient' program
	 * @param host the adress or hostname of the host hosting the BLAT server
	 * @param port the port of the respective BLAT index on host 'host'
	 * @param sequenceDir the sequence directory where sequences of the indices can be found
	 */
	public BLATQueryClient(String client, String host, String port, String sequenceDir){
		this.CLIENT = client;
		this.HOST = host;
		this.PORT = port;
		this.SEQDIR = sequenceDir;
		this.sessionID = RANDOM.nextInt();
	}
	
	private BufferedReader getPrimingsInternal(String name, String sequence) throws IOException, InterruptedException{
		
		//String[] cmdarray = new String[]{"printf", "%s\n",  ">" + name,  sequence};
		BufferedWriter writer = new BufferedWriter(new FileWriter("/tmp/3CPrimerDesignQuery-" + this.sessionID + "-" + sequence + ".tmp"));
		writer.write(">" + name + "\n" + sequence);
		writer.close();
		
		String[] cmdarray2 = new String[]{CLIENT, "query" ,HOST, PORT, "/tmp/3CPrimerDesignQuery-" + this.sessionID + "-" + sequence + ".tmp", ">" ,"/tmp/3CPrimerDesignQuery-result-" + this.sessionID + "-" + sequence + ".tmp"};
		//String[] cmdarray2 = new String[]{CLIENT, OPTIONS, HOST, PORT, SEQDIR, "/dev/stdin", "/dev/stdout"};
		
		System.out.println(Arrays.toString(cmdarray2));
		
		
		ProcessBuilder builder = new ProcessBuilder("sh", cmdarray2.toString());
		builder.redirectOutput(new File("/tmp/3CPrimerDesignQuery-result-" + this.sessionID + "-" + sequence + ".tmp"));
		builder.redirectError(new File("/tmp/3CPrimerDesignQuery-result-" + this.sessionID + "-" + sequence + ".error"));
		Process p = builder.start(); // throws IOException
		
//		Process proc2 = Runtime.getRuntime().exec(cmdarray2);
//		proc2.waitFor();	
		
		//FileReader reader = new FileReader(new File("/tmp/3CPrimerDesignQuery-result-" + this.sessionID + "-" + sequence + ".tmp"));

		//BufferedReader reader = new BufferedReader(new InputStreamReader(proc2.getErrorStream()));
//		String line;
//		while((line = reader.r != null){
//			System.out.println(line);
//			if(line.startsWith("Sorry, the BLAT/iPCR server seems to be down")) throw new IOException("BLAT server unavailable at host: " + this.HOST + " and port: " + this.PORT);
//		}

		return new BufferedReader(new FileReader("/tmp/3CPrimerDesignQuery-result-" + this.sessionID + "-" + sequence + ".tmp"));
		
		// START TEST
//		cmdarray = new String[]{CLIENT, OPTIONS, HOST, PORT, SEQDIR, "/dev/stdin", "/dev/stdout"};
//		proc = RUNTIME.exec(cmdarray);
//		stream = proc.getOutputStream();
//		writer = new BufferedWriter(new OutputStreamWriter(stream));
//		writer.write(">" + name + "\n" + sequence);
//		writer.close();
//		stream.close();
//		
//		proc.waitFor();
//		return new Result(new BufferedReader(new InputStreamReader(proc.getInputStream())), proc);
		// END TEST

		//System.out.print("generating query");
//		Process proc = Runtime.getRuntime().exec(cmdarray);
//		InputStream in = proc.getInputStream();		
		//proc.waitFor();
		
		//String[] cmdarray2 = new String[]{CLIENT, " " ,OPTIONS, " ", HOST, " ", PORT, " ", SEQDIR, " ", "/dev/stdin", " ", "/dev/stdout"};
		//String[] cmdarray2 = new String[]{CLIENT, OPTIONS, HOST, PORT, SEQDIR, "/dev/stdin", "/dev/stdout"};

		//Process proc2 = Runtime.getRuntime().exec(CLIENT + " " + OPTIONS + " " + HOST + " " + PORT + " " + SEQDIR + " " + "/dev/stdin" + " " + "/dev/stdout");
		
//		BufferedWriter writer2 = new BufferedWriter(new OutputStreamWriter(proc2.getOutputStream()));
//		writer2.write(">" + name + "\n");
//		writer2.write(sequence);
//		writer2.close();
//		
//		proc2.waitFor();
//		return new BufferedReader(new InputStreamReader(proc2.getInputStream()));
		
		//System.out.println(" - Processing result");
//		OutputStream out = proc2.getOutputStream();
		//InputStream result = proc2.getInputStream();
		
//		InputStream errors = proc2.getErrorStream();
//		BufferedReader errorReader = new BufferedReader(new InputStreamReader(errors));
		
//		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(out));
//		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
//		
//		String line;
//		while((line = reader.readLine()) != null){
//			writer.write(line + "\n");
//			//System.err.println(line);
//		}
//		in.close();
//		out.close();
//		writer.close();
//		reader.close();

		
//		proc.destroy();
		
//		String tmp;
//		while((tmp = errorReader.readLine()) != null){
//			System.err.println(tmp);
//		}
		
		//return new BufferedReader(new InputStreamReader(result));
		
		
//		ObjectArrayList hits = new ObjectArrayList();
//		String ln;
//		while((ln = reader2.readLine()) != null){
//			String[] values = ln.split("\t");
//			if(values.length == 21 && !values[0].equals("match")){
//				// prune hits
//				int queryLength = Integer.parseInt(values[10]);
//				int queryStart = Integer.parseInt(values[11]);
//				int queryEnd = Integer.parseInt(values[12]);
//				
//				// prune misprimings with NO end priming
//				if(queryEnd < queryLength - 5) continue;
//				else hits.add(new IndexHitImpl(new SimpleContig(values[13]), Integer.parseInt(values[15]), values[8].equals("+") ? true : false));
//				//System.out.println(ln);
//			}
//		}
//		System.out.println("Found " + hits.size() + " hits for query " + sequence);
//		return hits;
	}
	
	private boolean isDangerousMispriming(int length, int qEnd){
		if(qEnd < length - PRIMER_END_LENGTH) return false;
		else return true;
	}
	
	/**
	 * Returns wehther a BLAT server this client should connect to is actually up and running on 'host' and 'port'.
	 * 
	 * @return true iff a BLAT server this client should connect to is actually up and running on 'host' and 'port'
	 */
	public boolean isAvailable(){
		try{
			getPrimingsInternal("Test", "ATGATGATGATGATGATGATG");
		}
		catch(IOException e){
			return false;
		}
		catch(InterruptedException e){
			return false;
		}
		return true;
	}
	
	/**
	 * Queries a given BLAT index, returns all primings where the end (the last five basepairs) of the sequence are included in the priming.
	 * 
	 * (By convention of BLAT standard parameters, each priming has at lest length 11)
	 * 
	 * @param name the name of the sequence
	 * @param sequence the sequence to scan for primings
	 * 
	 * @return the number of primings of sequence 'sequence'
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public int getNumPrimings(String name, String sequence) throws IOException, InterruptedException{
		BufferedReader reader = getPrimingsInternal(name, sequence);
		String ln;
		int primings = 0;
		while((ln = reader.readLine()) != null){
			String[] values = ln.split("\t");
			if(values.length == 21 && !values[0].equals("match")){
				// prune hits
				int queryLength = Integer.parseInt(values[10]);
//				int queryStart = Integer.parseInt(values[11]);
				int queryEnd = Integer.parseInt(values[12]);
				
				// prune misprimings with NO end priming
				if(isDangerousMispriming(queryLength, queryEnd)) primings++;
				else continue;
			}
		}
		return primings;
	}
	
	/**
	 * Queries a given BLAT index, returns all primings where the end (the last five basepairs) of the sequence are included in the priming.
	 * 
	 * (By convention of BLAT standard parameters, each priming has at lest length 11)
	 * 
	 * @param name the name of the sequence
	 * @param sequence the sequence to scan for primings
	 * 
	 * @return a list of IndexHitImpl objects, each representing one priming
	 * 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public ObjectArrayList getPrimings(String name, String sequence) throws IOException, InterruptedException{
		BufferedReader reader = getPrimingsInternal(name, sequence);
		
		ObjectArrayList hits = new ObjectArrayList();
		String ln;
		while((ln = reader.readLine()) != null){
			String[] values = ln.split("\t");
			if(values.length == 21 && !values[0].equals("match")){
				// prune hits
				int queryLength = Integer.parseInt(values[10]);
//				int queryStart = Integer.parseInt(values[11]);
				int queryEnd = Integer.parseInt(values[12]);
				
				// prune misprimings with NO end priming
				if(isDangerousMispriming(queryLength, queryEnd)){
					//hits.add(new IndexHitImpl(new SimpleContig(values[13].replace('_', ' ').substring(0, values[13].length() - 1)), Integer.parseInt(values[15]), values[8].equals("+") ? true : false));
					hits.add(new IndexHitImpl(new SimpleSlimContig(values[13]), Integer.parseInt(values[15]), values[8].equals("+") ? true : false));
				}
				else continue;
			}
		}
		//System.out.println("Found " + hits.size() + " hits for query " + sequence);
		return hits;
	}
	
	public String toString(){
		StringBuffer buffy = new StringBuffer();
		buffy.append(this.CLIENT + " " + " " + this.HOST + " " + this.PORT + " " + this.SEQDIR);
		return buffy.toString();
	}
}
