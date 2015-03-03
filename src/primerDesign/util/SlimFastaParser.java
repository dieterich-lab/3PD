package primerDesign.util;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.text.NumberFormat;
import java.util.regex.Pattern;

/**
 * A simple, lightweight fasta parser.
 * 
 * @author Sebastian fršhler
 *
 */
public class SlimFastaParser {
	private StringBuffer currentSequence;
	private static final Pattern fastaStart = Pattern.compile(">.*");
	private BufferedReader reader;
	private String line;
	private SimpleContigImpl currentContig = null;
	private static final String HEADER_SEPARATOR = "\\s+";
	
	/**
	 * Creates a new fasta parser on a given file containing zero or more fasta sequences.
	 * 
	 * @param file the file containing the fasta sequences
	 * 
	 * @throws IOException
	 */
	public SlimFastaParser(File file) throws IOException{
		
		this.reader = new BufferedReader(new FileReader(file));
		this.line = reader.readLine();
	}
	
	/**
	 * Creates a new fasta parser parsing a given string containing one or more fasta sequences.
	 * 
	 * Each fasta header has to start in a new line, each sequence must follow its fasta header and may contain of several lines of sequence.
	 * 
	 * @param fastaSequence the sequence to parse - in fasta format
	 * 
	 * @throws IOException
	 */
	public SlimFastaParser(String fastaSequence) throws IOException{
		this.reader = new BufferedReader(new StringReader(fastaSequence));
		this.line = reader.readLine();
	}
	
	/**
	 * Parses the next contig.
	 * 
	 * The contig name is supposed to be the first 'word', additional comments etc. can be added but have to be separated by whitespace each!
	 * 
	 * All characters in the sequence are treated to be upper-case!
	 * 
	 * @return a simple contig representing this contig
	 * 
	 * @throws IOException
	 */
	public SimpleContig parseNextContig() throws IOException{
		if(hasNextContig()){
			this.currentSequence = new StringBuffer();
			this.currentContig = new SimpleContigImpl(this.line.trim().substring(1).split(HEADER_SEPARATOR)[0]);
			
			for(this.line = this.reader.readLine(); !(this.line == null) && !SlimFastaParser.fastaStart.matcher(this.line).matches(); this.line = this.reader.readLine()){
				this.currentSequence.append(line.trim());
			}
			char[] seq = new char[this.currentSequence.length()];
			for(int i=0; i<this.currentSequence.length(); i++){
				seq[i] = this.currentSequence.charAt(i);
			}
			this.currentContig.setSequence(seq);
			
			return this.currentContig;			
		}
		else throw new EOFException("There are no more Contigs in this file!");
	}
	
	/**
	 * Parses the next contig ignoring case.
	 * 
	 * The contig name is supposed to be the first 'word', additional comments etc. can be added but have to be separated by whitespace each!
	 * 
	 * All characters in the sequence are treated to be upper-case!
	 * 
	 * @return a simple contig representing this contig
	 * 
	 * @throws IOException
	 */
	public SimpleContigImpl parseNextContigIgnoreCase() throws IOException{
		if(hasNextContig()){
			this.currentSequence = new StringBuffer();
			this.currentContig = new SimpleContigImpl(line.trim().substring(1).split(HEADER_SEPARATOR)[0]);
			
			for(this.line = this.reader.readLine(); !(this.line == null) && !SlimFastaParser.fastaStart.matcher(this.line).matches(); this.line = this.reader.readLine()){
				this.currentSequence.append(line.trim().toUpperCase());
			}
			char[] seq = new char[this.currentSequence.length()];
			for(int i=0; i<this.currentSequence.length(); i++){
				seq[i] = this.currentSequence.charAt(i);
			}
			this.currentContig.setSequence(seq);
			
			return this.currentContig;			
		}
		else throw new EOFException("There are no more Contigs in this file!");
	}
	
	/**
	 * Checks whether more contigs can be parsed.
	 * 
	 * @return true iff more contigs are available to the parser
	 * 
	 * @throws IOException
	 */
	public boolean hasNextContig() throws IOException{
		while(this.line != null && this.line.equals("")){
			this.line = this.reader.readLine();
		}
		if(this.line != null && SlimFastaParser.fastaStart.matcher(this.line).matches()) return true;
		else return false;
	}
	
	public static void main(String[] args) throws IOException{
		File file = new File("/Users/froehler/TestGenomes/H.sapiens/TestHuman.fa");
		NumberFormat format = NumberFormat.getInstance();
		Runtime runtime = Runtime.getRuntime();
		
		System.err.println("Memory used: " + format.format(runtime.totalMemory()-runtime.freeMemory()) + " bytes");
		
		SlimFastaParser parser = new SlimFastaParser(file);
		
		System.err.println("Memory used: " + format.format(runtime.totalMemory()-runtime.freeMemory()) + " bytes");
		
		while(parser.hasNextContig()){
			SimpleContig contig = parser.parseNextContig();
			System.out.println("Read contig " + contig.getID() + " length " + format.format(contig.getSequenceLength()));
			System.err.println("Memory used: " + format.format(runtime.totalMemory()-runtime.freeMemory()) + " bytes");
			System.gc();
		}
		
//		SlimFastaParser parser = new SlimFastaParser(">test1\nATGC\nATGC\n>test2\nGCGC");
//		while(parser.hasNextContig()){
//			System.out.print(parser.parseNextContig().toFastaString());
//		}
	}
}
