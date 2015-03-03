package primerDesign.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;

import org.biojava.bio.program.ssaha.CompactedDataStoreFactory;
import org.biojava.bio.program.ssaha.DataStore;
import org.biojava.bio.program.ssaha.SearchListener;
import org.biojava.bio.program.ssaha.SequenceStreamer;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.DNAAmbPack;
import org.biojava.bio.symbol.Packing;
import org.biojavax.bio.seq.io.FastaFormat;

import primerDesign.dsc.DNASequenceIndex;
import primerDesign.util.Constants;

/**
 * test of the SSAHA algorithm.
 * 
 * @author froehler
 *
 */
public class SSAHATest implements DNASequenceIndex{
	private boolean includesScanRegion;
	private int wordSize;
	private DataStore store;
	private SymbolTokenization tokenization;

	public void createIndex(String sequence, int wordSize, boolean includesScanRegion) {
		int stepSize = 1;
		int threshold = Integer.MAX_VALUE;
		this.wordSize = wordSize;
		this.includesScanRegion = includesScanRegion;
		
		// 2Do: create custom alphabet
		try {
			this.tokenization = (SymbolTokenization) AlphabetManager.alphabetForName("DNA").getTokenization("token");
		} catch (Exception e) {
			e.printStackTrace();
		}
		SequenceStreamer.FileStreamer streamer = new SequenceStreamer.FileStreamer(new FastaFormat(), this.tokenization, new File(sequence));
		SequenceDB db;
		try {
			db = SeqIOTools.readFasta(new FileInputStream(sequence), DNATools.getDNA());
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		CompactedDataStoreFactory factory = new CompactedDataStoreFactory();
		File storeFile = new File("index-file.ssaha");
		Packing packing = new DNAAmbPack();
		try {
			//this.store = factory.buildDataStore(storeFile, db, packing, wordSize, threshold);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public Integer[] searchMatchPositionsInIndex(String searchString) {
		// TODO Auto-generated method stub
		return null;
	}

	public int searchNbMatchesInIndex(String searchString) {
		try {
			this.store.search("Primer", DNATools.createDNA(searchString), new SearchListener.Echo(new PrintStream(new FileOutputStream("Results.ssaha"))));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return 0;
	}

	public int getMaxWordSize() {
		return this.wordSize;
	}

	public String getSequence() {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean includesScanRegion() {
		return this.includesScanRegion;
	}
	
	public static void main(String[] args){
		SSAHATest test = new SSAHATest();
		test.createIndex("src/ensembl.txt", Constants.MAX_PRIMER_LENGTH, true);
	}
}
