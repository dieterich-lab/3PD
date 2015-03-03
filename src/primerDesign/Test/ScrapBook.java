/**
 * 
 */
package primerDesign.Test;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.dsc.RestrictionSite;


/**
 * @author froehler
 *
 */
public class ScrapBook {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws IllegalSymbolException 
	 * @throws IllegalAlphabetException 
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception{
//		int num=100;
//		SimpleTimer timer = new SimpleTimer();
//		System.out.println("Doing " + num + " queries");
//		Runtime runtime = Runtime.getRuntime();
//		for(int i=0; i<num; i++){
//			//doQuery();
//			runtime.exec(" /Users/froehler/Downloads/BLAT-bin/gfClient -minIdentity=20 -minScore=1 abt4-web 10110 /Users/froehler/Downloads/BLAT-bin /export-cdfs/Sebastian/PrimerDesign/TestPrimer.fa /dev/stdout");
//		}
//		System.out.println("Done in " + timer.getTimeString());
		
//		Vector<Primer> primers = new Vector<Primer>();
//		primers.add(new Primer("AACGCCATACAATGATAACAA" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("AGTACCTATACTACCCGTGAAACATATTCGAAAA" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("TTCCAAAAGAATACGAACGC" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//		
//		primers.add(new Primer("TCTCACTTTTCTCACTGAACT" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("TCAATTGTAAGGAGTACTCTGTGGGGAAA" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("TGTGATACAAATACAATAATACCGTG" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//		
//		primers.add(new Primer("TGTATGTCATTGGTCCCTAAT" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("CGTTCGACTCGTTTGAACACTCTTGAG" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("CCAACTTTTCCTCGTGATTTA" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//		
//		primers.add(new Primer("CAAAGTGGTGAAAGAATCAATC" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("CCGGTCATATGCGTGTTTTCTAAAGAGTTAC" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("TATTGACCCAATTATCGAAACAG" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//		
//		primers.add(new Primer("GGGGCAATAATCAATCATAGTT" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("CTCCCTCTATCGATCACTCATTTCGTGTTG" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("GTGAAGCCTTCAAATAAAGTAGAT" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//		
//		primers.add(new Primer("GGGAATATAGAGGGATGTCAT" , PrimerTypes.forwardPrimer, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("CAATGAAGAACGATGAAGAGGGAATCGTC" , PrimerTypes.hybridizationProbe, 0, 0, 0, new PrimerSearchParameters()));
//		primers.add(new Primer("ACGGTGGGAGTATTCTG" , PrimerTypes.reversePrimer, 0, 0, 0, new PrimerSearchParameters()));
//	
//		//System.out.println(primer.toString());
//		PrimerSearchParameters params = new PrimerSearchParameters();
//		params.setPRINT_DEBUG_LOG(true);
//		params.setTargetOrganism(TargetOrganisms.Ppacificus);
//		params.setEnzyme(new RestrictionEnzyme("EcoRI",DNATools.createDNA("gaattc"), 0,0));
//		params.setConfigFile(PrimerDesignWebProperites.parseConfigFile(new FileInputStream("/Users/froehler/Downloads/apache-tomcat-6.0.16/webapps/3CPrimerDesign/3CPrimerDesign-ServletConfig.xml")));
//		PrimerMisprimingCheck check = PrimerMisprimingCheckDeserializer.deserialize("/export-web/froehler/SEQUENCE_INDICES/SEQUENCES/Ppacificus.fa.SlimGenome.ser", params);
//		
//		for(Primer primer : primers){
//			System.out.println(primer.getSequence() + " has mirprimings: " + check.hasMisprimings(primer));
//		}
//		File file = new File(args[0]);
//		int tileSize = Integer.parseInt(args[1]);
//		
//		HashSet<String> set = new HashSet<String>();
//		SlimFastaParser parser = new SlimFastaParser(file);
//		SimpleContig contig;
//		while(parser.hasNextContig()){
//			contig = parser.parseNextContigIgnoreCase();
//			for(int i=0; i<contig.getSequenceLength(); i+=tileSize){
//				
//			}
//		}
//		SlimFastaParser parser = new SlimFastaParser(">Site1\nGAGAACGCCATACAATGATAACAACAACAGAAATATTCGCTTCATTTTTAATATTCGCGCTCTATTTTTCGAATATGT\nTTCACGGGTAGTATAGGTACTGAATTCTGGCGTTGTTCAATTGTGTCTTCAAATGGCATTCAATCTCATAGAAAAAAG\nTTTGTTCGCGTTCGTATTCTTTTGGAAATGCGCTGAAGAAATGAN\n>Site2\nTTTTCTATTGTACTTTTTCTTCCGAGACATTTCAGAGTGTTCTCACTTTTCTCACTGAACTTTCCCTTCTTTTCCCCA\nCAGAGTACTCCTTACAATTGAGAATTCCTTTTCTTCCCTCTCGTACAAGGTTGACTTCGCTCTTGCAAGGTGTCCCCT\nCTCCCTCCATCACGGTATTATTGTATTTGTATCACATAAAGATAN\n>Site3\nTTGTATGTCATTGGTCCCTAATCTTTCATGATCTACTCATCTCATTCTTTCCTCTTTCTCATTCTCCTCAAGAGTGTT\nCAAACGAGTCGAACGAAGAACGAATTCTCACCAAAAACGAGCGCGCTGTGATAATTTTCTCTCATAAATCACGAGGAA\nAAGTTGGAAGAAGTTCGACCTTGGTTTCTCCTCTTTCCTTCTTCN\n>Site4\nTTCGAAGAATCAACGATTCATCCTCAAAGTGGTGAAAGAATCAATCATTTATCGATAGTGAAAGCAGTAACTCTTTAG\nAAAACACGCATATGACCGGTCGAATTCACCTGATCAGGTGGAATGAGTTTTCGAGGAATGATTTTGTGTCGCCAGAAC\nTGTTTCGATAATTGGGTCAATATTTTTAGGCCGGCTGAAACTATN\n>Site5\nAATCTGTGAAGAGGAAAGGAGAAGAGGGGAGGGGCAATAATCAATCATAGTTCGTTGCCAACCACTCCCAACACGAAA\nTGAGTGATCGATAGAGGGAGAGAATTCAACAAGTGAGGGGACCAATACACTTTAGTGGTGCAAGATGTGGTTCATGGA\nATAATTTATCTATAGCTTATCTACTTTATTTGAAGGCTTCACGTN\n>Site6\nTAACCATCCATTATTCTTCGCTCTCGATCCGGGAATATAGAGGGATGTCATCGATTTGACGACGATTCCCTCTTCATC\nGTTCTTCATTGCTTCTGCGAAGAATTCCTCCACTTCCGATCTGCAATCGATTTATGACACGTTAACCAGAATACTCCC\nACCGTACCTACCGTATCCCTCTCCCTCTAATCAATGAAATGCCTN");
//		while(parser.hasNextContig()){
//			System.out.print(parser.parseNextContig().toFastaString());
//		}
//		System.out.println();
		
//		RestrictionEnzymeMapper mapper = new RestrictionEnzymeMapper();
//		int[] rssPositions;
//		String contigSeq = "TTTTCTATTGTACTTTTTCTTCCGAGACATTTCAGAGTGTTCTCACTTTTCTCACTGAACTTTCCCTTCTTTTCCCCACAGAGTACTCCTTACAATTGAGAATTCCTTTTCTTCCCTCTCGTACAAGGTTGACTTCGCTCTTGCAAGGTGTCCCCTCTCCCTCCATCACGGTATTATTGTATTTGTATCACATAAAGATA";
//		rssPositions = mapper.mapFeature(contigSeq, new RestrictionEnzyme("EcoRI", DNATools.createDNA("gaattc"), 0, 0));
//		int pos = contigSeq.indexOf("GAATTC");
//		
//		if(rssPositions[0] < 200/2 || contigSeq.length() - rssPositions[0] < 200/2){
//			throw new IllegalStateException();
//		}
		
//		System.out.println(Runtime.getRuntime().availableProcessors());
		
//		String primer = "TTCGCAGAAGCAATGAAGAACGATGAAGAGGGAAT";
//		String enzyme = "GAATTC";
//		char[] test = (primer + enzyme).toCharArray();
//		Arrays.copyOfRange(original, from, to);
//		Matcher pattern = Pattern.compile(".*" + enzyme + ".*").matcher(CharBuffer.wrap(test));
//
//		System.out.println(pattern.matches());
		
//		int numReads = Integer.parseInt(args[0]);
//		int minLength = Integer.parseInt(args[1]);
//		int maxLength = Integer.parseInt(args[2]);
////		byte qualValue = Byte.parseByte(args[3]);
//		
//		String currentName;
//		String currentSequence;
//		byte[] quality;
//		for(int i=0; i< numReads; i++){
//			currentName = "test " + i;
//			currentSequence = SeqTools.getRandomPrimerSequence(minLength, maxLength);
////			quality = new byte[currentSequence.length()];
////			Arrays.fill(quality, qualValue);
////			System.out.print((new SimpleContigWithQuality(currentName, currentSequence.toCharArray(), quality)).toFastaString());
//			System.out.print(new SimpleContigImpl(currentName, currentSequence.toCharArray()).toFastaString());
//		}
//		System.out.println();
		
//		ArrayList<String[]> pairs = new ArrayList<String[]>();
//		pairs.add(new String[]{"AACGCCATACAATGATAACAA", "AGTACCTATACTACCCGTGAAACATATTCGAAAA", "CTTCAGCGCATTTCCAAA"});
//		pairs.add(new String[]{"TCTCACTTTTCTCACTGAACT", "CAATTGTAAGGAGTACTCTGTGGGGAAAAGA", "GTGATACAAATACAATAATACCGTGA"});
//		pairs.add(new String[]{"TTGTATGTCATTGGTCCCTAA", "CTTCGTTCGACTCGTTTGAACACTCT", "TCGAACTTCTTCCAACTTTTC"});
//		pairs.add(new String[]{"AAGTGGTGAAAGAATCAATCAT", "CGGTCATATGCGTGTTTTCTAAAGAGTTACTG", "TGACCCAATTATCGAAACAG"});
//		pairs.add(new String[]{"GGGGAGGGGCAATAATC", "TCTATCGATCACTCATTTCGTGTTGGGAG", "GTGAAGCCTTCAAATAAAGTAGATAA"});
//		pairs.add(new String[]{"GGGAATATAGAGGGATGTCAT", "CAATGAAGAACGATGAAGAGGGAATCGTC", "GATACGGTAGGTACGGTG"});
//		
//		int i=1;
//		System.out.println("My pairs");
//		for(String[] element : pairs){
//			System.out.println("### pair " + i);
//			SequenceRegionAligner aligner = new SequenceRegionAligner();
//			SequenceRegionAlignment scores2 = aligner.alignSequenceRegions(element[0].toCharArray(), element[2].toCharArray(), 2, 4);
//			PrimerAlignmentScores scores = scores2.getGlobalAlignmentValues(0, element[0].length(), 0, element[2].length());
//			System.out.println("Up Down");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			
//			aligner = new SequenceRegionAligner();
//			scores2 = aligner.alignSequenceRegions(element[0].toCharArray(), element[1].toCharArray(), 2, 4);
//			scores = scores2.getGlobalAlignmentValues(0, element[0].length(), 0, element[1].length());
//			System.out.println("Up Probe");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			
//			aligner = new SequenceRegionAligner();
//			scores2 = aligner.alignSequenceRegions(element[2].toCharArray(), element[1].toCharArray(), 2, 4);
//			scores = scores2.getGlobalAlignmentValues(0, element[2].length(), 0, element[1].length());
//			System.out.println("Down Probe");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			i++;
//		}
//		
//		pairs = new ArrayList<String[]>();
//		pairs.add(new String[]{"CAACAACAGAAATATTCGCT", "AGTACCTATACTACCCGTGAAACATATTCGAAAA", "CGCATTTCCAAAAGAATAC"});
//		pairs.add(new String[]{"TTCTCACTTTTCTCACTGAACT", "CAATTGTAAGGAGTACTCTGTGGGGAAAAGA", "CAAATACAATAATACCGTGATG"});
//		pairs.add(new String[]{"TGTCATTGGTCCCTAATCT", "CTTCGTTCGACTCGTTTGAACACTCT", "AAAGAGGAGAAACCAAGGT"});
//		pairs.add(new String[]{"TGGTGAAAGAATCAATCATT", "CGGTCATATGCGTGTTTTCTAAAGAGTTACTG", "TTGACCCAATTATCGAAAC"});
//		pairs.add(new String[]{"TCTGTGAAGAGGAAAGGAG", "TCTATCGATCACTCATTTCGTGTTGGGAG", "GATAAATTATTCCATGAACCAC"});
//		pairs.add(new String[]{"TTATTCTTCGCTCTCGATC", "CAATGAAGAACGATGAAGAGGGAATCGTC", "TGATTAGAGGGAGAGGGAT"});
//		
//		i=1;
//		System.out.println("MuPleX pairs");
//		for(String[] element : pairs){
//			System.out.println("### pair " + i);
//			SequenceRegionAligner aligner = new SequenceRegionAligner();
//			SequenceRegionAlignment scores2 = aligner.alignSequenceRegions(element[0].toCharArray(), element[2].toCharArray(), 2, 4);
//			PrimerAlignmentScores scores = scores2.getGlobalAlignmentValues(0, element[0].length(), 0, element[2].length());
//			System.out.println("Up Down");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			
//			aligner = new SequenceRegionAligner();
//			scores2 = aligner.alignSequenceRegions(element[0].toCharArray(), element[1].toCharArray(), 2, 4);
//			scores = scores2.getGlobalAlignmentValues(0, element[0].length(), 0, element[1].length());
//			System.out.println("Up Probe");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			
//			aligner = new SequenceRegionAligner();
//			scores2 = aligner.alignSequenceRegions(element[2].toCharArray(), element[1].toCharArray(), 2, 4);
//			scores = scores2.getGlobalAlignmentValues(0, element[2].length(), 0, element[1].length());
//			System.out.println("Down Probe");
//			System.out.println("PA: " + scores.getPairScore() + " PEA: " + scores.getPairEndScore());
//			i++;
		
		
//		MultiSeqESAFatOpt esa = MultiSeqESAFatOpt.deserialize(new File("/Users/froehler/SEQUENCE_INDICES/Ppacificus-unmasked.fa.MultiSeqESAFatOpt.esaidx"));
//		System.out.println("Hits: " + (esa.findHitCount("TCCAAC") + esa.findHitCount("TCCGAC")));
//		System.out.println("Hits: " + (esa.findHitCount("GAATTC")));
		
		//System.out.println(SeqTools.revcompDNA("CGTTGGAttggTTTggttTCCAACGAATT".toCharArray()));
		
//		String fwPrimer = "TTATTCTTCGCTCTCGATC";
//		String revPrimer = "TCTGTGAAGAGGAAAGGAG";
////		
//		SequenceRegionAligner aligner = new SequenceRegionAligner();
//		SequenceRegionAlignment alignment = aligner.alignSequenceRegions(fwPrimer.toCharArray(), revPrimer.toCharArray(), 2, 4);
//		System.out.println(alignment.getGlobalAlignmentValues(0, fwPrimer.length(), 0, revPrimer.length()).getPairScore());
//		System.out.println(alignment.getGlobalAlignmentValues(0, fwPrimer.length(), 0, revPrimer.length()).getPairEndScore());
//		PrimerSearchParameters params = new PrimerSearchParameters();
//		params.setPickTaqManProbe(false);
//		
//		PrimerPair pair = new PrimerPair(new Primer("TTATTCTTCGCTCTCGATC", PrimerTypes.forwardPrimer, 89, 20, 12, params),
//				new Primer("TGATTAGAGGGAGAGGGAT", PrimerTypes.reversePrimer, 84, 8, 4, params),
//				new PrimerAlignmentScores(30,10), params);
//		System.out.println(pair.getDistanceToOptimalPrimerPair());

//		SantaLuciaTM tm = new SantaLuciaTM();
//		PrimerSearchParameters params = new PrimerSearchParameters();
//		
//		ArrayList<String> list = new ArrayList<String>();
//		list.add("TTATTCTTCGCTCTCGATC");
//		list.add("TGATTAGAGGGAGAGGGAT");
//		list.add("GCGCG");
//		list.add("ATATA");
//		
//		for(String primer : list){
//			System.out.println(primer.substring(primer.length() - 5) + ": " + tm.computeDeltaG(primer.substring(primer.length() - 5), 50E-09, 50E-03, 0, 0));
//		}

		
//		SlimFastaParser parser = new SlimFastaParser(new File("/Users/froehler/SEQUENCE_INDICES/Hsapiens/Hsapiens.fa"));
////		System.out.println(parser.parseNextContigIgnoreCase().toFastaString());
//		System.out.println(parser.parseNextContigIgnoreCase().toFastaString());
		
		RestrictionEnzyme enzyme = new RestrictionEnzyme("StyI", DNATools.createDNA("CCWWGG"), 0, 0);
		System.out.println(enzyme);
		System.out.println(enzyme.getForwardRegex());
		
		
	}
	
	public static void doQuery() throws IOException{
        URL url = new  URL("http://www.pristionchus.org/cgi-bin/webBlat");
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("POST");
        conn.setDoInput(true);
        conn.setDoOutput(true);
        conn.setRequestProperty ("Content-Type","application/x-www-form-urlencoded");

        DataOutputStream out = new DataOutputStream(conn.getOutputStream());
        String content = "wb_db=" + URLEncoder.encode ("P.pacificus genome")
        + "&wb_qType=" + URLEncoder.encode ("DNA")
        + "&wb_sort=" + URLEncoder.encode("query,score")
        + "&wb_output=" + URLEncoder.encode("psl no header")
        + "&wb_seq=" + URLEncoder.encode(">Test\nCAATGAAGAACGATGAAGAGGGAATCGTC");
        out.writeBytes (content);
		out.flush ();
		out.close ();

        conn.connect();
        
        InputStream stream = conn.getInputStream();
        BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
        String line;
	    while((line = reader.readLine()) != null){
	    	//System.out.println(line);
	    }
	}
	
	public static int getIndex(RestrictionSite site, int code){
		int pos = -1;
		for(int i=0; i<site.getPrimerPairs().size(); i++){
			if(code == site.getPrimerPair(i).getHashCode()){
				pos = i;
				break;
			}
		}
		return pos;
	}
}
