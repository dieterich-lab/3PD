package primerDesign.web;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Properties;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.servlet.ServletException;
import javax.servlet.UnavailableException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.seq.DNATools;

import primerDesign.algo.PrimerSearch;
import primerDesign.algo.SimpleGreedyPrimerPairPicking;
import primerDesign.dsc.PrimerPairSet;
import primerDesign.dsc.RestrictionSite;
import primerDesign.dsc.indexStructures.TargetOrganisms;
import primerDesign.dsc.indexStructures.blat.BLATQueryClient;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.Blat3CPrimerMisprimingScanPlusSeq;
import primerDesign.dsc.indexStructures.primerMisprimingCheck.PrimerMisprimingCheck;
import primerDesign.util.Constants;
import primerDesign.util.EmptyResultSetException;
import primerDesign.util.MyExtendedMath;
import primerDesign.util.NoRestrictionSitesFoundException;
import primerDesign.util.PrimerSearchParameters;
import primerDesign.util.RestrictionEnzymeMapper;
import primerDesign.util.SeqTools;
import primerDesign.util.SimpleContigImpl;
import primerDesign.util.SlimFastaParser;
import primerDesign.web.tools.EMailVerifier;
import primerDesign.web.tools.FeaturePositionPicDrawer;
import primerDesign.web.tools.Mail;
import primerDesign.web.tools.OptPrimerSetPrinter;
import primerDesign.web.tools.PrintParameters;
import cern.colt.list.ObjectArrayList;

/**
 * Implements a servlet for processing 3C primer design tasks coming from a webinterface.
 * 
 * For memory reasons, each job is put into a queue and executed serially. 
 * (This was a pre-requisite while using the ESA index and is kept until now despite a serial execution is not critical to this servlet anymore).
 * 
 * @author Sebastian Fr�hler
 *
 */
public class CCCPrimerDesign extends HttpServlet{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Properties configFile;
	private static final Pattern illegalFastaChars = Pattern.compile(".*[^AaTtGgCcUuMmRrWwSsYyKkVvHhDdBbNn]+.*");
	private static int jobID = 0;
	
	private static ObjectArrayList queue;
	private static QueueProcessor processor;
	private CCCEmailProperties mailProperties;
	private static final String pictureFilename = "/tmp/CCCPrimerDesign-results.png";
	
	/**
	 * Initializes the servlet by reading in the config file, reading in the RSS index and setting up the queue.
	 */
	public void init() throws UnavailableException{
		// read in config file with sequence index paths
		try {
			String configFileName = "/var/lib/tomcat7/webapps/3CPrimerDesign/3CPrimerDesign-ServletConfig.xml"; 
			this.configFile = PrimerDesignWebProperites.parseConfigFile(new FileInputStream(configFileName));
		} catch (IOException e) {
			e.printStackTrace();
			throw new UnavailableException("Error processing config file for this servlet!");
		}
		
		if(configFile.getProperty("MAX_NUM_PICKING_THREADS").equals("all")){
			Constants.MAX_NUM_PICKING_THREADS = Runtime.getRuntime().availableProcessors();
		}
		else{
			if(Integer.parseInt(configFile.getProperty("MAX_NUM_PICKING_THREADS")) < 1 || Integer.parseInt(configFile.getProperty("MAX_NUM_PICKING_THREADS")) > Runtime.getRuntime().availableProcessors()) throw new IllegalArgumentException("The number of picking threads on this system must be: 1<=x<=" + Runtime.getRuntime().availableProcessors() +"!");
			Constants.MAX_NUM_PICKING_THREADS = Integer.parseInt(configFile.getProperty("MAX_NUM_PICKING_THREADS"));
		}
		
		this.mailProperties = new CCCEmailProperties(configFile.getProperty("THIS_SERVICE_SENDER_EMAIL"), configFile.getProperty("ERROR_REPORT_EMAIL"),configFile.getProperty("SMTP_SERVER"));
		
		//this.rssIndex = RestrictionSitesIndex.deserialize(this.configFile.getProperty("BLAT_RSS_INDEX"));
		
		CCCPrimerDesign.queue = new ObjectArrayList();
		CCCPrimerDesign.processor = new QueueProcessor();
		
		getServletContext().log("Started 3C primer design server");
		processor.start();
	}

	/**
	 * Processes a 3C primer design request, creates a job and puts this job into the queue.
	 * 
	 * @param request the http request from the servelet engine
	 * @param response the http response to return to the servelet engine
	 */
	public void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException{
		PrintWriter out = response.getWriter();
		try{
			response.setContentType("text/html");
			
			int id = jobID++;
			
			// init and set search parameters
			PrimerSearchParameters searchParams = new PrimerSearchParameters();
			
			// check for errors, display error message?
			String value = request.getParameter("TargetOrganism");
			String targetOrganism;
			
			if(value.equals("C.elegans")){
				targetOrganism = TargetOrganisms.Celegans.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Celegans);
			}
			else if(value.equals("P.pacificus")){
				targetOrganism = TargetOrganisms.Ppacificus.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Ppacificus);
			}
			else if(value.equals("M.musculus")){
				targetOrganism = TargetOrganisms.Mmusculus.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Mmusculus);
			}
			else if(value.equals("D.melanogaster")){
				targetOrganism = TargetOrganisms.Dmelanogaster.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Dmelanogaster);
			}
			else if(value.equals("H.sapiens")){
				targetOrganism = TargetOrganisms.Hsapiens.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Hsapiens);
			}
			else if(value.equals("S.cerevisiae")){
				targetOrganism = TargetOrganisms.Scerevisiae.toString();
				searchParams.setTargetOrganism(TargetOrganisms.Scerevisiae);
			}
			else throw new IllegalArgumentException("Unsupported background index: " + value);
			BLATQueryClient client = new BLATQueryClient(this.configFile.getProperty("BLAT_gfClient"), this.configFile.getProperty("BLAT_host"), this.configFile.getProperty("BLAT_PORT_" + targetOrganism), this.configFile.getProperty("BLAT_SequenceDir"));
			PrimerMisprimingCheck misprimingScan = new Blat3CPrimerMisprimingScanPlusSeq(client, this.configFile.getProperty("SEQUENCE_PATH") + "/" + this.configFile.getProperty(targetOrganism + "-SequenceFile"), searchParams);
	
			if(!client.isAvailable()){
				out.println(getErrorMessage("This service is temporarily unavailable for the selected target organism!"));
				out.close();
				Mail.sendMail(this.mailProperties.getSender(), this.mailProperties.getSmtp(), this.mailProperties.getAdmin(), "3CPrimerDesign - job ERROR OCCURED", "Index for organism: " + targetOrganism + " is not available!");
				return;
			}
			
			// process restriction enzyme
			String[] values = request.getParameter("restrictionEnzyme").split(":");
			if(values.length != 4){
				throw new IllegalArgumentException("Unsupported restriction enzyme: " + request.getParameter("restrictionEnzyme"));
			}
			String enzymeName;
			String enzymeSite;
			int enzymeFwCutPos;
			int enzymeRevCutPos;
			
			try{
				enzymeName = values[0];
				enzymeSite = values[1];
				enzymeFwCutPos = Integer.parseInt(values[2]);
				enzymeRevCutPos = Integer.parseInt(values[3]);
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing restriction enzyme: " + request.getParameter("restrictionEnzyme"));
			}
			
			int nbPrimers = Integer.parseInt(request.getParameter("numberOfPrimerPairs"));
			if(nbPrimers <= 1) throw new IllegalArgumentException("The number of primer pairs must be >= 2! : " + nbPrimers);

			if(!request.getParameter("pickTaqManProbes").equals("true") && !request.getParameter("pickTaqManProbes").equals("false")){
				throw new IllegalArgumentException("Error parsing pickTaqManProbes - not boolean!");
			}
			searchParams.setPickTaqManProbe(Boolean.parseBoolean(request.getParameter("pickTaqManProbes")));
			
			try
			{
				searchParams.setMIN_RESTRICTION_FRAGMENT_LENGTH(Integer.parseInt(request.getParameter("restrictionFragmentMinSize")));
				searchParams.setMAX_RESTRICTION_FRAGMENT_LENGTH(Integer.parseInt(request.getParameter("restrictionFragmentMaxSize")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/max restriction fragment lengths: " + request.getParameter("restrictionFragmentMinSize") + "/" + request.getParameter("restrictionFragmentMaxSize"));
			}

			if(searchParams.getMIN_RESTRICTION_FRAGMENT_LENGTH() < 0 || searchParams.getMAX_RESTRICTION_FRAGMENT_LENGTH() < 0 || searchParams.getMAX_RESTRICTION_FRAGMENT_LENGTH() < searchParams.getMIN_RESTRICTION_FRAGMENT_LENGTH()){
				out.println(getErrorMessage("Error processing fragment length parameters: min restriction fragment length has to be <= max restriction fragment length!"));
				out.close();
				throw new IllegalArgumentException("Error processing fragment length parameters: min: '" + request.getParameter("restrictionFragmentMinSize") + "' max: '" + request.getParameter("restrictionFragmentMaxSize") + "'");
			}

			// primer parameters
			try{
				searchParams.setMIN_PRIMER_LENGTH(Integer.parseInt(request.getParameter("minPrimerLength")));
				searchParams.setOPT_PRIMER_LENGTH(Integer.parseInt(request.getParameter("optPrimerLength")));
				searchParams.setMAX_PRIMER_LENGTH(Integer.parseInt(request.getParameter("maxPrimerLength")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max primer lengths: " + request.getParameter("minPrimerLength") + "/" + request.getParameter("optPrimerLength") + "/" + request.getParameter("maxPrimerLength"));
			}
			
			if(searchParams.getMIN_PRIMER_LENGTH() > searchParams.getOPT_PRIMER_LENGTH() ||
			   searchParams.getOPT_PRIMER_LENGTH() > searchParams.getMAX_PRIMER_LENGTH()){
				out.println(getErrorMessage("Contraint violated: Min primer length <= Opt primer length <= Max primer length!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setMIN_TM(Double.parseDouble(request.getParameter("minPrimerTM")));
				searchParams.setOPT_TM(Double.parseDouble(request.getParameter("optPrimerTM")));
				searchParams.setMAX_TM(Double.parseDouble(request.getParameter("maxPrimerTM")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max primer TM: " + request.getParameter("minPrimerTM") + "/" + request.getParameter("optPrimerTM") + "/" + request.getParameter("minPrimerTM"));
			}
			
			if(searchParams.getMIN_TM() > searchParams.getOPT_TM() ||
			   searchParams.getOPT_TM() > searchParams.getMAX_TM()){
				out.println(getErrorMessage("Contraint violated: Min primer Tm <= Opt primer Tm <= Max primer Tm!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setMIN_GC(Double.parseDouble(request.getParameter("minPrimerGC")));
				searchParams.setOPT_GC(Double.parseDouble(request.getParameter("optPrimerGC")));
				searchParams.setMAX_GC(Double.parseDouble(request.getParameter("maxPrimerGC")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max primer GC: " + request.getParameter("minPrimerGC") + "/" + request.getParameter("optPrimerGC") + "/" + request.getParameter("maxPrimerGC"));
			}
			
			if(searchParams.getMIN_GC() > searchParams.getOPT_GC() ||
				searchParams.getOPT_GC() > searchParams.getMAX_GC()){
				out.println(getErrorMessage("Contraint violated: Min primer %GC <= Opt primer %GC <= Max primer %GC!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setMAX_PRIMER_TM_DIFFERENCE(Double.parseDouble(request.getParameter("maxPrimerTMDifference")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing max primer TM difference: " + request.getParameter("maxPrimerTMDifference"));
			}
			
			// probe parameters
			try{
				searchParams.setTAQMAN_MIN_PRIMER_LENGTH(Integer.parseInt(request.getParameter("minProbeLength")));
				searchParams.setTAQMAN_OPT_PRIMER_LENGTH(Integer.parseInt(request.getParameter("optProbeLength")));
				searchParams.setTAQMAN_MAX_PRIMER_LENGTH(Integer.parseInt(request.getParameter("maxProbeLength")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max probe lengths: " + request.getParameter("minProbeLength") + "/" + request.getParameter("optProbeLength") + "/" + request.getParameter("maxProbeLength"));
			}
			
			if(searchParams.getTAQMAN_MIN_PRIMER_LENGTH() > searchParams.getTAQMAN_OPT_PRIMER_LENGTH() ||
				searchParams.getTAQMAN_OPT_PRIMER_LENGTH() > searchParams.getTAQMAN_MAX_PRIMER_LENGTH()){
				out.println(getErrorMessage("Contraint violated: Min probe length <= Opt probe length <= Max probe length!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setTAQMAN_MIN_TM(Double.parseDouble(request.getParameter("minProbeTM")));
				searchParams.setTAQMAN_OPT_TM(Double.parseDouble(request.getParameter("optProbeTM")));
				searchParams.setTAQMAN_MAX_TM(Double.parseDouble(request.getParameter("maxProbeTM")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max probe TM: " + request.getParameter("minProbeTM") + "/" + request.getParameter("optProbeTM") + "/" + request.getParameter("maxProbeTM"));
			}
			
			if(searchParams.getTAQMAN_MIN_TM() > searchParams.getTAQMAN_OPT_TM() ||
				searchParams.getTAQMAN_OPT_TM() > searchParams.getTAQMAN_MAX_TM()){
				out.println(getErrorMessage("Contraint violated: Min probe Tm <= Opt probe Tm <= Max probe Tm!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setTAQMAN_MIN_GC(Double.parseDouble(request.getParameter("minProbeGC")));
				searchParams.setTAQMAN_OPT_GC(Double.parseDouble(request.getParameter("optProbeGC")));
				searchParams.setTAQMAN_MAX_GC(Double.parseDouble(request.getParameter("maxProbeGC")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max probe GC: " + request.getParameter("minProbeGC") + "/" + request.getParameter("optProbeGC") + "/" + request.getParameter("maxProbeGC"));
			}
			
			if(searchParams.getTAQMAN_MIN_GC() > searchParams.getTAQMAN_OPT_GC() ||
				searchParams.getTAQMAN_OPT_GC() > searchParams.getTAQMAN_MAX_GC()){
				out.println(getErrorMessage("Contraint violated: Min probe %GC <= Opt probe %GC <= Max probe %GC!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setMIN_TAQMAN_TM_DIFFERENCE(Double.parseDouble(request.getParameter("minProbeTMDifference")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min probe TM difference: " + request.getParameter("minProbeTMDifference"));
			}
			
			// amplicon properties
			try{
				searchParams.setMIN_AMPLICON_LENGTH(Integer.parseInt(request.getParameter("minAmpliconLength")));
				searchParams.setOPT_AMPLICON_LENGTH(Integer.parseInt(request.getParameter("optAmpliconLength")));
				searchParams.setMAX_AMPLICON_LENGTH(Integer.parseInt(request.getParameter("maxAmpliconLength")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing min/opt/max amplicon lengths: " + request.getParameter("minAmpliconLength") + "/" + request.getParameter("optAmpliconLength") + "/" + request.getParameter("maxAmpliconLength"));
			}
			
			if(searchParams.getMIN_AMPLICON_LENGTH() > searchParams.getOPT_AMPLICON_LENGTH() ||
				searchParams.getOPT_AMPLICON_LENGTH() > searchParams.getMAX_AMPLICON_LENGTH()){
				out.println(getErrorMessage("Contraint violated: Min amplicon length <= Opt amplicon length <= Max amplicon length!"));
				out.close();
				return;
			}
			
			try{
				searchParams.setSAFE_FALSE_POSITIVE_AMPLICON_LENGTH(Integer.parseInt(request.getParameter("safeFPAmpliconLength")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing safe FP amplicon length: " + request.getParameter("safeFPAmpliconLength"));
			}
			
			// alignment scores
			try{
				searchParams.setMAX_PRIMER_SELF_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("saMax")));
				searchParams.setMAX_PRIMER_SELF_END_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("seaMax")));
				searchParams.setMAX_PRIMER_PAIR_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("paMax")));
				searchParams.setMAX_PRIMER_PAIR_END_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("peaMax")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing primer saMax/seaMax/paMax/peaMax:: " + request.getParameter("saMax") + "/" + request.getParameter("seaMax") + "/" + request.getParameter("paMax") + "/" + request.getParameter("peaMax"));
			}
			
			try{
				searchParams.setMAX_TAQMAN_SELF_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("probeSaMax")));
				searchParams.setMAX_TAQMAN_SELF_END_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("probeSeaMax")));
				searchParams.setMAX_TAQMAN_PAIR_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("probePaMax")));
				searchParams.setMAX_TAQMAN_PAIR_END_ALIGNMENT_SCORE(Integer.parseInt(request.getParameter("probePeaMax")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing probe saMax/seaMax/paMax/peaMax: " + request.getParameter("probeSaMax") + "/" + request.getParameter("probeSeaMax") + "/" + request.getParameter("probePaMax") + "/" + request.getParameter("probePeaMax"));
			}
			
			// weights
			try{
				searchParams.setPRIMER_DELTA_TM_WEIGHT(Double.parseDouble(request.getParameter("deltaTMWeight")));
				searchParams.setPRIMER_DELTA_GC_WEIGHT(Double.parseDouble(request.getParameter("deltaGCWeight")));
				searchParams.setPRIMER_DELTA_LENGTH_WEIGHT(Double.parseDouble(request.getParameter("deltaLengthWeight")));
				searchParams.setPRIMER_DELTA_DISTANCE_TO_RSS_WEIGHT(Double.parseDouble(request.getParameter("deltaDRSSWeight")));
				searchParams.setSELF_ALIGNMENT_WEIGHT(Double.parseDouble(request.getParameter("deltaSAWeight")));
				searchParams.setSELF_END_ALIGNMENT_WEIGHT(Double.parseDouble(request.getParameter("deltaSEAWeight")));
				searchParams.setPAIR_ALIGNMENT_WEIGHT(Double.parseDouble(request.getParameter("deltaPAWeight")));
				searchParams.setPAIR_END_ALIGNMENT_WEIGHT(Double.parseDouble(request.getParameter("deltaPEAWeight")));
				searchParams.setPRIMER_FALSE_POSITIVES_WEIGHT(Double.parseDouble(request.getParameter("FPWeight")));
			}
			catch(NumberFormatException e){
				throw new IllegalArgumentException("Error parsing weights dTM/dGC/dL/dDRSS/dSA/dSEA/dPA/dPEA/FP: " + request.getParameter("deltaTMWeight") + "/" + request.getParameter("deltaGCWeight") + "/" + request.getParameter("deltaLengthWeight") + "/" + request.getParameter("deltaDRSSWeight") + "/" + request.getParameter("deltaSAWeight") + "/" + request.getParameter("deltaSEAWeight") + "/" + request.getParameter("deltaPAWeight") + "/" + request.getParameter("deltaPEAWeight") + "/" + request.getParameter("FPWeight"));
			}
			
			if(!EMailVerifier.isValidGlobalEmailAddress(request.getParameter("e-mail"))){
				out.println(getErrorMessage("Invalid E-Mail adress - please check syntax!"));
				out.close();
				return;
			}
			
			// process sequence(s)
			SlimFastaParser parser = new SlimFastaParser(request.getParameter("sequence"));
			Vector<SimpleContigImpl> contigs = new Vector<SimpleContigImpl>();
			while(parser.hasNextContig()){
				try{
					contigs.add(parser.parseNextContigIgnoreCase());
				}
				catch(IOException e){
					out.println(getErrorMessage("Error reading FASTA DNA sequence(s) - please check sequence(s)!"));
					out.close();
					return;
				}
			}
			
			// threaded queue execution
			RestrictionEnzyme enzyme = null;
			try{
				enzyme = new RestrictionEnzyme(enzymeName, DNATools.createDNA(enzymeSite), enzymeFwCutPos, enzymeRevCutPos);
			}catch(Exception e){
				e.printStackTrace();
				getServletContext().log("Error processing restriction enzyme: '" + enzymeName + "' '" + enzymeSite + "' '" + enzymeFwCutPos + "' '" + enzymeRevCutPos);
			}
			//searchParams.setSequence(sequence);
			searchParams.setContigs(contigs.toArray(new SimpleContigImpl[contigs.size()]));
			searchParams.setEnzyme(enzyme);
			// depending on operation mode, auto-set number of primer pairs to pick
			if(request.getParameter("Mode").equals("TargetedMode")) searchParams.setNumPrimers(contigs.size());
			else searchParams.setNumPrimers(nbPrimers);
			searchParams.setPrimerMisprimingCheck(misprimingScan);
			
			// check integrity of mail adress
	//		if(!EMailVerifier.isValidMailAdress(request.getParameter("e-mail"))){
	//			out.println(getErrorMessage("Errors in e-mail adress format: " + request.getParameter("e-mail")));
	//			out.close();
	//			return;
	//		}
			
			RestrictionEnzymeMapper mapper = new RestrictionEnzymeMapper();
			int[] rssPositions;
			String contigSeq;
			
			// #Contigs != #Pairs check
			if(contigs.size() > 1 && searchParams.getNumPrimers() != contigs.size()){
				out.println(getErrorMessage("Contraint violated: The number of primer pairs to pick must equal the number of sequences supplied in targeted search mode!"));
				out.close();
				return;
			}
			// sequence length checks
			for(int i=0; i<contigs.size(); i++){
				contigSeq = new String(contigs.get(i).getSequence());
				
				// scan region total length check
				if(contigs.get(i).getSequenceLength() < searchParams.getMAX_AMPLICON_LENGTH()){
					out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") must have length > 'max amplicon length'"));
					out.close();
					return;
				}
				// fasta sequence content check
				if(illegalFastaChars.matcher(contigSeq).matches()){
					out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") contains invalid characters!"));
					out.close();
					return;
				}
				// primers scan regions check - require at least most frequent RSS(s) to have sufficient sequence to screen for primers
				switch((rssPositions = mapper.mapFeature(contigSeq, searchParams.getEnzyme())).length){
					// no RSS at all
					case 0 : {					
						out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") does not contain a restriction site for enzyme " + enzymeName));
						out.close();
						return;
					}
					// two fragments (=one RSS)
					case 1 : {
						if(rssPositions[0] < searchParams.getMAX_AMPLICON_LENGTH()/2 || contigSeq.length() - rssPositions[0] < searchParams.getMAX_AMPLICON_LENGTH()/2){
							out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") must have sequence lengths > 'max amplicon length/2' at either side of the restriction site!"));
							out.close();
							return;
						}
						else break;
					}
					// three and more fragments (>= 2 RSSs)
					default : {
						if(rssPositions.length % 2 != 0){
							// if odd number of restriction sites, check most centred RSS only							
							if(rssPositions[Math.round(rssPositions.length/2)] < searchParams.getMAX_AMPLICON_LENGTH()/2 || contigSeq.length() - rssPositions[Math.round(rssPositions.length/2)] <= searchParams.getMAX_AMPLICON_LENGTH()/2){
								out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") must have sequence lengths > 'max amplicon length/2' at either side of the restriction site!"));
								out.close();
								return;
							}
						}
						else{
							// if even number of restriction sites, check bost most-centred RSSs
							if(Math.min(rssPositions[MyExtendedMath.floor(rssPositions.length/2)], contigSeq.length() - rssPositions[MyExtendedMath.floor(rssPositions.length/2)]) <= searchParams.getMAX_AMPLICON_LENGTH()/2 && Math.min(rssPositions[MyExtendedMath.ceil(rssPositions.length/2)], contigSeq.length() - rssPositions[MyExtendedMath.ceil(rssPositions.length/2)]) <= searchParams.getMAX_AMPLICON_LENGTH()/2){
								out.println(getErrorMessage("Contraint violated: sequence " + (i+1) + " (" + contigs.get(i).getID() + ") must have sequence lengths > 'max amplicon length/2' at either side of the restriction site!"));
								out.close();
								return;
							}
						}	
					}
				}
			}

			out.println("<html><head><title>Primer design for chromatin conformation capture assays</title></head>" +
					"<body><h1>Job submitted!</h1>\n" + 
					"Your job \'" + request.getParameter("jobID") + "\' is beeing processed!\n" + 
					"Results will be mailed to your e-mail adress: \'" + request.getParameter("e-mail") + "\' a.s.a.p! "); 
					out.println("</body></html>");
			out.close();
			
			getServletContext().log("Receiving 3C primer design request: " + id + " from: " + request.getRemoteAddr() + " for user: " + request.getParameter("e-mail"));
			getServletContext().log("PARAMS:\t" + id + "\t" + request.getParameter("jobID") + "\t" + request.getParameter("e-mail") + "\t" + searchParams.getLogString());
			
			synchronized (CCCPrimerDesign.queue) {
				CCCPrimerDesign.queue.add(new QueueElement(id, request.getParameter("jobID"), request.getParameter("e-mail"), this.mailProperties, searchParams));
				CCCPrimerDesign.queue.notify();
				Mail.sendMail(this.mailProperties.getSender(), this.mailProperties.getSmtp(), this.mailProperties.getAdmin(), "3CPrimerDesign - job " + id + " queued for user: " + request.getParameter("e-mail") + "\n\n", "");
			}
		}catch(NullPointerException e){
			e.printStackTrace();
			Mail.sendMail(this.mailProperties.getSender(), this.mailProperties.getSmtp(), this.mailProperties.getAdmin(), "3CPrimerDesign - job ERROR OCCURED", "Error in primer design job: " + (jobID - 1)  + " for user: " + request.getParameter("e-mail") + "\n\n" + returnStackTraceString(e));
			out.println(getErrorMessage("There was an error processing your job - this incident has been reported!"));
			out.close();
		}catch(Throwable e){
			Mail.sendMail(this.mailProperties.getSender(), this.mailProperties.getSmtp(), request.getParameter("e-mail"), "3CPrimerDesign - job ERROR OCCURED", "An error occured during primer design - this has been reported!");
			Mail.sendMail(this.mailProperties.getSender(), this.mailProperties.getSmtp(), this.mailProperties.getAdmin(), "3CPrimerDesign - job ERROR OCCURED", "Error in primer design job: " + (jobID - 1)  + " for user: " + request.getParameter("e-mail") + "\n\n" + returnStackTraceString(e));
			e.printStackTrace();
			out.println(getErrorMessage("There was an error processing your job - this incident has been reported!"));
			out.close();
		}
	}
	
	/**
	 * Processes a 3PD GET request returning the status of the servelet.
	 * 
	 * @param request the http servelet request
	 * @param response the http servelet response
	 */
	public void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException{
		if(request.getParameter("PrintStatus") != null && request.getParameter("PrintStatus").equals("true")){
			response.setContentType("text/html");
			PrintWriter out = response.getWriter();
			NumberFormat format = NumberFormat.getInstance();
			out.println("<html><head><title>Primer design for chromatin conformation capture assays - Status page</title></head>" +
					"<body><h1>Status:</h1>\n" + 
					"Job count: " + format.format(CCCPrimerDesign.jobID) + " including " + format.format(CCCPrimerDesign.queue.size()) + " jobs in the queue! </body></html>");
			out.close();
		}
		else doPost(request, response);
	}
	
	/**
	 * Returns the stack trace string of a throwable - for debugging and diagnostic purposes.
	 * 
	 * @param throwable the throwable
	 * 
	 * @return the associated stack trace
	 */
	private static String returnStackTraceString(Throwable throwable){
		StringBuffer buffy = new StringBuffer();
		buffy.append(throwable.getCause() + ": " + throwable.getLocalizedMessage());
		buffy.append("\n");
		
		StackTraceElement[] array = throwable.getStackTrace();
		for(StackTraceElement item : array){
			buffy.append(item.toString());
			buffy.append("\n");
		}
		return buffy.toString();
	}
	
	/**
	 * Returns an error message in case of an error the user can fix by changing some parameters of his request.
	 * 
	 * @param message the message pointing the user to his error
	 * 
	 * @return an error message in case of an error the user can fix by changing some parameters of his request
	 */
	public static String getErrorMessage(String message){
		return new String("<html><head><title>Primer design for chromatin conformation capture assays</title></head>" +
				"<body><h1>Error submitting job!</h1>\n" + 
				message + "<br><br>" + 
				"Press the return button in your browser to go back..." +
				"</body></html>");
	}
	
	/**
	 * A simple queue thread managing the queue of 3C primer design requests.
	 * 
	 * Currently - due to memory constraints - requests need to be processed serially.
	 * This thread always processes the 'oldest' request in the queue first (FIFO).
	 * 
	 * @author Sebastian Fr�hler
	 *
	 */
	private class QueueProcessor extends Thread{
		public void run(){
			while(true){
				while(!CCCPrimerDesign.queue.isEmpty()){
					// pop the 'oldest' element from the queue and delete it from the queue
					Runnable currentThread = (Runnable) CCCPrimerDesign.queue.get(0);
					synchronized (CCCPrimerDesign.queue) {
						CCCPrimerDesign.queue.remove(0);
					}
					// run current thread
					currentThread.run();
				}
				try {
					synchronized (CCCPrimerDesign.queue) {
						CCCPrimerDesign.queue.wait();
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Encapsulates a 3C primer design job to be put into the queue.
	 * 
	 * @author Sebastian Fr�hler
	 *
	 */
	private class QueueElement extends Thread{
		private int jobID;
		private String eMail;
		private CCCEmailProperties emailProperties;
		private String userJobID;
		private PrimerSearchParameters searchParams;
		
		/**
		 * Create a ne 3C primer design job.
		 * 
		 * @param jobID the internal job id of this service (the same as in the logfile!) 
		 * @param userJobID the job id as specified by the user
		 * @param eMail the mail adress to send the results to as specified by the user
		 * @param searchParams the search parameters of this 3C primer design job
		 */
		public QueueElement(int jobID, String userJobID, String eMail, CCCEmailProperties emailProperties, PrimerSearchParameters searchParams){
			this.jobID = jobID;
			this.userJobID = userJobID;
			this.eMail = eMail;
			this.emailProperties = emailProperties;
			this.searchParams = searchParams;
		}
		
		/**
		 * Executes the job and sends some result/error message to the user-specified e-mail adress.
		 */
		public void run(){
			try{
				// submit job
				Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.emailProperties.getAdmin(), "3CPrimerDesign - job " + jobID + " started for user: " + this.eMail + "\n\n", "");
				getServletContext().log("Starting job " + this.jobID + " picking oldest element from queue of size " + (CCCPrimerDesign.queue.size() + 1));
				
				String body = "";
				RestrictionSite[] optimalSites = null;
				PrimerPairSet bestPrimerPairSet = null;
				
				OptPrimerSetPrinter optSetOut = new OptPrimerSetPrinter();
				PrintParameters paramPrinter = new PrintParameters();
				try{	
					PrimerSearch search = new PrimerSearch();

					search.setEnzymePatterns(Pattern.compile(".*" + searchParams.getEnzyme().getRecognitionSite().seqString().toUpperCase() + ".*"), Pattern.compile(".*" + SeqTools.revcompDNA(searchParams.getEnzyme().getRecognitionSite().seqString().toUpperCase().toCharArray()) + ".*"));

					searchParams.setPrimerSearch(search);
					System.out.println("trying");

					optimalSites = search.naivePrimerSearch(searchParams);

					SimpleGreedyPrimerPairPicking picker = new SimpleGreedyPrimerPairPicking();
					picker.setPrintDebugInfo(true);
					picker.setPrintTimingStatusInfo(true);
					bestPrimerPairSet = picker.pickBestPrimerSet(optimalSites, searchParams);
				}catch(EmptyResultSetException e){
					body = e.getMessage() + "\n\n" + this.searchParams.getSearchStat().printWebserviceStat() + "\n\n" + this.searchParams.getPickingStat().printWebserviceStat() + "\n\n" + paramPrinter.print3CSearchparameters(searchParams);
				}
				catch(NullPointerException e){
					if(bestPrimerPairSet == null){
						body = "-- No valid primer pair set can be picked - check your selection parameters/ change enzyme!" + "\n\n" + this.searchParams.getSearchStat().printWebserviceStat() + "\n\n" + this.searchParams.getPickingStat().printWebserviceStat() + "\n\n" + paramPrinter.print3CSearchparameters(searchParams);
					}
					else e.printStackTrace();
				}
				catch(NoRestrictionSitesFoundException e){
					body = e.getMessage() + "\n\n" + this.searchParams.getSearchStat().printWebserviceStat() + "\n\n" + this.searchParams.getPickingStat().printWebserviceStat() + "\n\n" + paramPrinter.print3CSearchparameters(searchParams);
				}
				
				if(body.equals("")){
					body = "The optimal primer pairs for your search are:\n\n" + optSetOut.printBestPrimerSet(bestPrimerPairSet) + "\n\n" + paramPrinter.print3CSearchparameters(searchParams);
					int[] features = new int[bestPrimerPairSet.size()];
					for(int  i=0; i<bestPrimerPairSet.size(); i++){
						features[i] = bestPrimerPairSet.getPrimerPair(i).getForwardPrimer().getRestrictionSite().getPosition();
					}
					if(searchParams.getNumContigs() == 1){
						// create picture showing feature positions, mail results
						FeaturePositionPicDrawer.drawImage(pictureFilename, features, searchParams.getContigs()[0].getSequenceLength());
						Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.eMail, "3CPrimerDesign - job " + this.userJobID, body, pictureFilename);
						Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.emailProperties.getAdmin(), "3CPrimerDesign - job " + this.userJobID + " for user: " + this.eMail, body, pictureFilename);
					}
					else{
						// #Contigs > 1!!!
						Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.eMail, "3CPrimerDesign - job " + this.userJobID, body);
						Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.emailProperties.getAdmin(), "3CPrimerDesign - job " + this.userJobID + " for user: " + this.eMail, body);
					}
				}
				else{
					// mail results
					Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.eMail, "3CPrimerDesign - job " + this.userJobID, body);
					Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.emailProperties.getAdmin(), "3CPrimerDesign - job " + this.userJobID + " for user: " + this.eMail, body);
				}
				
				getServletContext().log("Completed job " + this.jobID + " successfully! - " + CCCPrimerDesign.queue.size() + " elements in queue left!");
			}catch(Throwable e){
				Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.eMail, "3CPrimerDesign - job ERROR OCCURED", "An error occured during primer design - this has been reported!");
				Mail.sendMail(this.emailProperties.getSender(), this.emailProperties.getSmtp(), this.emailProperties.getAdmin(), "3CPrimerDesign - job ERROR OCCURED", "Error in primer design job: " + this.jobID + " for user: " + this.eMail + "\n\n" + returnStackTraceString(e));
				e.printStackTrace();
				getServletContext().log("Completed job " + this.jobID + " with ERRORS! - " + CCCPrimerDesign.queue.size() + " elements in queue left!");
			}
		}
	}
}
