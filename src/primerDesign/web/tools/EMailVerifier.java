package primerDesign.web.tools;

import java.util.regex.Pattern;

import javax.mail.internet.AddressException;
import javax.mail.internet.InternetAddress;

/**
 * Simple class verifying the correctness of a mail address.
 * 
 * @author Sebastian Fršhler
 *
 */
public class EMailVerifier {
	
	static final Pattern topLevelDomains = Pattern.compile("[A-Za-z]{2}|com|org|net|gov|mil|biz|info|mobi|name|aero|jobs|museum|edu");
	
	/**
	 * Check whether an email address is valid according to RCF822.
	 * 
	 * Only mail address of the form <address>@<host>.<tld> are considered being valid global mail addresses!
	 * 
	 * @param aEmailAddress the address to verify
	 * 
	 * @return true iff an email address is valid according to RCF822
	 */
	public static boolean isValidGlobalEmailAddress(String aEmailAddress){
	    return isValidLocalEmailAdress(aEmailAddress) && hasNameAndDomain(aEmailAddress);
	}
	
	/**
	 * Checks whether an email address is valid according to RCF822 allowing local adresses.
	 * 
	 * @param aEmailAddress the address to verify
	 * 
	 * @return true iff an email address is valid according to RCF822
	 */
	public static boolean isValidLocalEmailAdress(String aEmailAddress){
		if (aEmailAddress == null) return false;
	    boolean result = true;
	    try {
	      InternetAddress emailAddr = new InternetAddress(aEmailAddress, true);
	      emailAddr.validate();
	    }
	    catch (AddressException ex){
	      result = false;
	    }
	    catch(Exception e){
	    	return false;
	    }
	    return result;
	}

	  private static boolean hasNameAndDomain(String aEmailAddress){
	    String[] tokens = aEmailAddress.split("@");
	    String[] tld = tokens[1].split(Pattern.quote("."));
	    return tokens.length == 2 && tokens[0].length() > 0 && tokens[1].length() >= 2 && topLevelDomains.matcher(tld[tld.length - 1]).matches();
	  }
	
	public static void main(String[] args){
		System.out.println(EMailVerifier.isValidGlobalEmailAddress("sebastian.froehler@tuebingen.mpg.de"));
		System.out.println(EMailVerifier.isValidGlobalEmailAddress("Adhideb.Ghosh@age.mpg.de"));
		System.out.println(EMailVerifier.isValidLocalEmailAdress("sebi"));
		System.out.println(EMailVerifier.isValidLocalEmailAdress("sebi@localhost"));		
		System.out.println(EMailVerifier.isValidLocalEmailAdress("MUMMIDI@uthscsa.edu"));	
		System.out.println(EMailVerifier.isValidGlobalEmailAddress("MUMMIDI@uthscsa.edu"));
	}
}
