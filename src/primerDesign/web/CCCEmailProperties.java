/**
 * 
 */
package primerDesign.web;

/**
 * Class storing mail preferences.
 * 
 * @author Sebastian Fršhler
 *
 */
public class CCCEmailProperties {
	private final String sender;
	private final String admin;
	private final String smtp;
	
	public CCCEmailProperties(String sender, String admin, String smtp){
		this.sender = sender;
		this.admin = admin;
		this.smtp = smtp;
	}

	/**
	 * @return the sender
	 */
	public String getSender() {
		return sender;
	}

	/**
	 * @return the admin
	 */
	public String getAdmin() {
		return admin;
	}

	/**
	 * @return the smtp
	 */
	public String getSmtp() {
		return smtp;
	}
}
