package primerDesign.web.tools;

import java.util.Date;
import java.util.Properties;

import javax.activation.DataHandler;
import javax.activation.DataSource;
import javax.activation.FileDataSource;
import javax.mail.BodyPart;
import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.NoSuchProviderException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeBodyPart;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.MimeMultipart;

/**
 * A class wrapping mail functionality for 3PD result mails.
 * 
 * This class tests the JavaMail API - which depends on the JavaBeans Activation Framework!
 * 
 * @author Sebastian Fr�hler
 */
public class Mail {
	public static void main(String[] args){
		//Mail.sendMail("froehler@tuebingen.mpg.de", "mailhost.tuebingen.mpg.de", "Sebastian.Froehler@tuebingen.mpg.de", "JavaMail test", "This is just a test of the JavaMail API!");
		Mail.sendMail("threepd@age.mpg.de", "mail.age.mpg.de", "threepd@age.mpg.de", "JavaMail test", "This is just a test of the JavaMail API!");
	}
	
	/**
	 * Convenience method to send a simple 'plain' text message encoded as UTF-8.
	 * 
	 * @param from the sender's adress
	 * @param smtp_host the smtp host to send messages with
	 * @param to the receiver's adress
	 * @param subject the subject of the message
	 * @param body the body of the message
	 */
	public static void sendMail(String from, String smtp_host, String to, String subject, String body){
		sendMail(from, smtp_host, to, subject, body, "UTF-8", "plain");
	}
	
	/**
	 * Convenience method to send a simple 'plain' text message encoded as UTF-8.
	 * 
	 * @param from the sender's adress
	 * @param smtp_host the smtp host to send messages with
	 * @param to the receiver's adress
	 * @param subject the subject of the message
	 * @param body the body of the message
	 * @param the filename to attach to the message
	 */
	public static void sendMail(String from, String smtp_host, String to, String subject, String body, String filename) throws Exception{
		sendMail(from, smtp_host, to, subject, body, filename, "UTF-8", "plain");
	}
	
	/**
	 * Detailed method for sending  a message.
	 * 
	 * @param from the sender's adress
	 * @param smtp_host the smtp host to send messages with
	 * @param to the receiver's adress
	 * @param subject the subject of the message
	 * @param body the body of the message
	 * @param encoding the character encoding to be used
	 * @param contentType the content type to be used
	 */
	public static void sendMail(String from, String smtp_host, String to, String subject, String body, String encoding, String contentType){
		
		
		
	    String host = "mail.age.mpg.de";
	    String username = "USER";
	    String password = "PW";
	    Properties props = new Properties();
	    props.put("mail.smtp.ssl.enable","TRUE");
	    props.put("mail.smpt.starttls.enable", "TRUE");
	    props.put("mail.smtp.port", "587");
	    props.put("mail.from", from);
	    props.put("mail.smtp.auth", "true");
	    props.put("mail.debug","false");   	    

	    // set any needed mail.smtps.* properties here
	    Session session = Session.getInstance(props);
	    MimeMessage msg = new MimeMessage(session);
        try {
			msg.setFrom();
	        msg.setRecipients(Message.RecipientType.TO, to);
	        msg.setSubject(subject);
	        msg.setSentDate(new Date());
	        msg.setText(body, encoding, contentType);

		} catch (MessagingException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

	    // set the message content here
	    Transport t = null;
		try {
			t = session.getTransport("smtps");
		} catch (NoSuchProviderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    try {
		try {
			t.connect(host, username, password);

		} catch (MessagingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			t.sendMessage(msg, msg.getAllRecipients());

		} catch (MessagingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    } finally {
		try {
			t.close();
		} catch (MessagingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    }
		
	}


	/**
	 * Detailed method for sending  a message.
	 * 
	 * @param from the sender's adress
	 * @param smtp_host the smtp host to send messages with
	 * @param to the receiver's adress
	 * @param subject the subject of the message
	 * @param body the body of the message
	 * @param encoding the character encoding to be used
	 * @param contentType the content type to be used
	 * 
	 * @throws Exception
	 */
	  public static void sendMail(String from, String smtp_host, String to, String subject, String body, String filename, String encoding, String contentType) throws Exception {
	    InternetAddress[] replyTo = {new InternetAddress(from)};
	    Properties mailProperties = new Properties();
	    mailProperties.put("mail.smtp.host", smtp_host);
	    Session session = Session.getDefaultInstance(mailProperties);

	    MimeMessage mimeMessage = new MimeMessage(session);
	    mimeMessage.setFrom(new InternetAddress(from));
	    mimeMessage.setReplyTo(replyTo);
	    mimeMessage.addRecipient(Message.RecipientType.TO, new InternetAddress(to));
	    mimeMessage.setSubject(subject);

	 // Create your new message part
	    BodyPart messageBodyPart = new MimeBodyPart();
	    String htmlText = "<img src=\"cid:memememe\">\n" + "<pre>" + body + "</pre>";
	    messageBodyPart.setContent(htmlText, "text/html");

	    // Create a related multi-part to combine the parts
	    MimeMultipart multipart = new MimeMultipart("related");
	    multipart.addBodyPart(messageBodyPart);

	    // Create part for the image
	    messageBodyPart = new MimeBodyPart();

	    // Fetch the image and associate to part
	    DataSource fds = new FileDataSource(filename);
	    messageBodyPart.setDataHandler(new DataHandler(fds));
	    messageBodyPart.setHeader("Content-ID","<memememe>");
	    messageBodyPart.setFileName("picture.png");

	    // Add part to multi-part
	    multipart.addBodyPart(messageBodyPart);
	    
//	    messageBodyPart = new MimeBodyPart();
//	    htmlText = body;
//	    messageBodyPart.setContent(htmlText, "text/html");
//	    multipart.addBodyPart(messageBodyPart);

	    // Associate multi-part with message
	    mimeMessage.setContent(multipart);

	    Transport.send(mimeMessage);
	  } 
}
