package primerDesign.Test;

import java.util.Properties;

import javax.activation.DataHandler;
import javax.activation.DataSource;
import javax.activation.FileDataSource;
import javax.mail.BodyPart;
import javax.mail.Message;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeBodyPart;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.MimeMultipart;

public class Temp {
	public Temp() {
	  }

	  public static void main(String[] args) throws Exception {
	    Temp mailSender = new Temp();
	    mailSender.sendMail();
	  }

	  public static Properties getMailProperties() throws Exception {
	    Properties mailProperties = new Properties();
	    mailProperties.put("mail.smtp.host", "mailhost.tuebingen.mpg.de");
	    return mailProperties;
	  }

	  public void sendMail() throws Exception {
	    InternetAddress[] replyTo = {new InternetAddress("froehler@tuebingen.mpg.de")};
	    Properties mailProperties = getMailProperties();
	    Session session = Session.getDefaultInstance(mailProperties);

	    MimeMessage mimeMessage = new MimeMessage(session);
	    mimeMessage.setFrom(new InternetAddress("froehler@tuebingen.mpg.de"));
	    mimeMessage.setReplyTo(replyTo);
	    mimeMessage.addRecipient(Message.RecipientType.TO, new InternetAddress("froehler@tuebingen.mpg.de"));
	    mimeMessage.setSubject("Testing image");

	 // Create your new message part
	    BodyPart messageBodyPart = new MimeBodyPart();
	    String htmlText = "<H1>Hello</H1>" + 
	      "<img src=\"cid:memememe\">";
	    messageBodyPart.setContent(htmlText, "text/html");

	    // Create a related multi-part to combine the parts
	    MimeMultipart multipart = new MimeMultipart("related");
	    multipart.addBodyPart(messageBodyPart);

	    // Create part for the image
	    messageBodyPart = new MimeBodyPart();

	    // Fetch the image and associate to part
	    DataSource fds = new FileDataSource("/Users/froehler/Desktop/myImage.png");
	    messageBodyPart.setDataHandler(new DataHandler(fds));
	    messageBodyPart.setHeader("Content-ID","<memememe>");
	    messageBodyPart.setFileName("Dist.png");

	    // Add part to multi-part
	    multipart.addBodyPart(messageBodyPart);

	    // Associate multi-part with message
	    mimeMessage.setContent(multipart);


	    Transport.send(mimeMessage);
	  } 
	
//	/**
//	 * @param args
//	 * @throws IOException 
//	 */
//	public static void main(String[] args) throws IOException {
//		int[] array = new int[]{10,6,3,8,11};
//		for(int i=0; i<5; i++) System.out.print(array[i] + " ");
//		System.out.println();
//		
//		bucketsort(5, 5 + 1, array);
//		for(int i=0; i<5; i++) System.out.print(array[i] + " ");
//		System.out.println();
//	}
//
//	private static void bucketsort(int n, int anzahlBuckets, int z[]) {
//		   // histogramm erstellen
//		   int buckets[] = new int[anzahlBuckets];
//		   for (int i=0; i<anzahlBuckets; i++) {
//		        buckets[i] = 0;
//		   }
//		   for (int i=0; i<n; i++) {
//		        buckets[z[i]]++;
//		   }
//		   // sortieren
//		   int x=0;
//		   for (int i=0; i<anzahlBuckets; i++) {
//		        while (buckets[i] > 0) {
//		               z[x++] = i;
//		               buckets[i]--;
//		        }
//		   }
//		 }
}
