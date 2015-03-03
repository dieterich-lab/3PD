/**
 * 
 */
package primerDesign.Test;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;

import primerDesign.util.SimpleTimer;

public class WebBlat {
	public static void main(String[] args) {

		String host = "localhost";
		String port = "80";
		String app = "cgi-bin/webBlat";

		SimpleTimer timer = new SimpleTimer();
		// for(int i=0; i<1000;i++){
		try {
			// Construct data
			// http://localhost/cgi-bin/webBlat?wb_db=Ppacificus%20%28unmasked%29&wb_qType=DNA&wb_sort=chrom,start&wb_output=psl%20no%20header&wb_seq=ATGATGATGATGATGATGATGATG
			String data = URLEncoder.encode("wb_db", "UTF-8") + "=" + URLEncoder.encode("Ppacificus (unmasked)", "UTF-8");
			data += "&" + URLEncoder.encode("wb_qType", "UTF-8") + "=" + URLEncoder.encode("DNA", "UTF-8");
			data += "&" + URLEncoder.encode("wb_sort", "UTF-8") + "=" + URLEncoder.encode("chrom,start", "UTF-8");
			data += "&" + URLEncoder.encode("wb_output", "UTF-8") + "=" + URLEncoder.encode("psl no header", "UTF-8");
//			data += "&" + URLEncoder.encode("wb_minScore", "UTF-8") + "=" + URLEncoder.encode("10", "UTF-8");
			data += "&" + URLEncoder.encode("wb_seq", "UTF-8") + "=" + URLEncoder.encode("GAGAACGCCATACAATGATAACAACAACAGAAATATTCGCTTCATTTTTAATATTCGCGCT", "UTF-8");

			// Send data
			URL url = new URL("http://" + host + ":" + port + "/" + app);
			URLConnection conn = url.openConnection();
			conn.setDoOutput(true);

			OutputStreamWriter wr = new OutputStreamWriter(conn
					.getOutputStream());
			wr.write(data);
			wr.flush();

			// Get the response
			BufferedReader rd = new BufferedReader(new InputStreamReader(conn.getInputStream()));
			String line;
			while ((line = rd.readLine()) != null) {
				// Process line...
				System.out.println(line); 
			}
			wr.close();
			rd.close();
		} catch (Exception e) {
			
		}
		// }
		System.out.println(timer.getTimeString());
	}
}