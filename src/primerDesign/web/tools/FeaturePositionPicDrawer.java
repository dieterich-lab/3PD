/**
 * 
 */
package primerDesign.web.tools;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Random;

import javax.imageio.ImageIO;

/**
 * A class drawing a picture visualizing the placement of the optimal primers on the sequence region supplied.
 * 
 * @author Sebastian Fršhler
 *
 */
public class FeaturePositionPicDrawer {
	static Random random = new Random();
	public static void drawImage(String file, int[] values, int length) throws IOException{
		int rectangleWidth = 350;
		NumberFormat format = NumberFormat.getInstance();
		
		// init image
		BufferedImage image = new BufferedImage(450, 60, BufferedImage.TYPE_INT_RGB);
		Graphics2D graphics = image.createGraphics();
		graphics.setBackground(Color.white);
		graphics.clearRect(0, 0, 450, 60);
		graphics.setColor(Color.BLACK);
		
		// draw 'sequence'
		graphics.draw(new Rectangle(25,20, rectangleWidth, 10));
		
		// draw position foreach coordinate
		int position;
		for(int value : values){
			graphics.setColor(Color.RED);
			position =  (int)(((double)rectangleWidth) * value / length);
			graphics.draw(new Line2D.Double(25 + position, 15, 25 + position, 35));
			graphics.drawString(format.format(value), 25 + position, 55);
		}
		
		File outFile = new File(file);
		file.replace(".png", random.nextInt() + ".png");
	    ImageIO.write(image , "png", outFile );
	}
	
	public static void main(String[] args) throws IOException{
		drawImage("/Users/froehler/Desktop/PairsPisitions.png", new int[]{17064, 28721, 62513, 53399, 8397, 36853}, 68643);
	}
}
