/**
 * 
 */
package primerDesign.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Class for transforming files in the LightCycler 480 csv format.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SimpleLC480csvTransformer {
	private static final int wellIndex = 0;
	private static final int measurementIndex = 2;
	private static final int cycleIndex = 4;
	private static final int intensityIndex = 7;
	
	public static void main(String[] args) throws IOException{
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		
		int cycles = Integer.parseInt(args[1]);
		int plateSize = Integer.parseInt(args[2]);
		
		String[][] matrix = new String[cycles + 1][plateSize + 1];
		String line;
		String[] values;
		int columns = 1;
		String lastWell = "";
		
		matrix[0][0] = "Cycles";
		
		for(int i=1; i<= cycles; i++){
			matrix[i][0] = Integer.toString(i);
		}
		
		while((line = reader.readLine()) != null){
			// if line is a measurement and NOT header
			if(line.matches("[A-H][0-9]+.*")){
				values = line.split("\t");
				if(values.length != 8) continue;
				if(lastWell == "") lastWell = values[wellIndex];
				else if(!lastWell.equals(values[wellIndex])){
					lastWell = values[wellIndex];
					columns++;
				}
				// if this line is a quantification measurement
				if(values[measurementIndex].equals("2")){
					// set value
					matrix[Integer.parseInt(values[cycleIndex])][columns] = values[intensityIndex];
					matrix[0][columns] = values[wellIndex];

				}
			}
		}
		
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[0].length; j++){
				System.out.print(matrix[i][j]);
				if(j < matrix[0].length - 1) System.out.print("\t");
			}
			System.out.println();
		}
	}
}
