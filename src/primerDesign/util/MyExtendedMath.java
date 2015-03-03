package primerDesign.util;

/**
 * Implements some math functions additionally to the Math class.
 * 
 * @author Sebastian Fršhler
 *
 */
public class MyExtendedMath {
	/**
	 * Returns the closest integer to value 'value' that is less than or equal than value.
	 * 
	 * Convenience method implementing Math.floor for integers
	 * 
	 * @param value the value to round down to an integer
	 * 
	 * @return the integer value that is less than or equal to value
	 */
	public static int floor(double value){
		return (int) Math.floor(value);
	}
	
	/**
	 * Returns the closest integer to value 'value' that is larger than value.
	 * 
	 * Convenience method implementing Math.ceil for integers
	 * 
	 * @param value the value to round up to an integer
	 * 
	 * @return the integer value that is larger than value
	 */
	public static int ceil(double value){
		return (int) Math.ceil(value);
	}
	
	/**
	 * Implements the Math.round method AND returning integer.
	 * 
	 * @param x the number to round
	 * 
	 * @return the 'closest' integer to x
	 */
	public static int round(double x){
		double round = x - MyExtendedMath.floor(x);
		if(round < 0.5) return MyExtendedMath.floor(x);
		else return MyExtendedMath.ceil(x);
	}
}
