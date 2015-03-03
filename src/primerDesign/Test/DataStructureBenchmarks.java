package primerDesign.Test;

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

import primerDesign.util.SimpleTimer;
import weka.core.FastVector;
import cern.colt.list.ObjectArrayList;

public class DataStructureBenchmarks {		
	
	public static void main(String[] args) throws IllegalAlphabetException, IllegalSymbolException{
		ObjectArrayList list = new ObjectArrayList();
		FastVector vector = new FastVector();
		int iterations = 10000000;
		SimpleTimer timer = new SimpleTimer();
		
		for(int i=0; i<iterations; i++){
			vector.addElement(i);
		}
		System.out.println("Adding " + iterations + " elements to vector took " + timer.getTimeString());
		
		for(int i=0; i<iterations; i++){
			list.add(i);
		}
		System.out.println("Adding " + iterations + " elements to list took " + timer.getTimeString());
		
		
		for(int i=0; i<iterations; i++){
			vector.elementAt(i);
		}
		System.out.println("Getting " + iterations + " elements from vector took " + timer.getTimeString());
		
		for(int i=0; i<iterations; i++){
			list.get(i);
		}
		System.out.println("Getting " + iterations + " elements from list took " + timer.getTimeString());
		
		for(int i=0; i<iterations; i++){
			list.getQuick(i);
		}
		System.out.println("Getting " + iterations + " elements from quick list took " + timer.getTimeString());
	}
}
