/**
 * 
 */
package primerDesign.dsc.indexStructures.rmi;

import java.io.File;
import java.rmi.Naming;

/**
 * @author Sebastian Fršhler
 *
 */
public class IndexRegistration {
	public static void main(String[] args){
	    //System.setSecurityManager(new RMISecurityManager()); 
	    try {
	      System.out.println("Reading in index");
	      RemoteIndexStructureSearchESAImpl index = new RemoteIndexStructureSearchESAImpl(new File(args[0]));
	      System.out.println("Performing sample lookup of 'toy' string ATG, found " + index.findHitPositions("ATG").size() + " matches");
	      System.out.println("Registering RemoteIndex in rmiregistry");
	      Naming.rebind("RemoteIndex", index);
	      System.out.println("  Done.");
	    } catch (Exception e) {
	      //System.err.println(e.toString());
	      e.printStackTrace();
	      System.exit(1);
	    }
	}
}
