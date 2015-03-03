package primerDesign.util;

/**
 * Thrown if a restriction site is used multiple times in different 'SequenceRegion's.
 * 
 * @author Sebastian Fršhler
 *
 */
public class DuplicateUseOfRestrictionSiteException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 5973755450725425002L;

	public DuplicateUseOfRestrictionSiteException() {
		super();
	}

	public DuplicateUseOfRestrictionSiteException(String arg0) {
		super(arg0);
	}

}
