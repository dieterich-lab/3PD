package primerDesign.util;

/**
 * Thrown if no hits can be found in an index but hits were expected.
 * +
 * @author Sebastian Fršhler
 *
 */
public class NoHitsInIndexException extends RuntimeException {

	private static final long serialVersionUID = -6746207910259496871L;

	public NoHitsInIndexException() {
		super();
	}

	public NoHitsInIndexException(String message) {
		super(message);
	}
}