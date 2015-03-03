package primerDesign.util;

/**
 * Thrown if the result set of a computation is empty.
 * +
 * @author Sebastian Fršhler
 *
 */
public class NoRestrictionSitesFoundException extends RuntimeException {

	private static final long serialVersionUID = 1579907120361003524L;

	public NoRestrictionSitesFoundException() {
		super();
	}

	public NoRestrictionSitesFoundException(String message) {
		super(message);
	}
}
