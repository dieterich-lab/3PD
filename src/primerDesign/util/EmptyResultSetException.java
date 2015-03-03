package primerDesign.util;

/**
 * Thrown if the result set of a computation is empty.
 * +
 * @author Sebastian Fršhler
 *
 */
public class EmptyResultSetException extends RuntimeException {

	private static final long serialVersionUID = 1579907120361003524L;

	public EmptyResultSetException() {
		super();
	}

	public EmptyResultSetException(String message) {
		super(message);
	}
}
