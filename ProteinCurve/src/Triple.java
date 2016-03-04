/**
 * A storage class to maintain 3 atoms that for a sphere with the query point.
 * This is more or less a place holder for a more complex data structure to
 * maintain which atoms are grouped together.
 * 
 * @author Kyle Diller
 *
 */
public class Triple {
	/**
	 * The three atoms that are grouped together.
	 */
	public Atom3D a1, a2, a3;

	/**
	 * Creates a grouping of three atoms with some initial references.
	 * 
	 * @param a
	 *            the first atom in the grouping.
	 * @param b
	 *            the second atom in the grouping.
	 * @param c
	 *            the third atom in ther grouping.
	 */
	public Triple(Atom3D a, Atom3D b, Atom3D c) {
		a1 = a;
		a2 = b;
		a3 = c;
	}
}
