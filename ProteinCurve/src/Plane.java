import javax.vecmath.Vector3d;

/**
 * A class to store the planes for the convex hull. Used to help determine if a
 * point is inside or outside the convex hull.
 * 
 * @author Kyle Diller
 *
 */
public class Plane {
	/**
	 * The point on the plane that is used as a reference for computing distance
	 * to the plane.
	 */
	private Vector3d p0;

	/**
	 * The normal vector to the plane that is used to help compute the distance
	 * to the plane.
	 */
	private Vector3d norm;

	/**
	 * One of the three atoms that make up the plane.
	 */
	public Atom3D a1, a2, a3;

	/**
	 * Creates a plane with a starting point, normal vector and 3 atoms that
	 * define it.
	 * 
	 * @param k
	 *            the starting point in the plane.
	 * @param a
	 *            the normal vector to the plane.
	 * @param b1
	 *            the first atom.
	 * @param b2
	 *            the second atom.
	 * @param b3
	 *            the third atom.
	 */
	public Plane(Vector3d k, Vector3d a, Atom3D b1, Atom3D b2, Atom3D b3) {
		norm = new Vector3d(a);
		p0 = new Vector3d(k);
		a1 = b1;
		a2 = b2;
		a3 = b3;
	}

	/**
	 * Calculates the closest distance from the plane to the point. <br>
	 * Follows: <br>
	 * distance = ((x, y, z) - p0) * norm
	 * 
	 * @param x
	 *            the x component of the point.
	 * @param y
	 *            the y component of the point.
	 * @param z
	 *            the z component of the point.
	 * @return the closest distance from the plane to the point.
	 */
	public double dis(double x, double y, double z) {
		Vector3d temp = new Vector3d(x, y, z);
		temp.sub(p0);
		double side = temp.dot(norm);
		return side;
	}

	/**
	 * Returns the closest point on the plane to an atom. <br>
	 * Follows: <br>
	 * p = atomLocation - atomRadius * norm
	 * 
	 * @param a
	 *            the atom to find the closest point to.
	 * @return the closest point of the plane to the atom.
	 */
	public Vector3d pointOnPlane(Atom3D a) {
		Vector3d adda = new Vector3d(norm);
		adda.x *= a.r;
		adda.y *= a.r;
		adda.z *= a.r;

		Vector3d temp = new Vector3d();
		temp.sub(a.getVector(), adda);
		return temp;
	}

	public String toString() {
		return p0 + " <> " + norm;
	}
}
