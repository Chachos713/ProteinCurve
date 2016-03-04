import java.util.ArrayList;

import javax.vecmath.Vector3d;

/**
 * A data structure used to store the location radius, whether an atom is buried
 * in the protein, and it's neighbors.
 * 
 * @author Kyle Diller
 *
 */
public class Atom3D {
	/**
	 * Reads an atom from a csv file.
	 * 
	 * @param str
	 *            a line from the csv file.
	 * @return the atom that is read from the line.
	 */
	public static Atom3D readAtom(String str) {
		String[] atm = str.split(",");

		return new Atom3D(Double.parseDouble(atm[1]),
				Double.parseDouble(atm[2]), Double.parseDouble(atm[3]),
				Double.parseDouble(atm[4]));
	}

	/**
	 * One of the values to describe an atom's location or the sace it take up.
	 */
	public double x, y, z, r;

	/**
	 * Used to help know if the atom is burried within the protein.
	 */
	public boolean buried;

	/**
	 * List of neighbors that could cause the atom to be buried or not.
	 */
	public ArrayList<Atom3D> neigh;

	/**
	 * Creates a new atom with x, y, z and radius component.
	 * 
	 * @param a
	 *            the x component.
	 * @param b
	 *            the y component.
	 * @param c
	 *            the z component.
	 * @param d
	 *            the radius of the atom.
	 */
	public Atom3D(double a, double b, double c, double d) {
		x = a;
		y = b;
		z = c;
		r = d;

		neigh = new ArrayList<Atom3D>();
	}

	public String toString() {
		return x + "," + y + "," + z + "," + r;
	}

	/**
	 * Calculates the distance from an atom to a 3 dimensional point in space.
	 * 
	 * @param x
	 *            the x component of the point.
	 * @param y
	 *            the y compoennt of the point.
	 * @param z
	 *            the z component of the point.
	 * @return the distance between the atom and the point.
	 */
	public double getDis3D(double x, double y, double z) {
		return Math.sqrt(disSQ(x, y, z));
	}

	/**
	 * Computes the distance between two atoms.
	 * 
	 * @param o
	 *            the other atom.
	 * @return the distance between this atom and another atom.
	 */
	public double dis3D(Atom3D o) {
		double dis = getDis3D(o.x, o.y, o.z);
		return dis;
	}

	/**
	 * Computes the distance squared between an atom and somep oint in space.
	 * 
	 * @param x
	 *            the x component of the point.
	 * @param y
	 *            the y component of the point.
	 * @param z
	 *            the z component of the point.
	 * @return the distance squared between an atom and a point.
	 */
	public double disSQ(double x, double y, double z) {
		return (this.x - x) * (this.x - x) + (this.y - y) * (this.y - y)
				+ (this.z - z) * (this.z - z);
	}

	/**
	 * Turns the location of the atom into a vector.
	 * 
	 * @return the vector representation of the atom's location.
	 */
	public Vector3d getVector() {
		return new Vector3d(x, y, z);
	}
}
