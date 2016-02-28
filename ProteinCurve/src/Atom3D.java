import java.util.ArrayList;

import javax.vecmath.Vector3d;

public class Atom3D implements Comparable<Atom3D> {
	public static Atom3D readAtom(String str) {
		String[] atm = str.split(",");

		return new Atom3D(Double.parseDouble(atm[1]),
				Double.parseDouble(atm[2]), Double.parseDouble(atm[3]),
				Double.parseDouble(atm[4]));
	}

	public double x, y, z, r;
	public boolean buried;
	public ArrayList<Atom3D> neigh;

	public Atom3D(double a, double b, double c, double d) {
		x = a;
		y = b;
		z = c;
		r = d;

		neigh = new ArrayList<Atom3D>();
	}

	public double getDis2D(double x, double y) {
		return Math.sqrt((this.x - x) * (this.x - x) + (this.y - y)
				* (this.y - y));
	}

	public String toString() {
		return x + "," + y + "," + z + "," + r;
	}

	public double dis2D(Atom3D o) {
		return getDis2D(o.x, o.y);
	}

	public double getDis3D(double x, double y, double z) {
		return Math.sqrt((this.x - x) * (this.x - x) + (this.y - y)
				* (this.y - y) + (this.z - z) * (this.z - z));
	}

	public double dis3D(Atom3D o) {
		double dis = getDis3D(o.x, o.y, o.z);
		return dis;
	}

	public double disSQ(double x, double y, double z) {
		return (this.x - x) * (this.x - x) + (this.y - y) * (this.y - y)
				+ (this.z - z) * (this.z - z);
	}

	public Vector3d getVector() {
		return new Vector3d(x, y, z);
	}

	@Override
	public int compareTo(Atom3D arg0) {
		return Double.compare(this.r, arg0.r);
	}
}
