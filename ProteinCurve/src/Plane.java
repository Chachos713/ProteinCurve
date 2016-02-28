import javax.vecmath.Vector3d;

public class Plane {
	private Vector3d p0, norm;
	public Atom3D a1, a2, a3;

	public Plane(Vector3d k, Vector3d a, Atom3D b1, Atom3D b2, Atom3D b3) {
		norm = new Vector3d(a);
		p0 = new Vector3d(k);
		a1 = b1;
		a2 = b2;
		a3 = b3;
	}

	public double dis(double x, double y, double z) {
		Vector3d temp = new Vector3d(x, y, z);
		temp.sub(p0);
		double side = temp.dot(norm);
		return side;
	}

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
		// return p0 + " <> " + norm;
		return p0.x + "," + p0.y + "," + p0.z + "," + norm.x + "," + norm.y
				+ "," + norm.z;
	}
}
