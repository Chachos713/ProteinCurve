import Jama.Matrix;

public class Curvature3D {
	public static Atom3D[] intersectingSphere(Atom3D a1, Atom3D a2, Atom3D a3,
			Atom3D a4) {
		try {
			Matrix r = new Matrix(3, 3);
			Matrix x = new Matrix(3, 3);
			Matrix y = new Matrix(3, 3);
			Matrix z = new Matrix(3, 3);

			r.set(0, 0, -2 * (a1.x - a2.x));
			r.set(1, 0, -2 * (a1.x - a3.x));
			r.set(2, 0, -2 * (a1.x - a4.x));
			r.set(0, 1, -2 * (a1.y - a2.y));
			r.set(1, 1, -2 * (a1.y - a3.y));
			r.set(2, 1, -2 * (a1.y - a4.y));
			r.set(0, 2, -2 * (a1.z - a2.z));
			r.set(1, 2, -2 * (a1.z - a3.z));
			r.set(2, 2, -2 * (a1.z - a4.z));

			x.set(0, 0, -2 * (a1.r - a2.r));
			x.set(1, 0, -2 * (a1.r - a3.r));
			x.set(2, 0, -2 * (a1.r - a4.r));
			x.set(0, 1, -2 * (a1.y - a2.y));
			x.set(1, 1, -2 * (a1.y - a3.y));
			x.set(2, 1, -2 * (a1.y - a4.y));
			x.set(0, 2, -2 * (a1.z - a2.z));
			x.set(1, 2, -2 * (a1.z - a3.z));
			x.set(2, 2, -2 * (a1.z - a4.z));

			y.set(0, 0, -2 * (a1.x - a2.x));
			y.set(1, 0, -2 * (a1.x - a3.x));
			y.set(2, 0, -2 * (a1.x - a4.x));
			y.set(0, 1, -2 * (a1.r - a2.r));
			y.set(1, 1, -2 * (a1.r - a3.r));
			y.set(2, 1, -2 * (a1.r - a4.r));
			y.set(0, 2, -2 * (a1.z - a2.z));
			y.set(1, 2, -2 * (a1.z - a3.z));
			y.set(2, 2, -2 * (a1.z - a4.z));

			z.set(0, 0, -2 * (a1.x - a2.x));
			z.set(1, 0, -2 * (a1.x - a3.x));
			z.set(2, 0, -2 * (a1.x - a4.x));
			z.set(0, 1, -2 * (a1.y - a2.y));
			z.set(1, 1, -2 * (a1.y - a3.y));
			z.set(2, 1, -2 * (a1.y - a4.y));
			z.set(0, 2, -2 * (a1.r - a2.r));
			z.set(1, 2, -2 * (a1.r - a3.r));
			z.set(2, 2, -2 * (a1.r - a4.r));

			Matrix temp = r;
			double det = Math.abs(r.det());

			double detTemp = Math.abs(x.det());
			if (det < detTemp) {
				det = detTemp;
				temp = x;
			}

			detTemp = Math.abs(y.det());
			if (det < detTemp) {
				det = detTemp;
				temp = y;
			}

			detTemp = Math.abs(z.det());
			if (det < detTemp) {
				det = detTemp;
				temp = z;
			}

			if (temp == z) {
				return solve4Z(z, a1, a2, a3, a4);
			} else if (temp == y) {
				return solve4Y(y, a1, a2, a3, a4);
			} else if (temp == x) {
				return solve4X(x, a1, a2, a3, a4);
			} else {
				return solve4R(r, a1, a2, a3, a4);
			}
		} catch (Exception e) {
		}

		return null;
	}

	private static Atom3D[] solve4Z(Matrix A, Atom3D a1, Atom3D a2, Atom3D a3,
			Atom3D a4) {
		Matrix B = A.transpose().times(A);
		Matrix C = new Matrix(3, 1);
		Matrix D = new Matrix(3, 1);

		D.set(0, 0, 2 * (a1.z - a2.z));
		D.set(1, 0, 2 * (a1.z - a3.z));
		D.set(2, 0, 2 * (a1.z - a4.z));

		C.set(0, 0, a2.x * a2.x - a1.x * a1.x + a2.y * a2.y - a1.y * a1.y
				- a2.r * a2.r + a1.r * a1.r);
		C.set(1, 0, a3.x * a3.x - a1.x * a1.x + a3.y * a3.y - a1.y * a1.y
				- a3.r * a3.r + a1.r * a1.r);
		C.set(2, 0, a4.x * a4.x - a1.x * a1.x + a4.y * a4.y - a1.y * a1.y
				- a4.r * a4.r + a1.r * a1.r);

		Matrix cR = B.inverse().times(A.transpose()).times(C);
		Matrix dR = B.inverse().times(A.transpose()).times(D);

		double Ax = cR.get(0, 0);
		double Ay = cR.get(1, 0);
		double Ar = cR.get(2, 0);

		double Bx = dR.get(0, 0);
		double By = dR.get(1, 0);
		double Br = dR.get(2, 0);

		double a = 1 + By * By + Bx * Bx - Br * Br;

		double b = -2 * a1.z + 2 * Ay * By - 2 * a1.y * By + 2 * Ax * Bx - 2
				* a1.x * Bx - 2 * Ar * Br + 2 * a1.r * Br;

		double c = a1.z * a1.z + Ay * Ay + Ax * Ax - Ar * Ar - 2 * a1.y * Ay
				- 2 * a1.x * Ax + 2 * a1.r * Ar;

		double[] radii = Util.quadratic(a, b, c);
		Matrix vari;

		if (radii == null) {
			return null;
		}

		Atom3D[] a3D = new Atom3D[radii.length];

		for (int i = 0; i < a3D.length; i++) {
			vari = dR.times(radii[i]).minus(cR);
			a3D[i] = new Atom3D(vari.get(0, 0), vari.get(1, 0), radii[i],
					vari.get(2, 0));
		}

		check(a1, a2, a3, a4, a3D);

		return a3D;
	}

	private static Atom3D[] solve4Y(Matrix A, Atom3D a1, Atom3D a2, Atom3D a3,
			Atom3D a4) {

		Matrix B = A.transpose().times(A);
		Matrix C = new Matrix(3, 1);
		Matrix D = new Matrix(3, 1);

		D.set(0, 0, 2 * (a1.y - a2.y));
		D.set(1, 0, 2 * (a1.y - a3.y));
		D.set(2, 0, 2 * (a1.y - a4.y));

		C.set(0, 0, a2.x * a2.x - a1.x * a1.x + a2.y * a2.y - a1.y * a1.y
				- a2.r * a2.r + a1.r * a1.r + a2.z * a2.z - a1.z * a1.z);
		C.set(1, 0, a3.x * a3.x - a1.x * a1.x + a3.y * a3.y - a1.y * a1.y
				- a3.r * a3.r + a1.r * a1.r + a3.z * a3.z - a1.z * a1.z);
		C.set(2, 0, a4.x * a4.x - a1.x * a1.x + a4.y * a4.y - a1.y * a1.y
				- a4.r * a4.r + a1.r * a1.r + a4.z * a4.z - a1.z * a1.z);

		Matrix cR = B.inverse().times(A.transpose()).times(C);
		Matrix dR = B.inverse().times(A.transpose()).times(D);

		double Ax = cR.get(0, 0);
		double Ar = cR.get(1, 0);
		double Az = cR.get(2, 0);

		double Bx = dR.get(0, 0);
		double Br = dR.get(1, 0);
		double Bz = dR.get(2, 0);

		double a = 1 + Bx * Bx + Bz * Bz - Br * Br;

		double b = -2 * a1.y + 2 * Ax * Bx - 2 * a1.x * Bx + 2 * Az * Bz - 2
				* a1.z * Bz - 2 * Ar * Br + 2 * a1.r * Br;

		double c = a1.y * a1.y + Ax * Ax + Az * Az - Ar * Ar - 2 * a1.x * Ax
				- 2 * a1.z * Az + 2 * a1.r * Ar;

		double[] radii = Util.quadratic(a, b, c);
		Matrix vari;

		if (radii == null) {
			return null;
		}

		Atom3D[] a3D = new Atom3D[radii.length];

		for (int i = 0; i < a3D.length; i++) {
			vari = dR.times(radii[i]).minus(cR);
			a3D[i] = new Atom3D(vari.get(0, 0), radii[i], vari.get(2, 0),
					vari.get(1, 0));
		}

		check(a1, a2, a3, a4, a3D);

		return a3D;
	}

	private static Atom3D[] solve4X(Matrix A, Atom3D a1, Atom3D a2, Atom3D a3,
			Atom3D a4) {

		Matrix B = A.transpose().times(A);
		Matrix C = new Matrix(3, 1);
		Matrix D = new Matrix(3, 1);

		D.set(0, 0, 2 * (a1.x - a2.x));
		D.set(1, 0, 2 * (a1.x - a3.x));
		D.set(2, 0, 2 * (a1.x - a4.x));

		C.set(0, 0, a2.x * a2.x - a1.x * a1.x + a2.y * a2.y - a1.y * a1.y
				- a2.r * a2.r + a1.r * a1.r + a2.z * a2.z - a1.z * a1.z);
		C.set(1, 0, a3.x * a3.x - a1.x * a1.x + a3.y * a3.y - a1.y * a1.y
				- a3.r * a3.r + a1.r * a1.r + a3.z * a3.z - a1.z * a1.z);
		C.set(2, 0, a4.x * a4.x - a1.x * a1.x + a4.y * a4.y - a1.y * a1.y
				- a4.r * a4.r + a1.r * a1.r + a4.z * a4.z - a1.z * a1.z);

		Matrix cR = B.inverse().times(A.transpose()).times(C);
		Matrix dR = B.inverse().times(A.transpose()).times(D);

		double Ar = cR.get(0, 0);
		double Ay = cR.get(1, 0);
		double Az = cR.get(2, 0);

		double Br = dR.get(0, 0);
		double By = dR.get(1, 0);
		double Bz = dR.get(2, 0);

		double a = 1 + By * By + Bz * Bz - Br * Br;

		double b = -2 * a1.x + 2 * Ay * By - 2 * a1.y * By + 2 * Az * Bz - 2
				* a1.z * Bz - 2 * Ar * Br + 2 * a1.r * Br;

		double c = a1.x * a1.x + Ay * Ay + Az * Az - Ar * Ar - 2 * a1.y * Ay
				- 2 * a1.z * Az + 2 * a1.r * Ar;

		double[] radii = Util.quadratic(a, b, c);
		Matrix vari;

		if (radii == null) {
			return null;
		}

		Atom3D[] a3D = new Atom3D[radii.length];

		for (int i = 0; i < a3D.length; i++) {
			vari = dR.times(radii[i]).minus(cR);
			a3D[i] = new Atom3D(radii[i], vari.get(1, 0), vari.get(2, 0),
					vari.get(0, 0));
		}

		check(a1, a2, a3, a4, a3D);

		return a3D;
	}

	private static Atom3D[] solve4R(Matrix A, Atom3D a1, Atom3D a2, Atom3D a3,
			Atom3D a4) {

		Matrix B = A.transpose().times(A);
		Matrix C = new Matrix(3, 1); // Constants
		Matrix D = new Matrix(3, 1); // Slope

		D.set(0, 0, 2 * (a1.r - a2.r));
		D.set(1, 0, 2 * (a1.r - a3.r));
		D.set(2, 0, 2 * (a1.r - a4.r));

		C.set(0, 0, a2.x * a2.x - a1.x * a1.x + a2.y * a2.y - a1.y * a1.y
				- a2.r * a2.r + a1.r * a1.r + a2.z * a2.z - a1.z * a1.z);
		C.set(1, 0, a3.x * a3.x - a1.x * a1.x + a3.y * a3.y - a1.y * a1.y
				- a3.r * a3.r + a1.r * a1.r + a3.z * a3.z - a1.z * a1.z);
		C.set(2, 0, a4.x * a4.x - a1.x * a1.x + a4.y * a4.y - a1.y * a1.y
				- a4.r * a4.r + a1.r * a1.r + a4.z * a4.z - a1.z * a1.z);

		Matrix cR = B.inverse().times(A.transpose()).times(C);
		Matrix dR = B.inverse().times(A.transpose()).times(D);

		double Ax = cR.get(0, 0);
		double Ay = cR.get(1, 0);
		double Az = cR.get(2, 0);

		double Bx = dR.get(0, 0);
		double By = dR.get(1, 0);
		double Bz = dR.get(2, 0);

		double a = Bx * Bx + By * By + Bz * Bz - 1;
		double b = 2 * Ax * Bx - 2 * a1.x * Bx + 2 * Ay * By - 2 * a1.y * By
				+ 2 * Az * Bz - 2 * a1.z * Bz - 2 * a1.r;
		double c = Ax * Ax + a1.x * a1.x - 2 * a1.x * Ax + Ay * Ay + a1.y
				* a1.y - 2 * a1.y * Ay + Az * Az + a1.z * a1.z - 2 * a1.z * Az
				- a1.r * a1.r;

		double[] radii = Util.quadratic(a, b, c);
		Matrix vari;

		if (radii == null) {
			return null;
		}

		Atom3D[] a3D = new Atom3D[radii.length];

		for (int i = 0; i < a3D.length; i++) {
			vari = dR.times(radii[i]).plus(cR);
			a3D[i] = new Atom3D(vari.get(0, 0), vari.get(1, 0), vari.get(2, 0),
					radii[i]);
		}

		check(a1, a2, a3, a4, a3D);

		return a3D;
	}

	private static boolean check(Atom3D a1, Atom3D a2, Atom3D a3, Atom3D a4,
			Atom3D[] a3d) {

		boolean good = false;

		for (int l = 0; l < a3d.length; l++) {
			good = true;
			if (Math.abs(a1.dis3D(a3d[l]) - (a1.r + a3d[l].r)) > Util.ERROR)
				good = false;
			if (Math.abs(a2.dis3D(a3d[l]) - (a2.r + a3d[l].r)) > Util.ERROR)
				good = false;
			if (Math.abs(a3.dis3D(a3d[l]) - (a3.r + a3d[l].r)) > Util.ERROR)
				good = false;
			if (Math.abs(a4.dis3D(a3d[l]) - (a4.r + a3d[l].r)) > Util.ERROR)
				good = false;

			if (!good) {
				a3d[l] = null;
			}
		}

		return true;
	}
}
