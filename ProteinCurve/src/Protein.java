import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

import javax.vecmath.Vector3d;

public class Protein {
	private ArrayList<Atom3D> atoms;
	private ArrayList<Plane> hull;

	private Protein() {
		atoms = new ArrayList<Atom3D>();
		hull = new ArrayList<Plane>();
	}

	public Protein(String protein) throws Exception {
		this();
		try {
			readProtein(protein);

			cleanburied();
			makeConvexHull();
			System.out.println(hull.size());
		} catch (Exception e) {
			throw new Exception(e);
		}
	}

	private void makeConvexHull() throws Exception {
		int a;
		Plane[] planes;
		Atom3D a1, a2, a3, temp;
		double d;

		int n = 0, m = 0;

		for (int i = 0; i < atoms.size(); i++) {
			a1 = atoms.get(i);

			for (int j = i + 1; j < atoms.size(); j++) {
				a2 = atoms.get(j);

				if (!canSee(atoms, a1, a2, 3)) {
					continue;
				}

				for (int k = j + 1; k < atoms.size(); k++) {
					a3 = atoms.get(k);

					if (!canSee(atoms, a1, a3, 3))
						continue;
					if (!canSee(atoms, a2, a3, 3))
						continue;

					planes = PlaneMath.getPlanes(a1, a2, a3);

					if (planes == null) {
						continue;
					}

					n += 2;

					for (int p = 0; p < planes.length; p++) {
						if (planes[p] == null)
							continue;

						for (a = 0; a < atoms.size(); a++) {
							temp = atoms.get(a);
							if (planes[p].a1 == temp || planes[p].a2 == temp
									|| planes[p].a3 == temp) {
								continue;
							}

							d = planes[p].dis(temp.x, temp.y, temp.z);

							if (d <= temp.r - Util.ERROR * 1000) {
								if (d > 0) {
									// System.out.println("Bad Mojo  " + d + " "
									// + temp + " >< "
									// + (Math.abs(temp.r - d)));
								}
								planes[p] = null;
								m++;

								break;
							}
						}
					}

					for (Plane p : planes) {
						if (p != null) {
							hull.add(p);
							// p.createPlane(Color.blue);
							// this.addChild(p);
						}
					}
				}
			}
		}

		System.out.println(n + " <> " + m);
	}

	private void cleanburied() {
		Atom3D query, temp;
		final double extraRadius = 1.5;

		for (int i = 0; i < atoms.size(); i++) {
			query = atoms.get(i);

			for (int j = 0; j < atoms.size(); j++) {
				if (j == i)
					continue;

				temp = atoms.get(j);

				if (temp.dis3D(query) <= temp.r + query.r + extraRadius * 2)
					query.neigh.add(temp);
			}
		}

		final double DX = 0.25;
		double myDphi, nPhi, phi, r, myDtheta, nTheta, theta, dSQ;
		double x, y, z;
		boolean buried;

		for (int a = 0; a < atoms.size(); a++) {
			query = atoms.get(a);

			myDphi = DX / (query.r + extraRadius);
			nPhi = (int) (Math.PI / myDphi) + 1;
			myDphi = Math.PI / (double) (nPhi);

			query.buried = true;

			start: for (int i = 0; i < nPhi; i++) {
				phi = (i + 0.5) * myDphi;
				r = (query.r + extraRadius) * Math.sin(phi);
				myDtheta = DX / r;
				nTheta = (int) (2 * Math.PI / myDtheta) + 1;
				myDtheta = 2.0 * Math.PI / (double) (nTheta);

				for (int j = 0; j < nTheta; j++) {
					theta = myDtheta * j;
					x = query.x + (query.r + extraRadius) * Math.sin(phi)
							* Math.cos(theta);
					y = query.y + (query.r + extraRadius) * Math.sin(phi)
							* Math.sin(theta);
					z = query.z + (query.r + extraRadius) * Math.cos(phi);

					buried = false;
					for (int b = 0; b < query.neigh.size(); b++) {
						temp = query.neigh.get(b);

						if (temp == query)
							continue;

						dSQ = temp.disSQ(x, y, z);

						if (dSQ < (temp.r + extraRadius)
								* (temp.r + extraRadius)) {
							buried = true;
							break;
						}
					}

					if (!buried) {
						query.buried = false;
						break start;
					}
				}
			}
		}

		for (int i = 0; i < atoms.size(); i++) {
			query = atoms.get(i);
			query.neigh = null;

			if (query.buried) {
				atoms.remove(i);
				i--;
			}
		}
	}

	private void readProtein(String protein) throws Exception {
		Scanner sc = new Scanner(new File(protein));
		sc.nextLine();
		Atom3D a;

		while (sc.hasNextLine()) {
			a = Atom3D.readAtom(sc.nextLine());

			atoms.add(a);
		}

		sc.close();
	}

	public Atom3D curvature(double x, double y, double z) {
		Atom3D query = new Atom3D(x, y, z, 0);

		for (Plane p : hull) {
			if (p.dis(x, y, z) <= 0)
				return new Atom3D(x, y, z, Double.POSITIVE_INFINITY);
		}

		PriorityQueue<Atom3D> pq = new PriorityQueue<Atom3D>();
		double d;

		for (Atom3D a : atoms) {
			d = a.getDis3D(x, y, z);

			if (d <= 20) {
				pq.add(a, d);
			}
		}

		ArrayList<Atom3D> atms = new ArrayList<Atom3D>();
		Atom3D a3d;

		while (pq.hasNext()) {
			a3d = pq.dequeue();

			if (this.canSee(atms, a3d, query, 1))
				atms.add(a3d);
		}

		double max = -1;
		Atom3D[] spheres;
		Atom3D a1, a2, a3, a4;
		Atom3D maxSph = null;
		a4 = new Atom3D(x, y, z, 0);
		Triple maxTrip = null;

		for (int i = 0; i < atms.size(); i++) {
			a1 = atms.get(i);
			for (int j = i + 1; j < atms.size(); j++) {
				a2 = atms.get(j);

				if (a2.dis3D(a1) > 20)
					continue;

				for (int k = j + 1; k < atms.size(); k++) {
					a3 = atms.get(k);

					if (a2.dis3D(a3) > 20 || a1.dis3D(a3) > 20)
						continue;

					spheres = Curvature3D.intersectingSphere(a1, a2, a3, a4);

					if (spheres == null)
						continue;

					for (Atom3D sph : spheres) {
						if (sph == null)
							continue;

						if (sph.r > max) {
							if (isGood(sph, atms)) {
								max = sph.r;

								maxTrip = new Triple(a1, a2, a3);
								maxSph = sph;
							} else {

							}
						}
					}
				}
			}
		}

		if (maxTrip == null) {
			System.err.println("ERROR:\n" + x + " <> " + y + " <> " + z);

			return new Atom3D(x, y, z, findClosest(x, y, z));
		}

		for (int i = 0; i < atms.size(); i++) {
			a4 = atms.get(i);

			if (a4 == maxTrip.a1 || a4 == maxTrip.a2 || a4 == maxTrip.a3)
				continue;

			spheres = Curvature3D.intersectingSphere(maxTrip.a1, maxTrip.a2,
					maxTrip.a3, a4);

			if (spheres == null)
				continue;

			for (Atom3D sph : spheres) {
				if (sph == null)
					continue;

				if (sph.r > max) {
					if (isGood(sph, atms)) {
						max = sph.r;
						maxSph = sph;
					}
				}
			}
		}

		return maxSph;
	}

	private double findClosest(double x, double y, double z) {
		double dist = Double.POSITIVE_INFINITY;
		double distTemp;

		for (int i = 0; i < atoms.size(); i++) {
			distTemp = atoms.get(i).getDis3D(x, y, z) - atoms.get(i).r;

			dist = Math.min(dist, distTemp);
		}

		return dist;
	}

	public boolean canSee(ArrayList<Atom3D> atms, Atom3D a1, Atom3D a2, int n) {
		Vector3d quest = a1.getVector();
		Vector3d pQuery = a2.getVector();
		Vector3d diff = new Vector3d();
		Vector3d norm = new Vector3d();
		double t, disSq, r;
		Atom3D a3D;
		Vector3d temp;
		Vector3d test = new Vector3d();

		diff.sub(pQuery, quest);
		norm.normalize(diff);

		r = diff.dot(norm);

		int a = 0;

		for (int i = 0; i < atms.size(); i++) {
			a3D = atms.get(i);

			if (a3D == a1 || a3D == a2)
				continue;

			temp = a3D.getVector();
			test.sub(temp, quest);
			t = test.dot(norm);

			if (t < 0 || t > r)
				continue;

			disSq = test.lengthSquared() - t * t;

			if (disSq < a3D.r * a3D.r) {
				a++;

				if (a >= n)
					return false;
			}
		}

		return true;
	}

	private boolean isGood(Atom3D sph, ArrayList<Atom3D> atms) {
		for (Atom3D a : atms) {
			if (a.dis3D(sph) - (a.r + sph.r) < -Util.ERROR) {
				// System.out.println(a + "\n" + sph);
				// System.out.println((a.dis3D(sph) - (a.r + sph.r)) + "\n");
				return false;
			}
		}

		return true;
	}
}
