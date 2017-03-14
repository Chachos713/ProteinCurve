import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

import javax.vecmath.Vector3d;

/**
 * A storage class to read, store, and compute the curvature of the protein.
 * 
 * @author Kyle Diller
 *
 */
public class Protein {
    /**
     * The list of atoms within the protein.
     */
    private ArrayList<Atom3D> atoms;

    /**
     * The list of planes that make up the convex hull.
     */
    private ArrayList<Plane> hull;

    /**
     * Creates an empty protein with no atoms or a convex hull.
     */
    private Protein() {
        atoms = new ArrayList<Atom3D>();
        hull = new ArrayList<Plane>();
    }

    /**
     * Creates a protein from a file (protein) that is read.
     * 
     * @param protein
     *            the file to read the protein from.
     * @throws Exception
     *             if there is a problem with the file.
     */
    public Protein(String protein) throws Exception {
        this();
        try {
            readProtein(protein);

            cleanburied();
            makeConvexHull();
            System.out.println(hull.size() + " <> " + atoms.size());
            System.out.println(hull.get(0));
        } catch (Exception e) {
            throw new Exception(e);
        }
    }

    /**
     * Creates the convex hull by using a semi-brute force method.
     * 
     * @throws Exception
     *             if there is a problem calculating a plane based on a singular
     *             matrix.
     */
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
                        }
                    }
                }
            }
        }

        System.out.println(n + " <> " + m);
    }

    /**
     * Cleans the atoms that are buried within the protein and have little to no
     * chance of being part of the curvature of the protein.
     */
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

    /**
     * Reads the file (protein) and stores the atoms that are read.
     * 
     * @param protein
     *            the file to read.
     * @throws Exception
     *             if there are problems with reading the file.
     */
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

    /**
     * Computes the curvature of the protein at a given point. <br>
     * Uses a semi-bruteforce algorithm to compute the curvature.
     * 
     * @param x
     *            the x component of the point.
     * @param y
     *            the y component of the point.
     * @param z
     *            the z component of the point.
     * @return the sphere (as an atom) that has the largest radius and contains
     *         the point.
     */
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

            // Checks if the point is in an atom
            if (d < a.r)
                return null;

            pq.add(a, d);
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

                if (sph.getDis3D(x, y, z) > sph.r)
                    continue;

                if (sph.r > max) {
                    if (isGood(sph, atms)) {
                        max = sph.r;
                        maxSph = sph;
                    }
                }
            }
        }

        return new Atom3D(x, y, z, maxSph.r);
    }

    /**
     * Finds the closest atom to a given point. <br>
     * Used if there is a problem computing the curvature of the atom.
     * 
     * @param x
     *            the x component of the point.
     * @param y
     *            the y component of the point.
     * @param z
     *            the z component of the point.
     * @return the distance between the nearest atom and the point.
     */
    private double findClosest(double x, double y, double z) {
        double dist = Double.POSITIVE_INFINITY;
        double distTemp;

        for (int i = 0; i < atoms.size(); i++) {
            distTemp = atoms.get(i).getDis3D(x, y, z) - atoms.get(i).r;

            dist = Math.min(dist, distTemp);
        }

        return dist;
    }

    /**
     * Counts how many atoms lie between two atoms.
     * 
     * @param atms
     *            the list of atoms to check against.
     * @param a1
     *            the first atom.
     * @param a2
     *            the second atom.
     * @param n
     *            the threshold of how many atoms make it invisible.
     * @return true if the number of atoms < n.
     */
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

    /**
     * Checks if a sphere does not enter the protein.
     * 
     * @param sph
     *            the sphere to check.
     * @param atms
     *            the list of atoms to check against.
     * @return true if the sphere does not enter the protein.
     */
    private boolean isGood(Atom3D sph, ArrayList<Atom3D> atms) {
        for (Atom3D a : atms) {
            if (a.dis3D(sph) - (a.r + sph.r) < -Util.ERROR) {
                return false;
            }
        }

        return true;
    }

    /**
     * Checks if a point is inside the protein or not
     * 
     * @param x
     *            the x coordinate of the point to check
     * @param y
     *            the y coordinate of the point to check
     * @param z
     *            the z coordinate of the point to check
     * @return true if the point is within the protein
     */
    public boolean isInside(double x, double y, double z) {
        for (Atom3D a : atoms) {
            if (a.getDis3D(x, y, z) <= a.r)
                return true;
        }

        return false;
    }
}
