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
     * Copy constructor
     */
    public Protein(Protein o) {
        this();
        for (int i = 0; i < o.atoms.size(); i++) {
            atoms.add(new Atom3D(o.atoms.get(i)));
        }

        for (int i = 0; i < o.hull.size(); i++) {
            hull.add(new Plane(o.hull.get(i)));
        }
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
        readProtein(protein);

        System.out.println("Atoms: " + atoms.size());
    }

    // Clean the protein to speed up curvature
    public void clean() {
        cleanburied();
        makeConvexHull();
        System.out.println("Hull: " + hull.size() + " <> Atoms: "
                + atoms.size());
    }

    /**
     * Creates the convex hull by using a semi-brute force method.
     * 
     * @throws Exception
     *             if there is a problem calculating a plane based on a singular
     *             matrix.
     */
    private void makeConvexHull() {
        int a;
        Plane[] planes;
        Atom3D a1, a2, a3, temp;
        double d;

        for (int i = 0; i < atoms.size(); i++) {
            a1 = atoms.get(i);

            for (int j = i + 1; j < atoms.size(); j++) {
                a2 = atoms.get(j);

                for (int k = j + 1; k < atoms.size(); k++) {
                    a3 = atoms.get(k);

                    planes = PlaneMath.getPlanes(a1, a2, a3);

                    if (planes == null) {
                        continue;
                    }

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
                                planes[p] = null;

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

                theta: for (int j = 0; j < nTheta; j++) {
                    theta = myDtheta * j;
                    x = query.x + (query.r + extraRadius) * Math.sin(phi)
                            * Math.cos(theta);
                    y = query.y + (query.r + extraRadius) * Math.sin(phi)
                            * Math.sin(theta);
                    z = query.z + (query.r + extraRadius) * Math.cos(phi);

                    for (int b = 0; b < query.neigh.size(); b++) {
                        temp = query.neigh.get(b);

                        if (temp == query)
                            continue;

                        dSQ = temp.disSQ(x, y, z);

                        if (dSQ < (temp.r + extraRadius)
                                * (temp.r + extraRadius)) {
                            continue theta;
                        }
                    }

                    query.buried = false;
                    break start;
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

                for (int k = j + 1; k < atms.size(); k++) {
                    a3 = atms.get(k);

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
                            }
                        }
                    }
                }
            }
        }

        if (maxTrip == null) {
            System.err.println("ERROR:\n" + x + " <> " + y + " <> " + z);

            // atms.get(0) should be the closest atom
            // return new Atom3D(x, y, z, findClosest(x, y, z));

            Atom3D a = findPair(x, y, z, atms);
            return new Atom3D(x, y, z, a.r);
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

    // Finds the largest sphere when the point is between two atoms
    // Has 3 Tests
    // 1. Point is within convex hull
    // 2. Is the point linear with the two atoms
    // 3. When it is not
    public Atom3D findPair(double x, double y, double z, ArrayList<Atom3D> atms) {
        Atom3D best = new Atom3D(x, y, z, 0);
        Atom3D q = new Atom3D(x, y, z, 0);

        for (int i = 0; i < atms.size(); i++) {
            for (int j = i + 1; j < atms.size(); j++) {
                // Check if the point is not in the convex hull of the two atoms
                if (!inAtoms(atms.get(i), atms.get(j), q)) {
                    continue;
                }

                // Are the atoms in a singular line
                if (isSingular(atms.get(i), atms.get(j), q)) {
                    // Solve for the singular

                    System.err.println("Expected a linear set for the atoms");

                    // TODO Check to see if there is a valid sphere
                } else {
                    Atom3D[] spheres = Curvature2D.nonSingular(atms.get(i),
                            atms.get(j), q);

                    check: for (Atom3D a : spheres) {
                        // Check if the sphere is valid
                        // Check if the radius is larger

                        if (a.r <= best.r) {
                            continue;
                        }

                        for (int k = 0; k < atms.size(); k++) {
                            if (a.dis3D(atms.get(k)) < a.r + atms.get(k).r)
                                continue check;
                        }

                        best = a;
                    }
                }
            }
        }

        return best;
    }

    // Checks if the atom is in the convex hull of the other two atoms
    private boolean inAtoms(Atom3D a1, Atom3D a2, Atom3D q) {
        Vector3d v1 = a1.getVector();
        Vector3d v2 = a2.getVector();
        Vector3d vq = q.getVector();

        Vector3d line = new Vector3d();
        line.sub(v2, v1);

        Vector3d qline = new Vector3d();
        qline.sub(vq, v1);

        double dist = line.dot(qline);

        if (dist < 0 || dist > 1)
            return false;

        double rq = Math.sqrt(qline.dot(qline) - dist * dist);

        double rp = (a2.r - a1.r) * dist + a1.r;

        return rq <= rp;
    }

    // Check if the two atoms and point is in a straight line
    private boolean isSingular(Atom3D a1, Atom3D a2, Atom3D q) {
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();

        v1.sub(a1.getVector(), q.getVector());
        v2.sub(a2.getVector(), q.getVector());

        double d1 = Math.abs(v1.dot(v2));
        double d2 = Math.abs(v2.length() * v1.length());

        return d1 == d2;
    }
}
