import javax.vecmath.Vector3d;

/**
 * A utility class used to calculate the planes defined by 3 atoms. Mainly used
 * to calculate the convex hull of the protein.
 * 
 * @author Kyle Diller
 *
 */
public class PlaneMath {

    /**
     * Calculates the 2 planes that are defined by 3 atoms.
     * 
     * @param a1
     *            The first atom.
     * @param a2
     *            The second atom.
     * @param a3
     *            The third atom.
     * @return the planes that are defined by the 3 atoms.
     */
    public static Plane[] getPlanes(Atom3D a1, Atom3D a2, Atom3D a3) {
        Vector3d p1, p2, p3;
        p1 = a1.getVector();
        p2 = a2.getVector();
        p3 = a3.getVector();
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();
        Vector3d v3 = new Vector3d();
        Vector3d temp = new Vector3d();

        v1.sub(p2, p1);
        temp.sub(p3, p1);
        double check1 = Math.abs(v1.dot(temp));
        double check2 = Math.abs(v1.length() * temp.length());

        // Checks that the three atoms are not in a straight line
        // check1 = check2 * cos(angle btw v1 and temp)
        if (Math.abs(check1 - check2) < Util.ERROR) {
            // System.err.println("ERROR: Singular");
            return null;
        }

        // Get the normal vector
        double dis21 = v1.length();
        v1.x /= dis21;
        v1.y /= dis21;
        v1.z /= dis21;

        if (Math.abs(1 - v1.length()) > Util.ERROR) {
            // System.err.println("ERROR: v1 Not Normalized");
            return null;
        }

        v2.sub(p3, p1);
        double dot = v2.dot(v1);
        temp = new Vector3d(v1);
        temp.x *= dot;
        temp.y *= dot;
        temp.z *= dot;

        v2.sub(temp);
        double dis31 = v2.length();
        v2.x /= dis31;
        v2.y /= dis31;
        v2.z /= dis31;

        if (Math.abs(1 - v2.length()) > Util.ERROR) {
            // System.err.println("ERROR: v2 Not Normalized");
            return null;
        }

        v3.cross(v1, v2);

        if (Math.abs(1 - v3.length()) > Util.ERROR) {
            // System.err.println("ERROR: v3 Not Normalized");
            return null;
        }

        if (Math.abs(v1.dot(v2)) > Util.ERROR
                || Math.abs(v1.dot(v3)) > Util.ERROR
                || Math.abs(v2.dot(v3)) > Util.ERROR) {
            // System.err.println("ERROR: v1, v2, or v3 is not orthoganol");
            return null;
        }

        double alpha1 = (a2.r - a1.r) / dis21;

        temp.sub(p3, p1);
        dot = temp.dot(v1);

        double alpha2 = a3.r - a1.r - dot * alpha1;
        alpha2 /= dis31;

        if ((alpha1 * alpha1 + alpha2 * alpha2) - Util.ERROR > 1) {
            // System.err.println("ERROR: alpha1 and alpha2 are greater than 1");
            return null;
        }

        double[] alpha3 = Util.quadratic(1, 0, alpha1 * alpha1 + alpha2
                * alpha2 - 1);

        Vector3d v1temp = new Vector3d(v1);
        Vector3d v2temp = new Vector3d(v2);
        Vector3d v3temp, adda, knot, addatemp;

        v1temp.x *= alpha1;
        v1temp.y *= alpha1;
        v1temp.z *= alpha1;

        v2temp.x *= alpha2;
        v2temp.y *= alpha2;
        v2temp.z *= alpha2;

        if (alpha3 == null)
            return null;

        double alp3;
        Plane[] planes = new Plane[alpha3.length];

        for (int i = 0; i < alpha3.length; i++) {
            alp3 = alpha3[i];
            v3temp = new Vector3d(v3);
            v3temp.x *= alp3;
            v3temp.y *= alp3;
            v3temp.z *= alp3;

            adda = new Vector3d();
            adda.add(v1temp, v2temp);
            adda.add(v3temp);

            addatemp = new Vector3d(adda);
            addatemp.x *= a1.r;
            addatemp.y *= a1.r;
            addatemp.z *= a1.r;

            knot = new Vector3d();
            knot.sub(p1, addatemp);

            planes[i] = new Plane(knot, adda, a1, a2, a3);

            planes[i] = checkPlane(planes[i], a1, a2, a3);
        }

        return planes;
    }

    /**
     * Checks that a plane is touching all 3 spheres.
     * 
     * @param plane
     *            the plane to check.
     * @param a1
     *            the first atom.
     * @param a2
     *            the second atom.
     * @param a3
     *            the third atom.
     * @return true if all 3 atoms are touching the plane.
     */
    private static Plane checkPlane(Plane plane, Atom3D a1, Atom3D a2, Atom3D a3) {
        double disTemp = plane.dis(a1.x, a1.y, a1.z);

        if (disTemp < 0 || Math.abs(disTemp - a1.r) > Util.ERROR) {
            // System.err.println("ERROR: plane no good1:  " + disTemp + "  "
            // + a1.r);
            return null;
        }

        disTemp = plane.dis(a2.x, a2.y, a2.z);

        if (disTemp < 0 || Math.abs(disTemp - a2.r) > Util.ERROR) {
            // System.err.println("ERROR: plane no good2" + disTemp + "  " +
            // a2.r);
            return null;
        }

        disTemp = plane.dis(a3.x, a3.y, a3.z);

        if (disTemp < 0 || Math.abs(disTemp - a3.r) > Util.ERROR) {
            // System.err.println("ERROR: plane no good3" + disTemp + "  " +
            // a3.r);
            return null;
        }

        return plane;
    }
}