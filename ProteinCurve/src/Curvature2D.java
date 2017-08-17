import javax.vecmath.Vector3d;

import Jama.Matrix;

public class Curvature2D {

    public static Atom3D[] nonSingular(Atom3D a1, Atom3D a2, Atom3D q) {
        // First create the plane that the atoms lie on
        // q becomes the origin
        Vector3d X = new Vector3d();
        X.sub(a1.getVector(), q.getVector());
        X.normalize();

        Vector3d temp = new Vector3d(X);

        Vector3d Y = new Vector3d();
        Y.sub(a2.getVector(), q.getVector());
        temp.scale(X.dot(Y));
        Y.sub(Y, temp);
        Y.normalize();

        temp.sub(a1.getVector(), q.getVector());

        // X and Y axes are created with q being the origin
        Atom3D a12d = new Atom3D(temp.dot(X), temp.dot(Y), 0, a1.r);

        temp.sub(a2.getVector(), q.getVector());
        Atom3D a22d = new Atom3D(temp.dot(X), temp.dot(Y), 0, a2.r);

        // Solve for the radius
        Matrix A = new Matrix(2, 2);
        A.set(0, 0, a12d.x);
        A.set(0, 1, a12d.y);
        A.set(1, 0, a22d.x);
        A.set(1, 1, a22d.y);
        Matrix AI = A.inverse();

        Matrix C = new Matrix(2, 1);
        C.set(0, 0, a12d.r);
        C.set(1, 0, a22d.r);

        Matrix D = new Matrix(2, 1);
        D.set(0, 0, .5 * (a12d.getVector().dot(a12d.getVector()) - a12d.r
                * a12d.r));
        D.set(1, 0, .5 * (a22d.getVector().dot(a22d.getVector()) - a22d.r
                * a22d.r));

        // Found the matrices for A, B, and C
        // Using the quarratic equation, it can be found the solution to the
        // radius

        Matrix AIC = AI.times(C);
        Matrix AID = AI.times(D);

        double a = AIC.get(0, 0) * AIC.get(0, 0) + AIC.get(1, 0)
                * AIC.get(1, 0);
        double b = -2
                * (AID.get(0, 0) * AIC.get(0, 0) + AID.get(1, 0)
                        * AIC.get(1, 0));
        double c = AID.get(0, 0) * AID.get(0, 0) + AID.get(1, 0)
                * AID.get(1, 0);

        double[] sols = Util.quadratic(a, b, c);

        Atom3D[] atms = new Atom3D[sols.length];

        for (int i = 0; i < sols.length; i++) {
            Matrix P = AI.times(D.minus(C.times(sols[i])));
            double x = P.get(0, 0);
            double y = P.get(1, 0);
            atms[i] = new Atom3D(q.x + X.x * x + Y.x * y, q.y + X.y * y + Y.y
                    * y, q.z + X.z * x + Y.z * y, sols[i]);

            if (!check(atms[i], a1, a2, q))
                atms[i] = null;
        }
        return atms;
    }

    public static boolean check(Atom3D sph, Atom3D a1, Atom3D a2, Atom3D q) {
        if (sph.dis3D(a1) < a1.r + sph.r) {
            System.err.println("ERROR ON SPHERE");
            return false;
        }

        if (sph.dis3D(a2) < a2.r + sph.r) {
            System.err.println("ERROR ON SPHERE");
            return false;
        }

        if (sph.dis3D(q) < q.r + sph.r) {
            System.err.println("ERROR ON SPHERE");
            return false;
        }

        return true;
    }
}
