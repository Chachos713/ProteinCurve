/**
 * A standard utility class with some values and functions that are needed
 * through out the program.
 * 
 * @author Kyle Diller
 *
 */
public class Util {
	/**
	 * The standard +/- range for an error to be considered acceptable. This is
	 * used due to the fact of rounding errors when dealing with float and
	 * double arithmetic.
	 */
	public static final double ERROR = 0.001;

	/**
	 * Computes the two possible values for x based on the quadratic formula.
	 * 
	 * @param a
	 *            the coefficient to the squared portion of the quadratic
	 *            equation.
	 * @param b
	 *            the coefficient to the linear portion of the quadratic
	 *            equation.
	 * @param c
	 *            the constant at the end of the quadratic equation.
	 * @return The possible values of x in the equation y = a * x^2 + b * x + c. <br>
	 *         null - if the x values are imaginary. <br>
	 *         double[] - if the x values are real.
	 */
	public static double[] quadratic(double a, double b, double c) {
		if (b * b - 4 * a * c < 0)
			return null;
		else if (b * b - 4 * a * c == 0) {
			double[] num = new double[1];
			num[0] = -b / (2 * a);
			return num;
		} else {
			double rad = Math.sqrt(b * b - 4 * a * c);
			double[] num = new double[2];
			num[0] = (-b + rad) / (2 * a);
			num[1] = (-b - rad) / (2 * a);

			return num;
		}
	}
}
