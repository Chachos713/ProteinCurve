public class Util {
	public static final double ERROR = 0.001;
	
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
