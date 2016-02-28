import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;

public class Main {

	public static void main(String[] args) throws Exception {
		if (args.length != 3) {
			usage();
			return;
		}

		Protein p = new Protein(args[0]);

		Scanner sc = new Scanner(new File(args[1]));
		FileWriter os = new FileWriter(new File(args[2]));

		sc.nextLine();
		os.write("Point,X,Y,Z,R");

		String line;
		Atom3D temp;
		double x, y, z;
		String[] split;

		while (sc.hasNextLine()) {
			line = sc.nextLine();
			split = line.split(",");
			x = Double.parseDouble(split[1]);
			y = Double.parseDouble(split[2]);
			z = Double.parseDouble(split[3]);

			temp = p.curvature(x, y, z);

			os.write("\n" + split[0] + "," + temp.x + "," + temp.y + ","
					+ temp.z + "," + temp.z);
		}

		sc.close();
		os.close();
	}

	private static void usage() {
		System.err
				.println("java -jar ProteinCurve.jar <protein file> <mol file> <output file>");
		System.err.println("Note: all files are to be in csv format");
	}
}
