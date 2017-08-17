import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;

/**
 * The main entry point into the program. Can change according to what it is
 * needed to do.
 * 
 * @author Kyle Diller
 *
 */
public class Main {

    /**
     * This version of the program reads in a protein and a molecule, then
     * computes the curvature of the protein at each of the molecules atoms.
     * 
     * @param args
     *            the list of files to read and write to.
     * @throws Exception
     *             thrown if there is some problem reading a file.
     */
    public static void main(String[] args) throws Exception {
        if (args.length != 3) {
            usage();
            return;
        }
        
        System.out.println("Starting Program");

        long start = System.currentTimeMillis();
        Protein p = new Protein(args[0]);
        System.out
                .println("Make Time: " + (System.currentTimeMillis() - start));

        Scanner sc = new Scanner(new File(args[1]));
        FileWriter os = new FileWriter(new File(args[2]));

        sc.nextLine();
        os.write("Point,X,Y,Z,R");

        String line;
        Atom3D temp;
        double x, y, z;
        String[] split;

        long tot = 0;
        int i = 0;

        while (sc.hasNextLine()) {
            line = sc.nextLine();
            split = line.split(",");
            x = Double.parseDouble(split[1]);
            y = Double.parseDouble(split[2]);
            z = Double.parseDouble(split[3]);

            start = System.currentTimeMillis();
            temp = p.curvature(x, y, z);
            tot += System.currentTimeMillis() - start;
            i++;

            os.write("\n" + split[0] + "," + temp.x + "," + temp.y + ","
                    + temp.z + "," + temp.r);
        }

        System.out.println("Avg Time: " + (tot / i));

        sc.close();
        os.close();
    }

    /**
     * A standard error usage message telling the user how to run this program.
     */
    private static void usage() {
        System.err
                .println("java -jar ProteinCurve.jar <protein file> <mol file> <output file>");
        System.err.println("Note: all files are to be in csv format");
    }
}
