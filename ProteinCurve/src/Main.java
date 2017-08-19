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
    public static void main(String[] args) throws Exception {
        Main m = new Main();
        m.run();
    }

    public synchronized void run() throws Exception {
        File file = new File("Curvature");
        String[] names = file.list();
        String last = "";
        Node<String> pros = null;

        final int[] using = { 0, 0 };

        for (int i = 0; i < names.length; i++) {
            String cur = names[i].substring(0, names[i].indexOf('-'));
            if (!last.equals(cur)) {
                pros = new Node<String>(cur, pros);
                last = cur;
                using[1]++;
            }
        }
        final int max = 4;

        while (pros != null) {
            while (using[0] >= max)
                wait();

            final String name = pros.getData();

            Thread t = new Thread(() -> {
                System.err.println("Starting " + name + " <> " + using[1]);

                try {
                    main1(new String[] { "Curvature\\" + name + "-sphere.csv",
                            "Curvature\\" + name + "-ligand.csv",
                            "Curvature\\" + name + "-out.csv" });
                } catch (Exception e) {
                }
                using[0]--;
                returnCore();
            });
            t.start();
            using[0]++;
            using[1]--;

            pros = pros.getNext();
        }

        System.err.println("FINISHED WITH EVERYTHING");
    }
    
    public synchronized void returnCore(){
        notifyAll();
    }

    /**
     * This version of the program reads in a protein and a molecule, then
     * computes the curvature of the protein at each of the molecules atoms.
     * 
     * @param args
     *            the list of files to read and write to.
     * @throws Exception
     *             thrown if there is some problem reading a file.
     */
    public static void main1(String[] args) throws Exception {
        if (args.length != 3) {
            usage();
            return;
        }

        System.out.println("Starting Program");

        Protein p = new Protein(args[0]);

        System.out.println("Cleaning");
        long start = System.currentTimeMillis();
        p.clean();
        System.out
                .println("Make Time: " + (System.currentTimeMillis() - start));

        Scanner sc = new Scanner(new File(args[1]));
        FileWriter os = new FileWriter(new File(args[2]));

        sc.nextLine();
        os.write("Point,X,Y,Z,R Approx,R Num,Time");

        String line;
        Atom3D temp;
        double x, y, z;
        String[] split;

        long tot = 0;
        int i = 0;

        System.out.println("Calculating Curvature");

        while (sc.hasNextLine()) {
            line = sc.nextLine();
            split = line.split(",");
            x = Double.parseDouble(split[3]);
            y = Double.parseDouble(split[4]);
            z = Double.parseDouble(split[5]);

            start = System.currentTimeMillis();
            temp = p.curvature(x, y, z);
            long time = System.currentTimeMillis() - start;
            tot += time;
            i++;

            if (i % 100 == 0) {
                System.out.println(i + " <> " + (tot / i));
            }

            os.write("\n" + split[0] + "," + temp.x + "," + temp.y + ","
                    + temp.z + "," + temp.r + "," + split[7] + "," + time);
            os.flush();
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
