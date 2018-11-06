package TN93;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static java.lang.Math.log;

public class TN93 {
    public static void tn93Fasta(File inputFile, File outputFile) {
        try {
            LinkedList<Seq> seqs = read_fasta(inputFile);
            double[][] dist = tn93(seqs);
            PrintWriter f = new PrintWriter(outputFile);
            for (int i = 1; i < dist.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    f.println(String.format("%s,%s,%f", seqs.get(i).name, seqs.get(j).name, dist[i][j]));
                }
            }
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static double[][] tn93(LinkedList<Seq> seqs) {
        double[][] dist = new double[seqs.size()][seqs.size()];
        for (int i = 1; i < dist.length; ++i) {
            for (int j = 0; j < i; ++j) {
                dist[i][j] = dist[j][i] = tn93(seqs.get(i).seq, seqs.get(j).seq);
            }
        }
        return dist;
    }

    static double tn93(String s1, String s2) {
        int A=0,C=1,G=2,T=3;

        Map<Character, Integer> nucl = new HashMap<Character, Integer>();
        nucl.put('A', A);
        nucl.put('C', C);
        nucl.put('G', G);
        nucl.put('T', T);

        List<Integer> s1_enc = new LinkedList<Integer>();
        List<Integer> s2_enc = new LinkedList<Integer>();
        for(int i=0; i<s1.length(); ++i) {
            char c1 = s1.charAt(i);
            char c2 = s2.charAt(i);
            if(!(nucl.containsKey(c1) && nucl.containsKey(c2))) continue;
            s1_enc.add(nucl.get(c1));
            s2_enc.add(nucl.get(c2));
        }

        int n = 0;
        int[] nucl_counts = new int[4];
        int[][] nucl_pair_count = new int[4][4];

        for(int i=0; i<s1_enc.size(); ++i) {
            int c1 = s1_enc.get(i);
            int c2 = s2_enc.get(i);
            ++nucl_counts[c1];
            ++nucl_counts[c2];
            ++nucl_pair_count[c1][c2];
            ++nucl_pair_count[c2][c1];
            ++n;
        }
        double[] nucl_freq = new double[4];
        for(int i=0; i<4; ++i) nucl_freq[i] = (double) nucl_counts[i]/2/n;
        double p1 = (double) nucl_pair_count[A][G]/n;
        double p2 = (double) nucl_pair_count[C][T]/n;
        double q = ((double) nucl_pair_count[A][C]+nucl_pair_count[A][T]+nucl_pair_count[C][G]+
                nucl_pair_count[C][T]+nucl_pair_count[G][T])/n;
        double g_r = nucl_freq[A]+nucl_freq[G];
        double g_y = nucl_freq[A]+nucl_freq[G];
        double k_ag = 2*nucl_freq[A]*nucl_freq[G]/g_r;
        double k_tc = 2*nucl_freq[T]*nucl_freq[C]/g_y;
        double k_ry = 2*g_r*g_y;

        return -k_ag*log(1-p1/k_ag-q/(2*g_r))
                -k_tc*log(1-p2/k_tc-q/(2*g_y))
                -(k_ry-k_ag*g_y-k_tc*g_r)*log(1-q/k_ry);
    }

    public static LinkedList<Seq> read_seqs(Scanner sc) {
        LinkedList<Seq> seqs = new LinkedList<Seq>();
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if (line.charAt(0) == '>') {
                seqs.add(new Seq(line.substring(1), ""));
            }
            else {
                seqs.getLast().seq = seqs.getLast().seq.concat(line);
            }
        }
        return seqs;
    }

    private static LinkedList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        return read_seqs(sc);
    }
}
