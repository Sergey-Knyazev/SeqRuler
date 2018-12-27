package TN93;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;

import static java.lang.Math.log;

public class TN93 {
    public static void tn93Fasta(File inputFile, File outputFile, float edgeThreshold) {
        PrintWriter f = null;
        try {
            LinkedList<Seq> seqs = read_fasta(inputFile);
            double[][] dist = tn93(seqs);
            f = new PrintWriter(outputFile);
            f.println("Source,Target,Dist");
            for (int i = 1; i < dist.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    if (dist[i][j] <= edgeThreshold) f.println(
                            String.format("%s,%s,%f", seqs.get(i).getName(), seqs.get(j).getName(), dist[i][j]));
                }
            }
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }

    public static double[][] tn93(LinkedList<Seq> seqs) {
        double[][] dist = new double[seqs.size()][seqs.size()];
        long startTime = System.nanoTime(), estimatedTime;
        for (int i = 1; i < dist.length; ++i) {
            if(i % 1000 == 0) {
                estimatedTime = System.nanoTime() - startTime;
                System.out.print(String.format("%d out of %d for ", i, dist.length));
                System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                System.out.println(" sec");
            }
            for (int j = 0; j < i; ++j) {
                dist[i][j] = dist[j][i] = tn93(seqs.get(i).getSeq_enc(), seqs.get(j).getSeq_enc());
            }
        }
        return dist;
    }

    private static double tn93(int[] s1, int[] s2) {
        int total_non_gap = 0;
        int[] nucl_counts = new int[4];
        int[][] nucl_pair_count = new int[4][4];

        int length = Math.min(s1.length, s2.length);

        for(int i=0; i<length; ++i) {
            int c1 = s1[i];
            int c2 = s2[i];
            if(c1==Seq.N || c2==Seq.N) continue;
            ++nucl_counts[c1];
            ++nucl_counts[c2];
            ++nucl_pair_count[c1][c2];
            ++nucl_pair_count[c2][c1];
            ++total_non_gap;
        }
        double[] nucl_freq = new double[4];
        for(int i=0; i<4; ++i) nucl_freq[i] = (double) nucl_counts[i]/2/total_non_gap;
        double p1 = (double) nucl_pair_count[Seq.A][Seq.G]/total_non_gap;
        double p2 = (double) nucl_pair_count[Seq.C][Seq.T]/total_non_gap;
        double q = ((double) nucl_pair_count[Seq.A][Seq.C]+nucl_pair_count[Seq.A][Seq.T]+nucl_pair_count[Seq.C][Seq.G]+
                nucl_pair_count[Seq.G][Seq.T])/total_non_gap;
        double g_r = nucl_freq[Seq.A]+nucl_freq[Seq.G];
        double g_y = nucl_freq[Seq.C]+nucl_freq[Seq.T];
        boolean useK2P = false;
        //for(int i=0; i<4; ++i) if (nucl_freq[i] == 0) useK2P = true;
        //if (useK2P) {
        //    //TODO:
        //    System.out.println("Hi");
        //    return 0;
        //}
        double k_ag = nucl_freq[Seq.A]*nucl_freq[Seq.G]/g_r;
        double k_tc = nucl_freq[Seq.T]*nucl_freq[Seq.C]/g_y;
        double k_ry = g_r*g_y;

        return -2*k_ag*log(1-p1/(2*k_ag)-q/(2*g_r))
                -2*k_tc*log(1-p2/(2*k_tc)-q/(2*g_y))
                -2*(k_ry-k_ag*g_y-k_tc*g_r)*log(1-q/(2*k_ry));
    }

    public static LinkedList<Seq> read_seqs(Scanner sc) {
        LinkedList<Seq> seqs = new LinkedList<Seq>();
        String name="", seq="";
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            if(line.charAt(0)=='>') {
                if (name.length()!=0) seqs.add(new Seq(name, seq));
                name = line.substring(1);
                seq="";
            }
            else seq=seq.concat(line);
        }
        if(name.length()!=0) seqs.add(new Seq(name, seq));
        return seqs;
    }

    private static LinkedList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        LinkedList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
}
