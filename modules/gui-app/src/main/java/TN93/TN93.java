package TN93;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.Observable;

import static java.lang.Math.log;

public class TN93 extends Observable {
    private static final double TN_93_MAX_DIST = 1000.;
    private float edgeThreshold = 1;
    private File inputFile;
    private File outputFile;

    public void setEdgeThreshold(float edgeThreshold) {
        this.edgeThreshold = edgeThreshold;
    }


    public void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(File outputFile) {
        this.outputFile = outputFile;
    }

    public void tn93Fasta() {
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

    public double[][] tn93(LinkedList<Seq> seqs) {
        double[][] dist = new double[seqs.size()][seqs.size()];
        long startTime = System.nanoTime(), estimatedTime;
        int pairs_count = (dist.length * dist.length - dist.length)/2;
        int current_pair = 0;
        for (int i = 1; i < dist.length; ++i) {
            for (int j = 0; j < i; ++j) {
                current_pair++;
                if(current_pair % (pairs_count/100) == 0) {
                    estimatedTime = System.nanoTime() - startTime;
                    int percCompleted = current_pair*100/pairs_count;
                    System.out.print(String.format("%d%% completed for ", percCompleted));
                    System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                    System.out.println(" sec");
                    setChanged();
                    notifyObservers(percCompleted);
                }
                dist[i][j] = dist[j][i] = tn93(seqs.get(i).getSeq_enc(), seqs.get(j).getSeq_enc());
            }
        }
        setChanged();
        notifyObservers(100);
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
        for(int i=0; i<4; ++i) if (nucl_freq[i] == 0) useK2P = true;

        double k_ag = nucl_freq[Seq.A]*nucl_freq[Seq.G]/g_r;
        double k_tc = nucl_freq[Seq.T]*nucl_freq[Seq.C]/g_y;
        double k_ry = g_r*g_y;

        double dist;

        if(useK2P) {
            double l1 = 1.-2.*(p1+p2)-q;
            double l2 = 1.-2.*q;
            if(l1>0.&&l2>0.) dist = log(p1)/2-log(p2)/4;
            else dist=TN_93_MAX_DIST;
        }
        else {
            double l_ag = 1-p1/(2*k_ag)-q/(2*g_r);
            double l_tc = 1-p2/(2*k_tc)-q/(2*g_y);
            double l_ry = 1-q/(2*k_ry);
            if(l_ag>0. && l_tc>0. && l_ry>0.)
                dist=-2*k_ag*log(l_ag)-2*k_tc*log(l_tc)-2*(k_ry-k_ag*g_y-k_tc*g_r)*log(l_ry);
            else dist=TN_93_MAX_DIST;
        }
        return dist <= 0. ? 0. : dist;
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
