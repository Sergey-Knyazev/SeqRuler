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
    public static final String[] ambiguityModeList = {"resolve", "average", "skip"};
    private float edgeThreshold = 1;
    private File inputFile;
    private File outputFile;
    private String ambiguityMode;

    public void setAmbiguityMode(String ambiguityMode) {
        this.ambiguityMode = ambiguityMode;
    }

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
        catch(Exception e) {
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
                dist[i][j] = dist[j][i] = tn93(seqs.get(i).getSeqEnc(), seqs.get(j).getSeqEnc());
            }
        }
        setChanged();
        notifyObservers(100);
        return dist;
    }
    private class PairStat {
        int total_non_gap;
        double[] nucl_counts;
        double[][] nucl_pair_count;

        PairStat() {
            total_non_gap = 0;
            nucl_counts = new double[4];
            nucl_pair_count = new double[4][4];
        }
        void updateStat(int c1, int c2) {
            ++nucl_counts[c1];
            ++nucl_counts[c2];
            ++nucl_pair_count[c1][c2];
            ++nucl_pair_count[c2][c1];
        }
        void resolveMatch(int match) {
            int count=0;
            int m = match;
            for(int i=0; i!=4; ++i) {
                if((m&1)!=0) ++count;
                m=m>>1;
            }
            double frac =1.0/count;
            for(int i=0; i!=4; ++i) {
                if((match&1)!=0) {
                    nucl_counts[i]+=2.0*frac;
                    nucl_pair_count[i][i]+=2.0*frac;
                }
                match=match>>1;
            }
        }
        void averageAmbig(int mask1, int mask2) {
            int count1=0, count2=0;
            int m1=mask1, m2=mask2;
            for(int i=0; i!=4; ++i) {
                if((m1&1)!=0) ++count1;
                if((m2&1)!=0) ++count2;
                m1=m1>>1;
                m2=m2>>1;
            }
            double frac =1.0/(count1*count2);
            for(int i=0; i!=4; ++i) {
                if((mask1&1)!=0)
                    for(int j=0; j!=4; ++j) {
                        if((mask2&1)!=0) {
                            nucl_counts[i]+=frac;
                            nucl_counts[j]+=frac;
                            nucl_pair_count[i][j]+=frac;
                            nucl_pair_count[j][i]+=frac;
                        }
                        mask2=mask2>>1;
                    }
                mask1=mask1>>1;
            }

        }
    }
    private PairStat getPairStat(int[] s1, int[] s2) {
        PairStat a = new PairStat();

        int length = Math.min(s1.length, s2.length);

        for(int i=0; i<length; ++i) {
            int c1 = s1[i];
            int c2 = s2[i];
            if(c1==Seq.N || c2==Seq.N) continue;
            ++a.total_non_gap;
            if(c1<=Seq.T && c2<=Seq.T) {
                a.updateStat(c1, c2);
            }
            else {
                int m1 = Seq.getMask(c1);
                int m2 = Seq.getMask(c2);
                switch(ambiguityMode) {
                    case "resolve":
                        int match = m1&m2;
                        if(match != 0) {
                            a.resolveMatch(match);
                            break;
                        }
                    case "average":
                        a.averageAmbig(c1,c2);
                        break;
                    case "skip":
                        break;
                }
            }
        }
        return a;
    }

    private double tn93(int[] s1, int[] s2) {
        PairStat ps = getPairStat(s1, s2);

        double[] nucl_freq = new double[4];
        for(int i=0; i<4; ++i) nucl_freq[i] = ps.nucl_counts[i]/2.0/ps.total_non_gap;
        double p1 = ps.nucl_pair_count[Seq.A][Seq.G]/ps.total_non_gap;
        double p2 = ps.nucl_pair_count[Seq.C][Seq.T]/ps.total_non_gap;
        double q = (ps.nucl_pair_count[Seq.A][Seq.C]+ps.nucl_pair_count[Seq.A][Seq.T]+ps.nucl_pair_count[Seq.C][Seq.G]+
                ps.nucl_pair_count[Seq.G][Seq.T])/ps.total_non_gap;
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
        while (sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if (line.length() == 0) continue;
            if (line.charAt(0) == '>') {
                if (name.length() != 0)
                    seqs.add(new Seq(name, seq));
                name = line.substring(1);
                seq = "";
            } else seq = seq.concat(line);
        }
        if (name.length() != 0)
            seqs.add(new Seq(name, seq));
        return seqs;
    }

    private static LinkedList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        LinkedList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
}
