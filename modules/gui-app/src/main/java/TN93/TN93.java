package TN93;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.Observable;

//Parallelization
import java.util.concurrent.*;


import static java.lang.Math.log;

public class TN93 extends Observable {
    // Used for parallelization
    private static class Triplet<A, B, C> {
        final A first;
        final B second;
        final C third;
    
        public Triplet(A first, B second, C third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }
    }


    private float edgeThreshold = 1;
    private File inputFile;
    private File outputFile;
    private String ambiguityHandling;
    private int cores = 1;
    private double max_ambiguity_fraction = -1; // -1 means no limit (default)
    private Map<Seq, Double> ambig_fractions = new HashMap<>();

    public void setEdgeThreshold(float edgeThreshold) {
        this.edgeThreshold = edgeThreshold;
    }

    public void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(File outputFile) {
        this.outputFile = outputFile;
    }

    public void setAmbiguityHandling(String ambiguityHandling) {
        this.ambiguityHandling = ambiguityHandling;
    }

    public void setMaxAmbiguityFraction(double max_ambiguity_fraction) {
        this.max_ambiguity_fraction = max_ambiguity_fraction;
    }

    public void setCores(int cores) {
        this.cores = cores;
    }

    final static int[][] resolutions = {
        // A,C,G,T
        {1, 0, 0, 0},  // A             -> A (0) (Adenine)
        {0, 1, 0, 0},  // C             -> C (1) (Cytosine)
        {0, 0, 1, 0},  // G             -> G (2) (Guanine)
        {0, 0, 0, 1},  // T             -> T (3) (Thymine)
        {0, 0, 0, 1},  // T             -> U (4) (Uracil)
        {1, 0, 1, 0},  // A | G         -> R (5) (Either Purine)
        {0, 1, 0, 1},  // C | T         -> Y (6) (Either Pyrimidine)
        {0, 1, 1, 0},  // C | G         -> S (7)
        {1, 0, 0, 1},  // A | T         -> W (8)
        {0, 0, 1, 1},  // G | T         -> K (9)
        {1, 1, 0, 0},  // A | C         -> M (10)
        {0, 1, 1, 1},  // C | G | T     -> B (11) (Not Adenine)
        {1, 0, 1, 1},  // A | G | T     -> D (12) (Not Cytosine)
        {1, 1, 0, 1},  // A | C | T     -> H (13) (Not Guanine)
        {1, 1, 1, 0},  // A | C | G     -> V (14) (Not Thymine)
        {1, 1, 1, 1},  // A | C | G | T -> N (15)
        {1, 1, 1, 1},  // A | C | G | T -> ? (16)
        {0, 0, 0, 0},  // GAP           -> - (17)
    };

    final static double[] resolutionsCount = {
        1.0,  // A
        1.0,  // C
        1.0,  // G
        1.0,  // T
        1.0,  // U
        1.0 / 2.0,  // R
        1.0 / 2.0,  // Y
        1.0 / 2.0,  // S
        1.0 / 2.0,  // W
        1.0 / 2.0,  // K
        1.0 / 2.0,  // M
        1.0 / 3.0,  // B
        1.0 / 3.0,  // D
        1.0 / 3.0,  // H
        1.0 / 3.0,  // V
        1.0 / 4.0,  // N
        1.0 / 4.0,  // ?
        0.0,  // GAP
    };

    public void tn93Fasta() {
        PrintWriter f = null;
        try {
            System.out.println("Reading input file...");
            ArrayList<Seq> seqs = read_fasta(inputFile);
            System.out.println("Calculating distances...");
            tn93(seqs);
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }

    public void tn93(ArrayList<Seq> seqs){
        if ("resolve".equals(ambiguityHandling)) {
            ambig_fractions = count_ambiguities(seqs);
        }

        if (cores >= seqs.size()) {
            cores = seqs.size()-1;
        }

        System.out.println("Creating thread pool with " + cores + " threads...");
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        List<Future<Triplet<Integer, Integer, Double>>> futures = new ArrayList<>();
        AtomicReference<PrintWriter> writerRef = new AtomicReference<>();
        try {
            writerRef.set(new PrintWriter(outputFile));
        } catch (IOException e) {
            e.printStackTrace();
        }

        writerRef.get().println("Source,Target,Distance");


        long pairs_count = (seqs.size() * seqs.size() - seqs.size())/2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        System.out.println("Submitting jobs...");
        for (int i = 1; i < seqs.size(); ++i) {
            System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                final int row = i;
                final int col = j;
                final Seq seq1 = seqs.get(row);
                final Seq seq2 = seqs.get(col);
                
                futures.add(executor.submit( () -> {
                    double d = tn93(seq1, seq2);
                    if (d < this.edgeThreshold)
                        writerRef.get().println(String.format("%s,%s,%f", seq1.getName(), seq2.getName(), d));
                    return new Triplet<>(row, col, d);
                }));

                if (futures.size() > 1000 || pairs_count < 1000 ) {
                    for (Future<Triplet<Integer, Integer, Double>> future : futures) {
                        try {
                            future.get();         
                        }
                        catch (InterruptedException | ExecutionException e) {
                            e.printStackTrace();
                        }
        
                        current_pair++;
                        // There is some issue here when input has exactly 2 sequences
                        if (pairs_count < 100 || current_pair % (pairs_count / 100) == 0 ) {
                            estimatedTime = System.nanoTime() - startTime;
                            int percCompleted = (int) (current_pair*100/pairs_count);
                            System.out.print(String.format("%d%% completed in ", percCompleted));
                            System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                            System.out.println(" sec                                ");
                            setChanged();
                            notifyObservers(percCompleted);
                        } 
                    }
                    futures.clear();
                }
            }
        }

        for (Future<Triplet<Integer, Integer, Double>> future : futures) {
            try {
                future.get();         
            }
            catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }

            current_pair++;
            if (pairs_count < 100 || current_pair % (pairs_count / 100) == 0) {
                estimatedTime = System.nanoTime() - startTime;
                int percCompleted = (int) (current_pair*100/pairs_count);
                System.out.print(String.format("%d%% completed in ", percCompleted));
                System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                System.out.println(" sec                                ");
                setChanged();
                notifyObservers(percCompleted);
            } 
        }

        executor.shutdown();
        writerRef.get().flush();
        writerRef.get().close();
        setChanged();
        notifyObservers(100);
        return;
    }

    private HashMap<Seq,Double> count_ambiguities(ArrayList<Seq> seqs) {
        HashMap<Seq,Double> ambig_fractions = new HashMap<>();
        int seq_length = seqs.get(0).getSeq_enc().length;
        for (Seq seq : seqs) {
            int ambigs_count = 0;
            for (int i = 0; i < seq.getSeq_enc().length; ++i) {
                int enc = seq.getSeq_enc()[i];
                if (enc > 4 && enc != 17) {
                    ambigs_count++;
                }
            }
            ambig_fractions.put(seq, ambigs_count / (double) seq_length);
        }
        return ambig_fractions;
    }

    private double tn93(Seq s1, Seq s2) {
        double[][] nucl_pair_counts = countPairwiseNucl(s1, s2);
        double[] nucl_counts = getNuclCounts(nucl_pair_counts);

        double total_non_gap = 2.0 / Arrays.stream(nucl_counts).sum();

        double[] nucl_freq = new double[4];
        double auxd = 1.0 / Arrays.stream(nucl_counts).sum();
        for(int i=0; i<4; ++i) 
            nucl_freq[i] = nucl_counts[i] * auxd;

        double dist = 0.0;
        double AG_counts = nucl_pair_counts[Seq.A][Seq.G] + nucl_pair_counts[Seq.G][Seq.A];
        double AG = AG_counts * total_non_gap;
        double CT_counts = nucl_pair_counts[Seq.C][Seq.T] + nucl_pair_counts[Seq.T][Seq.C];
        double CT = CT_counts * total_non_gap;
        double matching = (nucl_pair_counts[Seq.A][Seq.A] + nucl_pair_counts[Seq.C][Seq.C] + nucl_pair_counts[Seq.T][Seq.T] + nucl_pair_counts[Seq.G][Seq.G]) * total_non_gap;
        double tv = 1 - (AG + CT + matching);

        boolean useK2P = false;
        for(int i=0; i<4; ++i) if (nucl_freq[i] == 0.0) useK2P = true;

        if (useK2P) {
            AG = 1 - 2 * (AG + CT) - tv;
            CT = 1 - 2 * tv;
            if (AG > 0 && CT > 0) {
                dist = -0.5 * log(AG) - 0.25 * log(CT);
            } else {
                dist = 1.0;
            }
        } else {
            double fR = nucl_freq[Seq.A] + nucl_freq[Seq.G];
            double fY = nucl_freq[Seq.C] + nucl_freq[Seq.T];
            double K1 = 2 * nucl_freq[Seq.A] * nucl_freq[Seq.G] / fR;
            double K2 = 2 * nucl_freq[Seq.C] * nucl_freq[Seq.T] / fY;
            double K3 = 2 * ( fR * fY - nucl_freq[Seq.A] * nucl_freq[Seq.G] * fY / fR - nucl_freq[Seq.C] * nucl_freq[Seq.T] * fR / fY);
            AG = 1 - AG / K1 - 0.5 * tv / fR;
            CT = 1 - CT / K2 - 0.5 * tv / fY;
            tv = 1 - 0.5 * tv / fY / fR;
            dist = -K1 * log(AG) - K2 * log(CT) - K3 * log(tv);
        }
        return dist;
    }


    private double[][] countPairwiseNucl(Seq s1, Seq s2) {
        double nucl_pair_counts[][] = new double[4][4];
        if ("resolve".equals(ambiguityHandling)) {
            if (max_ambiguity_fraction != -1)
                if (ambig_fractions.get(s1) > max_ambiguity_fraction || ambig_fractions.get(s2) > max_ambiguity_fraction)
                    return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
            return countNucl_resolve(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        }
        else if ("average".equals(ambiguityHandling))
            return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else if ("gapmm".equals(ambiguityHandling))
            return countNucl_gapmm(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else if ("skip".equals(ambiguityHandling))
            return countNucl_skip(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else
        {
            if (max_ambiguity_fraction != -1)
                if (ambig_fractions.get(s1) > max_ambiguity_fraction || ambig_fractions.get(s2) > max_ambiguity_fraction)
                    return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
            return countNucl_resolve(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        }
    }


    private static double[] getNuclCounts(double[][] nucl_pair_counts) {
        double[] nucl_freq = new double[4];
        for(int i=0; i<4; ++i) {
            for(int j=0; j<4; ++j) {
                nucl_freq[i] += nucl_pair_counts[i][j];
                nucl_freq[j] += nucl_pair_counts[i][j];
            }
        }
        return nucl_freq;
    }


    private double[][] countNucl_resolve(int[] s1, int[] s2, double[][]nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for(int i=0; i<L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];
            if (c1 == 17 && c2 == 17) continue;                 // both are gaps; continue
            if (c1 < 4 && c2 < 4) {                             // Neither is ambiguous
                nucl_pair_counts[c1][c2]++;
            }
            else {
                if (c1 < 4) {                                   // if c1 is not ambiguous, c2 is
                    if (resolutionsCount[c2] > 0) {             // if c2 is not a gap
                        if (resolutions[c2][c1] == 1) {         // if c2 can resolve to c1               
                            nucl_pair_counts[c1][c1]++;          // Resolve c2 to c1 
                            continue;
                        }
                        for (int j=0; j<4; j++) {               // else: average over all possible resolutions
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4){                               // c2 is not ambiguous, c1 is
                    if (resolutionsCount[c1] > 0) {
                        if (resolutions[c1][c2] == 1) {
                            nucl_pair_counts[c2][c2]++;
                            continue;
                        }
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                } 
                else {                                          // Both c1 and c2 are ambiguous
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) { // if both are not gaps
                        int matched = 0;
                        boolean[] positive_match = new boolean[4];
                        for (int j=0; j<4; j++) { 
                            if (resolutions[c1][j] == 1 && resolutions[c2][j] == 1) {
                                positive_match[j] = true;
                                matched++;
                            }
                        }
                        if (matched > 0) {
                            double norm2 = 1.0/matched;
                            for (int k=0; k<4; k++) {
                                if (positive_match[k]) {
                                    nucl_pair_counts[k][k] += norm2;
                                }
                            }
                            continue;
                        }
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k] == 1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                } 
            }
        }
        return nucl_pair_counts;
    }


    private double[][] countNucl_average(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for(int i=0; i<L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 == 17 || c2 == 17) continue;   

            if (c1 < 4 && c2 < 4) {
                nucl_pair_counts[c1][c2]++;
            }
            else {
                if (c1 < 4) {
                    if (resolutionsCount[c2] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4){
                    if (resolutionsCount[c1] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                } 
                else {
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j]==1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k]==1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                } 
            }
        }
        return nucl_pair_counts;
    }


    private static double[][] countNucl_gapmm(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for (int i = 0; i < L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 == 17 && c2 == 17) continue;

            if (c1 < 4 && c2 < 4) 
                nucl_pair_counts[c1][c2]++;
            else {
                if (c1 == 17 || c2 == 17) {
                    if (c1 == 17) 
                        c1 = 15;
                    else 
                        c2 = 15;
                }
                if (c1 < 4) {
                    if (resolutionsCount[c2] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4) {
                    if (resolutionsCount[c1] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                }
                else {
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j]==1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k]==1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return nucl_pair_counts;
    }


    private static double[][] countNucl_skip(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);

        for (int i = 0; i < L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 < 4 && c2 < 4) {
                nucl_pair_counts[c1][c2]++;
            }
        }
        return nucl_pair_counts;
    }


    public static ArrayList<Seq> read_seqs(Scanner sc) {
        ArrayList<Seq> seqs = new ArrayList<Seq>();
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


    private static ArrayList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        ArrayList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
}
