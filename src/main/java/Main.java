import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static java.lang.Math.log;

public class Main {
    public static void main(String[] args) {
        String s1 = "ACGTN";
        String s2 = "TGCAA";

        int A=0,C=1,G=2,T=3;

        Map<Character, Integer> nucl = new HashMap<>();
        nucl.put('A', A);
        nucl.put('C', C);
        nucl.put('G', G);
        nucl.put('T', T);

        List<Integer> s1_enc = new LinkedList<>();
        List<Integer> s2_enc = new LinkedList<>();
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

        double d = -k_ag*log(1-p1/k_ag-q/(2*g_r))
               -k_tc*log(1-p2/k_tc-q/(2*g_y))
               -(k_ry-k_ag*g_y-k_tc*g_r)*log(1-q/k_ry);
        System.out.println(q/k_ry);
    }
}
