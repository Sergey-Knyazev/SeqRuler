import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) {
        String a = "ACGTN";
        String b = "TGCAA";

        String[] DNAs = new String[2];
        DNAs[0] = a;
        DNAs[1] = b;

        int A = 0; int C = 1; int G = 2; int T = 3;

        int[][] nucl_counts = new int[DNAs.length][4];
        int[] unumbig_nucl_counts = new int[DNAs.length];
        List<List<List<Integer>>> nucl_groups = new ArrayList<>(DNAs.length);

        for(int i=0; i!=DNAs.length; ++i) {
            nucl_groups.add(new ArrayList<List<Integer>>(4));
            for(int j=0; j!=4; ++j) nucl_groups.get(i).add(new ArrayList<Integer>());
            for(int j=0; j!=DNAs[i].length(); ++j) {
                char c = DNAs[i].charAt(j);
                if (c == 'A')      {++nucl_counts[i][A]; nucl_groups.get(i).get(A).add(j);}
                else if (c == 'C') {++nucl_counts[i][C]; nucl_groups.get(i).get(C).add(j);}
                else if (c == 'G') {++nucl_counts[i][G]; nucl_groups.get(i).get(G).add(j);}
                else if (c == 'T') {++nucl_counts[i][T]; nucl_groups.get(i).get(T).add(j);}
                else continue;
                ++unumbig_nucl_counts[i];
            }
        }
    }
}
