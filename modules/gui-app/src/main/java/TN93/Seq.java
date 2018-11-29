package TN93;

import java.util.HashMap;
import java.util.Map;

public class Seq {
    final static int A=0,C=1,G=2,T=3,N=4;
    static Map<Character, Integer> nucl;
    static {
        nucl = new HashMap<Character, Integer>();
        nucl.put('A', A);
        nucl.put('C', C);
        nucl.put('G', G);
        nucl.put('T', T);
        nucl.put('N', N);
    }

    private String name;
    private String seq;
    private int[] seq_enc;
    Seq(String name, String seq) {
        this.name = name;
        this.seq = seq;
        seq_enc = new int[seq.length()];
        for(int i=0;i<seq.length();++i) {
            char c = seq.charAt(i);
            if(!nucl.containsKey(c)) seq_enc[i] = nucl.get('N');
            else seq_enc[i] = nucl.get(c);
        }
    }
    public String getName() {
        return name;
    }

    public String getSeq() {
        return seq;
    }

    public int[] getSeq_enc() {
        return seq_enc;
    }
}