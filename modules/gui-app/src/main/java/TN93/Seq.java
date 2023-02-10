package TN93;

import java.util.HashMap;
import java.util.Map;

public class Seq {
    final static int A=0, C=1, G=2, T=3, Gap=17;
    static Map<Character, Integer> nucl;
    static {
        nucl = new HashMap<Character, Integer>();
        nucl.put('A', 0);
        nucl.put('C', 1);
        nucl.put('G', 2);
        nucl.put('T', 3);
        nucl.put('U', 4);
        nucl.put('R', 5);
        nucl.put('Y', 6);
        nucl.put('S', 7);
        nucl.put('W', 8);
        nucl.put('K', 9);
        nucl.put('M', 10); 
        nucl.put('B', 11); 
        nucl.put('D', 12); 
        nucl.put('H', 13); 
        nucl.put('V', 14); 
        nucl.put('N', 15); 
        nucl.put('?', 16);
        nucl.put('-', 17); 
    }

    private String name;
    private String seq;
    private int[] seq_enc;
    Seq(String name, String seq) {
        this.name = name;
        this.seq = seq.toUpperCase();
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