package TN93;

public class Seq {
    private final static String NUCLEOTIDES = "ACGTRYKMSWBDHVN";
    private static int[] CIPHER;
    private final static int OFFSET = Character.getNumericValue('A');
    private final static int GAP = '-';
    private final static int GAP_CODE = NUCLEOTIDES.length()+1;
    static {

        CIPHER = new int[Character.getNumericValue('Y')- OFFSET];
        for(int i = 0; i!= NUCLEOTIDES.length(); ++i) {
            CIPHER[Character.getNumericValue(NUCLEOTIDES.charAt(i))- OFFSET]=i;
        }
    }
    private final static int[] masks = {
            0x0001, //A
            0x0010, //C
            0x0100, //G
            0x1000, //T
            0x0101, //R
            0x1010, //Y
            0x1100, //K
            0x0011, //M
            0x0110, //S
            0x1001, //W
            0x1110, //B
            0x1101, //D
            0x1011, //H
            0x0111, //V
            0x1111, //N
            0x0000, //-

    };

    private String name;
    private String seq;
    private int[] seq_enc;

    private static int get_code(char c) {
        return c==GAP? GAP_CODE: CIPHER[Character.getNumericValue(c) - OFFSET];
    }

    Seq(String name, String seq) {
        this.name = name;
        this.seq = seq.toUpperCase();
        seq_enc = new int[seq.length()];
        for(int i=0;i<seq.length();++i) {
            char c = seq.charAt(i);
            seq_enc[i] = get_code(c);
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