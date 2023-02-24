import TN93.TN93;
import picocli.CommandLine;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Observable;
import java.util.Observer;
import javax.swing.filechooser.FileNameExtensionFilter;

import static javax.swing.JOptionPane.showMessageDialog;

@CommandLine.Command(name = "tn93", mixinStandardHelpOptions = true, version = "0.1")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input file with sequences",
            paramLabel = "FILE")
    private File inputFile;
    @CommandLine.Option(names={"-o", "--outFile"},
            description="output file with distances",
            paramLabel = "FILE")
    private File outputFile;
    @CommandLine.Option(names={"-s", "--server"}, description="run jetty server")
    private boolean is_server;
    @CommandLine.Option(names={"-t", "--edge-threshold"},
            description="edges above the threshold are not reported in output")
    private String edgeThresholdString;
    @CommandLine.Option(names={"-a", "--ambiguity", "--ambiguities"},
            description="How to handle ambiguous nucleotides. One of [resolve, average, gapmm, skip]")
    private String ambiguityHandling;
    @CommandLine.Option(names={"-g", "--fraction"},
            description="Maximum allowable fraction of ambiguities allowed for 'resolve' mode. If exceeded, use 'average' mode.")
    private float max_ambiguity_fraction;

    public void run() {
        if(is_server) {
            try {
                run_server();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        if(inputFile == null || outputFile == null) {
            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    createAndShowGUI();
                }
            });
        }
        else {
            TN93 tn93 = new TN93();
            tn93.setEdgeThreshold(Float.parseFloat(edgeThresholdString));
            tn93.setInputFile(inputFile);
            tn93.setOutputFile(outputFile);
            System.out.println(ambiguityHandling);
            tn93.setAmbiguityHandling(ambiguityHandling);
            tn93.setMaxAmbiguityFraction(max_ambiguity_fraction);
            tn93.tn93Fasta();
        }
    }
    private static void run_server() throws InterruptedException {
        System.out.println("To stop the server press Ctrl-C");
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                System.out.println("Stopping the server...");
            }
        });
        //TODO: Run Server
        while(true) Thread.sleep(1000);
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("TN93");
        JPanel mainPane = new JPanel();
        mainPane.setLayout(new BoxLayout(mainPane, BoxLayout.Y_AXIS));
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        TN93_Panel panel = new TN93_Panel();
        mainPane.add(panel);
        panel.setOpaque(true); //content panes must be opaque
        frame.setContentPane(mainPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
}

class TN93_Panel extends JPanel implements ActionListener, Observer {
    private float edgeThreshold;
    private JButton inBut, outBut, runBut;
    private JTextField fastaTextField, edgeListTextField, edgeThresholdField;
    private JProgressBar progress;
    private ButtonGroup ambiguityHandlingGroup;
    private JRadioButton resolveBut, averageBut, gapmmBut, skipBut;
    private String ambiguityHandling;

    private File fastaFile, edgeListFile;
    private TN93 tn93;


    TN93_Panel() {
        tn93 = new TN93();
        tn93.addObserver(this);

        inBut = new JButton("Load Fasta");
        outBut = new JButton("Save as: Edge List CSV File");
        runBut = new JButton("Run TN93");

        ambiguityHandlingGroup = new ButtonGroup();
        resolveBut = new JRadioButton("Resolve");
        averageBut = new JRadioButton("Average");
        gapmmBut = new JRadioButton("Gapmm");
        skipBut = new JRadioButton("Skip");

        ambiguityHandlingGroup.add(resolveBut);
        ambiguityHandlingGroup.add(averageBut);
        ambiguityHandlingGroup.add(gapmmBut);
        ambiguityHandlingGroup.add(skipBut);

        resolveBut.setSelected(true);

        fastaTextField = new JTextField(20);
        edgeListTextField = new JTextField(20);
        edgeThresholdField = new JTextField("0.015");
        progress = new JProgressBar(0, 100);
        progress.setStringPainted(true);

        inBut.setActionCommand("loadFasta");
        outBut.setActionCommand("specifyEdgeListFile");
        runBut.setActionCommand("runTN93");

        inBut.addActionListener(this);
        outBut.addActionListener(this);
        runBut.addActionListener(this);

        setLayout(new BorderLayout());

        JPanel fastaPanel = new JPanel();
        fastaPanel.setLayout(new GridLayout(4, 2));
        fastaPanel.add(new JLabel("Maximum edge length:", JLabel.RIGHT));
        fastaPanel.add(edgeThresholdField);
        fastaPanel.add(fastaTextField);
        fastaPanel.add(inBut);
        fastaPanel.add(edgeListTextField);
        fastaPanel.add(outBut);
        fastaPanel.add(progress);
        fastaPanel.add(runBut);
        add(fastaPanel, BorderLayout.NORTH);

        JPanel ambigsPanel = new JPanel();        
        ambigsPanel.setLayout(new GridLayout(1, 4));
        ambigsPanel.add(resolveBut);
        ambigsPanel.add(averageBut);
        ambigsPanel.add(gapmmBut);
        ambigsPanel.add(skipBut);
        ambigsPanel.setBorder(BorderFactory.createTitledBorder("Ambiguity Handling"));
        add(ambigsPanel, BorderLayout.SOUTH);
    }

    public void actionPerformed(ActionEvent e) {
        if("loadFasta".equals(e.getActionCommand())) {
            JFileChooser fileopen =  new JFileChooser();
            FileNameExtensionFilter filter = new FileNameExtensionFilter("FASTA FILES", "fa", "fas", "fasta");
            fileopen.addChoosableFileFilter(filter);
            if(JFileChooser.APPROVE_OPTION == fileopen.showDialog(null, "Open Fasta file")) {
                fastaFile = fileopen.getSelectedFile();
                fastaTextField.setText(fastaFile.getName());
            }
        }
        else if("specifyEdgeListFile".equals(e.getActionCommand())) {
            JFileChooser fileopen =  new JFileChooser();
            //fileopen.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            if(JFileChooser.APPROVE_OPTION == fileopen.showDialog(null, "Save as: Edge List CSV File")) {
                edgeListFile = fileopen.getSelectedFile();
                edgeListTextField.setText(edgeListFile.getName());
            }
        }
        else if("runTN93".equals(e.getActionCommand())) {
            if(fastaFile == null) {
                showMessageDialog(null, "Specify an input Fasta file!");
                return;
            }
            if(edgeListFile == null) {
                showMessageDialog(null, "Specify an output file using Save as!");
                return;
            }
            if(!parseEdgeThreshold()) {
                showMessageDialog(null, "Threshold should be number!");
            }
            inactivatePanel();
            tn93.setEdgeThreshold(edgeThreshold);
            tn93.setInputFile(fastaFile);
            tn93.setOutputFile(edgeListFile);
            
            if(resolveBut.isSelected()) 
                tn93.setAmbiguityHandling("resolve");
            else if(averageBut.isSelected()) 
                tn93.setAmbiguityHandling("average");
            else if(gapmmBut.isSelected()) 
                tn93.setAmbiguityHandling("gapmm");
            else if(skipBut.isSelected()) 
                tn93.setAmbiguityHandling("skip");
            else 
                tn93.setAmbiguityHandling("resolve");
            
            SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
                @Override
                protected Void doInBackground() {
                    tn93.tn93Fasta();
                    return null;
                }
                @Override
                protected void done() {
                    activatePanel();
                }
            };
            worker.execute();
        }
    }
    private void inactivatePanel() {
        runBut.setEnabled(false);
        inBut.setEnabled(false);
        outBut.setEnabled(false);
        progress.setValue(0);
    }
    private void activatePanel() {
        runBut.setEnabled(true);
        inBut.setEnabled(true);
        outBut.setEnabled(true);
    }

    private boolean parseEdgeThreshold() {
        try {
            edgeThreshold = Float.parseFloat(edgeThresholdField.getText());
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
    public void update(Observable obj, Object arg) {
        this.progress.setValue((Integer) arg);
    }

}

class StatusPanel extends JPanel {
    protected JTextArea ta;

    StatusPanel() {
        ta = new JTextArea(10, 40);
        add(new JScrollPane(ta), BorderLayout.PAGE_START);
    }
}
