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
    @CommandLine.Option(names={"-c", "--cores"},
            description="Number of cores to use for parallel processing.")
    private int cores;

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
            tn93.setCores(cores);
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
        TN93_Panel panel = new TN93_Panel(frame);
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
    private JTextField maxAmbiguityFractionField;
    private JLabel maxAmbiguityFractionLabel;
    private JCheckBox useMaxCoresCheckbox;
    private JTextField numCoresField;

    private File fastaFile, edgeListFile;
    private TN93 tn93;
    private JFrame frame;


    TN93_Panel(JFrame frame) {
        this.frame = frame;
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
        ambigsPanel.setLayout(new GridLayout(2, 1));
        ambigsPanel.setBorder(BorderFactory.createTitledBorder("Ambiguity Handling"));

        JPanel radioButtonsPanel = new JPanel();
        radioButtonsPanel.setLayout(new GridLayout(1, 4));
        radioButtonsPanel.add(resolveBut);
        radioButtonsPanel.add(averageBut);
        radioButtonsPanel.add(gapmmBut);
        radioButtonsPanel.add(skipBut);
        ambigsPanel.add(radioButtonsPanel, BorderLayout.NORTH);

        JPanel maxAmbigsPanel = new JPanel();
        maxAmbigsPanel.setLayout(new GridLayout(1, 2));
        maxAmbiguityFractionLabel = new JLabel("Maximum ambiguity fraction: ", JLabel.CENTER);
        maxAmbigsPanel.add(maxAmbiguityFractionLabel);
        maxAmbiguityFractionField = new JTextField("0.05");
        maxAmbigsPanel.add(maxAmbiguityFractionField, BorderLayout.SOUTH);
        
        ambigsPanel.add(maxAmbigsPanel);
        add(ambigsPanel, BorderLayout.CENTER);

        resolveBut.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ambigsPanel.add(maxAmbigsPanel);
                ambigsPanel.setLayout(new GridLayout(2, 1));
                ambigsPanel.revalidate();
                ambigsPanel.repaint();
                frame.pack();
            }
        });

        ActionListener hideMaxAmbiguityListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ambigsPanel.remove(maxAmbigsPanel);
                ambigsPanel.setLayout(new GridLayout(1, 1));
                ambigsPanel.revalidate();
                ambigsPanel.repaint();
                frame.pack();
            }
        };

        averageBut.addActionListener(hideMaxAmbiguityListener);
        gapmmBut.addActionListener(hideMaxAmbiguityListener);
        skipBut.addActionListener(hideMaxAmbiguityListener);

        JPanel coresPanel = new JPanel();
        coresPanel.setLayout(new GridLayout(1, 2));
        useMaxCoresCheckbox = new JCheckBox("Use all available cores", true);
        numCoresField = new JTextField(Integer.toString(Runtime.getRuntime().availableProcessors()));
        numCoresField.setEnabled(false);

        useMaxCoresCheckbox.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (useMaxCoresCheckbox.isSelected()) {
                    numCoresField.setEnabled(false);
                    numCoresField.setText(Integer.toString(Runtime.getRuntime().availableProcessors()));
                } else {
                    numCoresField.setEnabled(true);
                }
            }
        });

        coresPanel.add(new JLabel("Number of cores to use:", JLabel.CENTER));
        coresPanel.add(numCoresField);
        coresPanel.add(useMaxCoresCheckbox);

        coresPanel.setBorder(BorderFactory.createTitledBorder("Processing"));
        add(coresPanel, BorderLayout.SOUTH);
    }

    public void actionPerformed(ActionEvent e) {
        if("loadFasta".equals(e.getActionCommand())) {
            FileDialog fileDialog = new FileDialog(new Frame(), "Open Fasta file", FileDialog.LOAD);
            fileDialog.setFilenameFilter((dir, name) -> name.endsWith(".fa") || name.endsWith(".fas") || name.endsWith(".fasta"));
            fileDialog.setVisible(true);
        
            String selectedFileDirectory = fileDialog.getDirectory();
            String selectedFileName = fileDialog.getFile();
        
            if (selectedFileDirectory != null && selectedFileName != null) {
                fastaFile = new File(selectedFileDirectory, selectedFileName);
                fastaTextField.setText(fastaFile.getName());
            }
        }
        
        else if ("specifyEdgeListFile".equals(e.getActionCommand())) {
            FileDialog fileDialog = new FileDialog(new Frame(), "Save as: Edge List CSV File", FileDialog.SAVE);
            fileDialog.setFilenameFilter((dir, name) -> name.endsWith(".csv"));
            fileDialog.setVisible(true);
        
            String selectedFileDirectory = fileDialog.getDirectory();
            String selectedFileName = fileDialog.getFile();
        
            if (selectedFileDirectory != null && selectedFileName != null) {
                edgeListFile = new File(selectedFileDirectory, selectedFileName);
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
            
            if(resolveBut.isSelected()) {
                tn93.setAmbiguityHandling("resolve");
                tn93.setMaxAmbiguityFraction(Double.parseDouble(maxAmbiguityFractionField.getText()));
            } else if(averageBut.isSelected()) 
                tn93.setAmbiguityHandling("average");
            else if(gapmmBut.isSelected()) 
                tn93.setAmbiguityHandling("gapmm");
            else if(skipBut.isSelected()) 
                tn93.setAmbiguityHandling("skip");
            else 
                tn93.setAmbiguityHandling("resolve");

            if(useMaxCoresCheckbox.isSelected()) {
                tn93.setCores(Runtime.getRuntime().availableProcessors());
            } else {
                try {
                    tn93.setCores(Integer.parseInt(numCoresField.getText()));
                } catch (NumberFormatException ex) {
                    showMessageDialog(null, "Number of cores should be number!");
                }
            }
        
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
