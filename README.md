# SeqRuler

GUI interface for calculating TN93 distance between all sequences in the input fasta file.

## Download or Compilation

- Download:
[tn93.jar](https://github.com/CDCgov/SeqRuler/releases/download/v1.0/tn93.jar)


- Compilation:
```bash
mvn clean install
```

## Running

```bash
java -jar tn93.jar
```

## Help

```bash
java -jar tn93.jar -h

Usage: tn93 [-hsV] [-i=FILE] [-o=FILE] [-t=<edgeThresholdString>]
  -h, --help           Show this help message and exit.
  -i, --inFile=FILE    input file with sequences
  -o, --outFile=FILE   output file with distances
  -s, --server         run jetty server
  -t, --edge-threshold=<edgeThresholdString>
                       edges above the threshold are not reported in output
  -V, --version        Print version information and exit.
```
