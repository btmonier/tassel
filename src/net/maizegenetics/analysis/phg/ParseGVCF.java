package net.maizegenetics.analysis.phg;

/**
 * @author Terry Casstevens Created June 28, 2017
 */

import com.google.common.collect.ImmutableList;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;


public class ParseGVCF {

    private static final Logger myLogger = Logger.getLogger(ParseGVCF.class);

    private static final int NUM_LINES_PER_BLOCK = 10;

    private ParseGVCF() {
        // utility
    }

    public static BlockingQueue<Future<ProcessLines>> parse(String filename) {
        myLogger.info("ParseGVCF: filename: " + filename);
        ExecutorService pool = Executors.newWorkStealingPool();
        BlockingQueue<Future<ProcessLines>> queue = new LinkedBlockingQueue<>();
        pool.submit(new ReadLines(filename, queue, pool));
        return queue;
    }

    public static class ProcessLines implements Callable<ProcessLines> {

        private final List<String> myLines;
        private final int myStartLine;
        private final int myNumLines;
        private final boolean myIsHeader;
        private List<GVCFLine> myProcessedLines = null;

        // first one containing header lines (starts with #)
        public ProcessLines(List<String> headerLines) {
            myLines = Collections.unmodifiableList(headerLines);
            myStartLine = 1;
            myNumLines = headerLines.size();
            myIsHeader = true;
        }

        // data lines following header lines
        public ProcessLines(int lineNum, List<String> lines) {
            myLines = Collections.unmodifiableList(lines);
            myStartLine = lineNum;
            myNumLines = lines.size();
            myIsHeader = false;
        }

        // indicates last element in queue
        public ProcessLines() {
            myLines = null;
            myStartLine = -1;
            myNumLines = 0;
            myIsHeader = false;
        }

        @Override
        public ProcessLines call() throws Exception {

            if (isFinal() || isHeader()) {
                return this;
            }

            List<GVCFLine> temp = new ArrayList<>();
            int currentLineNum = myStartLine;
            for (String current : myLines) {
                temp.add(new GVCFLine(currentLineNum, current));
                currentLineNum++;
            }
            myProcessedLines = temp;

            return this;

        }

        public List<String> lines() {
            return myLines;
        }

        public int startLine() {
            return myStartLine;
        }

        public int numLines() {
            return myNumLines;
        }

        public boolean isHeader() {
            return myIsHeader;
        }

        public boolean isFinal() {
            return myStartLine == -1;
        }

        public List<GVCFLine> processedLines() {
            return myProcessedLines;
        }

    }

    public static class GVCFLine {

        private final String myLine;
        private final int myLineNum;
        // Chromosome
        private String myChromosome = null;
        // Starting physical position (inclusive)
        private int myStartPosition = -1;
        // Ending physical position (inclusive)
        private int myEndPosition = -1;
        // Sequence length
        private int mySeqLength = 0;
        // DP value - Total depth
        private int myDepth = 0;
        // AD values - Allele depths
        private List<Integer> myAlleleDepths = null;
        // Is Homozygous
        private boolean myIsHomozygous = true;
        // Reference (REF)
        private String myReference = null;
        // Alternates (ALT)
        private List<String> myAlternates = null;
        // GT values - genotype
        private List<Integer> myGenotypeIndices = null;
        // GT indices translated to allele string
        private List<String> myGenotypes = null;
        // Whether phased (/ : genotype unphased | : genotype phased)
        private boolean myIsPhased = false;
        // Ploidy
        private int myPloidy = -1;

        public GVCFLine(int lineNum, String line) {
            myLineNum = lineNum;
            myLine = line;
            parseLine();
        }

        /**
         * Line number in the .g.vcf file
         *
         * @return line number
         */
        public int lineNum() {
            return myLineNum;
        }

        @Override
        public String toString() {
            return myLine;
        }

        public String chromosome() {
            return myChromosome;
        }

        /**
         * Start physical position (inclusive)
         *
         * @return start physical position
         */
        public int startPosition() {
            return myStartPosition;
        }

        /**
         * End physical position (inclusive)
         *
         * @return end physical position
         */
        public int endPosition() {
            return myEndPosition;
        }

        /**
         * Length of sequence represented by this line
         *
         * @return sequence length
         */
        public int seqLength() {
            return mySeqLength;
        }

        /**
         * Total depth
         *
         * @return depth
         */
        public int depth() {
            return myDepth;
        }

        public List<Integer> alleleDepths() {
            return myAlleleDepths;
        }

        public boolean isHomozygous() {
            return myIsHomozygous;
        }

        /**
         * @return reference
         *
         * <a href="http://samtools.github.io/hts-specs/VCFv4.2.pdf">http://samtools.github.io/hts-specs/VCFv4.2.pdf</a>
         *
         * REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
         * The value in the POS field refers to the position of the first base in the String. For simple insertions and
         * deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
         * Strings must include the base before the event (which must be reflected in the POS field), unless the event
         * occurs at position 1 on the contig in which case it must include the base after the event; this padding base
         * is not required (although it is permitted) for e.g. complex substitutions or other events where all alleles
         * have at least one base represented in their Strings. If any of the ALT alleles is a symbolic allele (an
         * angle-bracketed ID String) then the padding base is required and POS denotes the coordinate of the base
         * preceding the polymorphism. Tools processing VCF files are not required to preserve case in the allele
         * Strings. (String, Required).
         */
        public String reference() {
            return myReference;
        }

        /**
         * @return alternate base(s)
         *
         * <a href="http://samtools.github.io/hts-specs/VCFv4.2.pdf">http://samtools.github.io/hts-specs/VCFv4.2.pdf</a>
         *
         * ALT - alternate base(s): Comma separated list of alternate non-reference alleles. These alleles do not have
         * to be called in any of the samples. Options are base Strings made up of the bases A,C,G,T,N,*, (case
         * insensitive) or an angle-bracketed ID String or a breakend replacement string as described in the section on
         * breakends. The ‘*’ allele is reserved to indicate that the allele is missing due to a upstream deletion. If
         * there are no alternative alleles, then the missing value should be used. Tools processing VCF files are not
         * required to preserve case in the allele String, except for IDs, which are case sensitive. (String; no
         * whitespace, commas, or angle-brackets are permitted in the ID String itself)
         */
        public List<String> alternates() {
            return myAlternates;
        }

        /**
         * @return genotype
         *
         * <a href="http://samtools.github.io/hts-specs/VCFv4.2.pdf">http://samtools.github.io/hts-specs/VCFv4.2.pdf</a>
         *
         * GT : genotype, encoded as allele values separated by either of / or |. The allele values are 0 for the
         * reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele
         * list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc. For haploid calls, e.g.
         * on Y, male nonpseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call
         * might look like 0/0/1. If a call cannot be made for a sample at a given locus, ‘.’ should be specified for
         * each missing allele in the GT field (for example ‘./.’ for a diploid genotype and ‘.’ for haploid genotype).
         * The meanings of the separators are as follows (see the PS field below for more details on incorporating
         * phasing information into the genotypes): / : genotype unphased | : genotype phased
         */
        public List<Integer> genotypeIndices() {
            return myGenotypeIndices;
        }

        public List<String> genotypes() {
            return myGenotypes;
        }

        public boolean phased() {
            return myIsPhased;
        }

        public int ploidy() {
            return myPloidy;
        }

        private void parseLine() {

            if (myLine == null || myLine.startsWith("#")) {
                myLogger.error(myLineNum + ": " + myLine);
                throw new IllegalStateException("GVCFLine: line shouldn't be null or start with #");
            }

            // 0        1       2       3       4       5       6       7       8       9
            // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  W22
            // 8       17639   .       C       <NON_REF>       .       .       END=17789       GT:DP:GQ:MIN_DP:PL      0:12:99:8:0,140
            // 8       17790   .       T       A,<NON_REF>     175.00  .       DP=5;MLEAC=1,0;MLEAF=1.000,0.000;RAW_MQ=18000.00        GT:AD:DP:GQ:PL:SB       1:0,5,0:5:99:205,0,205:0,0,2,3
            String[] tokens = myLine.split("\t");

            myChromosome = tokens[0];
            myStartPosition = Integer.parseInt(tokens[1]);
            myEndPosition = myStartPosition;

            String id = tokens[2];
            myReference = tokens[3];
            String alt = tokens[4];
            ImmutableList.Builder<String> altBuilder = new ImmutableList.Builder<>();
            altBuilder.add(alt.split(","));
            myAlternates = altBuilder.build();
            String qual = tokens[5];
            String filter = tokens[6];

            String[] info = tokens[7].split(";");
            for (String current : info) {
                String[] keyValue = current.split("=");
                if (keyValue[0].equals("END")) {
                    myEndPosition = Integer.parseInt(keyValue[1]);
                    if (myEndPosition < myStartPosition) {
                        myLogger.error(myLineNum + ": " + myLine);
                        throw new IllegalStateException("GVCFLine: End position less than start position.");
                    }
                }
            }

            mySeqLength = myEndPosition - myStartPosition + 1;

            String[] formats = tokens[8].split(":");
            String[] values = tokens[9].split(":");

            int numFormats = formats.length;
            if (numFormats != values.length) {
                myLogger.error(myLineNum + ": " + myLine);
                throw new IllegalArgumentException("GVCFLine: number of formats not equal to number of values.");
            }

            boolean isHomoBasedOnGenotype = false;
            boolean isAlt = false;
            int alleleDepthSum = 0;
            for (int i = 0; i < numFormats; i++) {
                if (formats[i].equals("DP")) {
                    myDepth = Integer.parseInt(values[i]);
                } else if (formats[i].equals("AD")) {
                    String[] alleleDepths = values[i].split(",");
                    List<Integer> temp = new ArrayList<>();
                    for (String alleleDepth : alleleDepths) {
                        int depthValue = Integer.parseInt(alleleDepth);
                        alleleDepthSum += depthValue;
                        temp.add(depthValue);
                    }
                    myAlleleDepths = Collections.unmodifiableList(temp);
                    if (alleleDepths.length > 1 && !alleleDepths[0].equals("0") && !alleleDepths[1].equals("0")) {
                        myIsHomozygous = false;
                    }
                } else if (formats[i].equals("GT")) {

                    String[] alleles = null;
                    if (values[i].contains("/")) {
                        myIsPhased = false;
                        if (values[i].contains("|")) {
                            throw new IllegalStateException("ParseGVCF: genotypes (GT) can't have both / and |");
                        }
                        alleles = values[i].split("/");
                        if (alleles[0].equals(alleles[1])) {
                            isHomoBasedOnGenotype = true;
                        }
                        if (!alleles[0].equals("0") || !alleles[1].equals("0")) {
                            isAlt = true;
                        }
                    } else if (values[i].contains("|")) {
                        myIsPhased = true;
                        alleles = values[i].split("|");
                        if (alleles[0].equals(alleles[1])) {
                            isHomoBasedOnGenotype = true;
                        }
                        if (!alleles[0].equals("0") || !alleles[1].equals("0")) {
                            isAlt = true;
                        }
                    } else {
                        myPloidy = 1;
                        if (!values[i].equals("0")) {
                            isAlt = true;
                        }
                        alleles = new String[]{values[i]};
                    }

                    ImmutableList.Builder<Integer> genotypeIndices = new ImmutableList.Builder<>();
                    ImmutableList.Builder<String> genotypes = new ImmutableList.Builder<>();

                    for (String current : alleles) {
                        int alleleIndex = -1;
                        try {
                            alleleIndex = Integer.valueOf(current);
                        } catch (NumberFormatException ne) {
                            myLogger.error(myLineNum + ": " + myLine);
                            throw new IllegalStateException("ParseGVCF: GT index not a number.");
                        }
                        try {
                            if (alleleIndex == 0) {
                                genotypeIndices.add(0);
                                genotypes.add(myReference);
                            } else {
                                genotypeIndices.add(alleleIndex);
                                genotypes.add(myAlternates.get(alleleIndex - 1));
                            }
                        } catch (ArrayIndexOutOfBoundsException ae) {
                            myLogger.error(myLineNum + ": " + myLine);
                            throw new IllegalStateException("ParseGVCF: GT index not defined.");
                        }
                    }

                    myGenotypeIndices = genotypeIndices.build();
                    myGenotypes = genotypes.build();

                }
            }

            // Compare sum of AD values equals DP
            if (alleleDepthSum != 0 && alleleDepthSum != myDepth) {
                myLogger.error(myLineNum + ": " + myLine);
                throw new IllegalStateException("GVCFLine: DP: " + myDepth + " not equal to AD sum: " + alleleDepthSum);
            }

            if ((myPloidy >= 2) && (myIsHomozygous != isHomoBasedOnGenotype)) {
                myLogger.error("Line: " + myLine);
                throw new IllegalStateException("GVCFLine: Homozygous doesn't match based on DP and GT");
            }

            if (myDepth > 0) {
                if (myIsHomozygous) {
                    if (isAlt) {
                        // numHomozygousAlt += length;
                        if (myDepth == 1) {
                            // numHomozygousAlt1 += length;
                        } else if (myDepth == 2) {
                            // numHomozygousAlt2 += length;
                        } else {
                            // numHomozygousAlt3 += length;
                        }
                    } else {
                        // numHomozygous += length;
                    }
                } else {
                    // numHeterozygous += length;
                }
            }

        }

    }

    private static class ReadLines implements Runnable {

        private final String myFilename;
        private final BlockingQueue<Future<ProcessLines>> myQueue;
        private final ExecutorService myPool;

        public ReadLines(String filename, BlockingQueue<Future<ProcessLines>> queue, ExecutorService pool) {
            myFilename = filename;
            myQueue = queue;
            myPool = pool;
        }

        @Override
        public void run() {

            try (BufferedReader reader = Utils.getBufferedReader(myFilename)) {

                int lineNum = 1;
                List<String> temp = new ArrayList<>();
                String line = reader.readLine();

                while (line != null && line.startsWith("#")) {
                    temp.add(line);
                    line = reader.readLine();
                }
                myQueue.add(myPool.submit(new ProcessLines(temp)));
                lineNum += temp.size();

                int count = 0;
                temp = new ArrayList<>();
                while (line != null) {
                    temp.add(line);
                    count++;
                    if (count == NUM_LINES_PER_BLOCK) {
                        myQueue.add(myPool.submit(new ProcessLines(lineNum, temp)));
                        temp = new ArrayList<>();
                        count = 0;
                        lineNum += NUM_LINES_PER_BLOCK;
                    }
                    line = reader.readLine();
                }

                if (!temp.isEmpty()) {
                    String[] lines = new String[temp.size()];
                    myQueue.add(myPool.submit(new ProcessLines(lineNum, temp)));
                }

                myQueue.add(myPool.submit(new ProcessLines()));

            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("ParseGVCF: ReadLines: problem reading file: " + myFilename);
            }

        }

    }

}
