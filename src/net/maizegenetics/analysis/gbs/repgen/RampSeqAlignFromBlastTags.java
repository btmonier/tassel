/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.log4j.Logger;
import org.jfree.data.general.Dataset;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * This class takes a rAmpSeq database populated with tags,a reference Genome, a filtered
 * BLAST file output and primers, then:
 *   (1) creates  reference tags using the blast alignment output
 *   (2) runs tag/tag alignment for each tag in the DB,
 *   (3) runs tag/refTag alignment for each tag against each ref tag
 *   (4) runs refTag/refTag alignment for each refTag in the DB>
 *   (5) all alignments are stored in the db tagAlignments table.
 * 
 * Alignments are performed and stored in groups to prevent overwhelming
 * the DB with massive load commands.
 * 
 * Blast was run on CBSU using these parameters:
 * Make the reference files:
 *    makeblastdb -dbtype nucl -in <agpv4 referencene genome fasta> -parse_seqids -out maizeAGPV4.db
 *    
 * Run blast using maizeAGPV4.db from command above:
 *   blastn -num_threads 24 -db maizeAGPV4.db -query anp68R1Tags.fasta -evalue 1e-60 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out anp68TagsR1Result/blastANP68_R1.txt
 *   
 * The blast output file was filtered using the 3 commands below.  The first filters
 * identity down to 98 %,  the second gets alignment lengths that were at least 148,
 * the 3rd filters it down to just the chrom, start, end positions:
 * 
 * awk '$3 >= 98.000 {print $0}' blastANP68_R1.txt > blastANP68_R1_98per.txt
 * awk '$4 >= 148 {print $0}' blastANP68_R1_98per.txt > blastANP68_R1_98per_148align.txt
 * awk {'printf ("%s\t%s\t%s\n", $2, $9, $10)'} blastANP68_R1_98per_148align.txt > blastANP68_R1_98per_148align_3cols.txt
 * 
 * It is the last file from awk, blastANP68_R1_98per_148align_3cols.txt, that is given as a parameter here.
 * 
 * @author lcj34
 *
 */
public class RampSeqAlignFromBlastTags extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(RepGenPhase2AlignerPlugin.class);
    
    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<String> refGenome = new PluginParameter.Builder<String>("ref", null, String.class).guiName("Reference Genome File").required(true)
            .description("Referemce Genome File for aligning against ").build();
    private PluginParameter<Integer> minIdentity = new PluginParameter.Builder<Integer>("min_identity", 98, Integer.class).guiName("Minimum identity")
            .description("Parameter to filter the BLAST results by minimum identity").build();
    private PluginParameter<String> blastFile = new PluginParameter.Builder<String>("blastFile", null, String.class).guiName("Blast File").required(true).inFile()
            .description("Tab delimited Blast file output with NO header line, that contains only the data from columns for chrom, start postion, end position. \nThis data should be filtered to contain only entries whose identiy value was 98% or greater.").build();
    
    static GenomeSequence myRefSequence = null;

    public RampSeqAlignFromBlastTags() {
        super(null, false);
    }

    public RampSeqAlignFromBlastTags(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RampSeqAlignFromBlastTags(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public void postProcessParameters() {

        if (myDBFile.isEmpty() || !Files.exists(Paths.get(dBFile()))) {
            throw new IllegalArgumentException("RampSeqAlignFromBlastTags: postProcessParameters: Input DB not set or found");
        }
        if (!refGenome.isEmpty()) {
            myRefSequence = GenomeSequenceBuilder.instance(refGenome());
        } else {
            throw new IllegalArgumentException("RampSeqAlignFromBlastTags: postProcessParameters: reference genome not set or found");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
        RepGenDataWriter repGenData=null;
 
        try {           
            System.out.println("RampSeqAlignFromBlastTags:processData begin"); 
            repGenData=new RepGenSQLite(dBFile());
            time = System.nanoTime();
            repGenData.addMappingApproach("ReferenceBLAST");
            Multimap<Tag,Position> reftagPosMap = createRefTagsFromBlast(blastFile());
            repGenData.putTagPositionMapping(reftagPosMap,refGenome(), "ReferenceBLAST");
            ((RepGenSQLite)repGenData).close();

            myLogger.info("Finished RampSeqAlignFromBlastTags\n");
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        
        
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return DataSet.getDataSet(repGenData);
    }

    public Multimap<Tag,Position> createRefTagsFromBlast(String chromStartEndFile) {
        System.out.println("Begin createRefTagsFromBlast");
        final Multimap<Tag,Position> reftagPosMap=HashMultimap.create();
        try {
            Files.lines(Paths.get(chromStartEndFile))
                    .map(line -> line.split("\t"))
                    .filter(token -> Double.parseDouble(token[2])>minIdentity())
                    .forEach(token -> {
                        Chromosome chr=new Chromosome(token[1]);
                        int subjectStart=Integer.parseInt(token[8]);
                        int subjectEnd=Integer.parseInt(token[9]);
                        byte strand=(subjectStart<subjectEnd)?Position.STRAND_PLUS:Position.STRAND_MINUS;
                        Position position=new GeneralPosition.Builder(chr,subjectStart).strand(strand).build();
                        Tag tag;
                        if(strand==Position.STRAND_PLUS) {
                            tag=TagBuilder.instance(NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.chromosomeSequence(chr,subjectStart, subjectEnd))).build();
                        } else {
                            tag=TagBuilder.instance(NucleotideAlignmentConstants.nucleotideBytetoString(myRefSequence.chromosomeSequence(chr,subjectEnd, subjectStart))).build();
                            tag=TagBuilder.reverseComplement(tag).build();
                        }
                        System.out.println(tag.toCSVString()+" position:"+position.toString());
                        reftagPosMap.put(tag,position);
                    });

            System.out.println("Read " + reftagPosMap.size() + " lines from BLAST input file.");
        } catch (IOException exc) {
            exc.printStackTrace();
        }
        return reftagPosMap;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    public static void main(String[] args) {
        GeneratePluginCode.generate(RampSeqAlignFromBlastTags.class);
    }

     
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    public RepGenDataWriter runPlugin(DataSet input) {
        return (RepGenDataWriter) performFunction(input).getData(0).getData();
    }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String dBFile() {
        return myDBFile.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags dBFile(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Referemce Genome File for aligning against
     *
     * @return Reference Genome File
     */
    public String refGenome() {
        return refGenome.value();
    }

    /**
     * Set Reference Genome File. Referemce Genome File for
     * aligning against
     *
     * @param value Reference Genome File
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags refGenome(String value) {
        refGenome = new PluginParameter<>(refGenome, value);
        return this;
    }

    /**
     * Parameter to filter the BLAST results by minimum identity
     *
     * @return Minimum identity
     */
    public Integer minIdentity() {
        return minIdentity.value();
    }

    /**
     * Set Minimum identity. Parameter to filter the BLAST
     * results by minimum identity
     *
     * @param value Minimum identity
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags minIdentity(Integer value) {
        minIdentity = new PluginParameter<>(minIdentity, value);
        return this;
    }

    /**
     * Tab delimited Blast file output with NO header line,
     * that contains only the data from columns for chrom,
     * start postion, end position.
     * This data should be filtered to contain only entries
     * whose identiy value was 98% or greater.
     *
     * @return Blast File
     */
    public String blastFile() {
        return blastFile.value();
    }

    /**
     * Set Blast File. Tab delimited Blast file output with
     * NO header line, that contains only the data from columns
     * for chrom, start postion, end position.
     * This data should be filtered to contain only entries
     * whose identiy value was 98% or greater.
     *
     * @param value Blast File
     *
     * @return this plugin
     */
    public RampSeqAlignFromBlastTags blastFile(String value) {
        blastFile = new PluginParameter<>(blastFile, value);
        return this;
    }
}
