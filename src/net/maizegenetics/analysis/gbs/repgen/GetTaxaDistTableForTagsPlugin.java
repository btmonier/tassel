/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

/**
 * This plugin takes a rAmpSeq database file and pulls all the chrom/pos and taxa distribution
 * data for each tag that aligns perfectly to a reference tag.
 * 
 * Output is written to 4 tab delimited files:
 *   - One is a summary file that give stats for each tag telling the number of Samples (taxa), the number of Samples
 *     that contain this tag, and the total depth for the tag.
 *   - The second file is a tag/taxa-distribution file.  It shows tag counts for each tag in each taxa
 *   
 *   The last 2 files look like: 
 *   Each file has a header line:  Tag Chromosome Posiiton Sample1 .... SampleN 
 *   -one that shows single-mapping tags
 *   -second that shows tags that mapped to multiple places on the genome.
 *   -The output directory should end with a "/" followed by any prefix to occur before the set file name as below
 *   The files are named outputPrefix() + singleMappedTags.txt
 *                       outputPrefix() + multipleMappedTags.txt
 * @author lcj34
 *
 */
public class GetTaxaDistTableForTagsPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(GetTaxaDistTableForTagsPlugin.class);
    private PluginParameter<String> inputDB = new PluginParameter.Builder<>("db", null, String.class).guiName("Input Database File").required(true).inFile()
            .description("Input Database File").build();
    private PluginParameter<String> outputPrefix = new PluginParameter.Builder<>("o", null, String.class).guiName("").required(true)
            .description("Directory and any prefix for writing output.\nThere are 4 files, one for tags aligning to a single position, and one for tags aigning to multiple positions.\nIn addition a tag-taxadistribution file is created as well as a tag-taxadistribution summary file.").build();
    
    public GetTaxaDistTableForTagsPlugin() {
        super(null, false);
    }

    public GetTaxaDistTableForTagsPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public GetTaxaDistTableForTagsPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        System.out.println("BEGIN CreateHistogramsForTagMetrics: get tags/refTags from DB");
        RepGenDataWriter repGenData=new RepGenSQLite(inputDB());
        
        // FIrst get all tag/taxadist .  This will include tags
        // that are NOT ref tags as well as ones that are.
        BufferedWriter fileWriter = null;
        BufferedWriter summaryWriter = null;
        StringBuilder tagTaxaDistSB = new StringBuilder();
        StringBuilder summarySB = new StringBuilder();
        List<Taxon> taxa = repGenData.getTaxaList();
        
        // Create tagTaxaDist header line
        tagTaxaDistSB.append("Tag");
        taxa.stream().forEach(item -> { // column names are the taxon names
            tagTaxaDistSB.append("\t");
            tagTaxaDistSB.append(item.getName());
        });
        tagTaxaDistSB.append("\n");
   
        String summaryHeader = "Tag\tNumber of Samples\tSamples With Tag\tTotal Tag Depth\n";
        summarySB.append(summaryHeader);
        
        int numSamples = taxa.size();
        Set<Tag> myTags = repGenData.getTags();
        int tagcount = 0;
        int experimentalTags = 0; // these are all tags with taxadist (ie, the tags loaded from fastq)
        int numRefTagsNoTaxaDist = 0;
        for (Tag myTag: myTags) {
            tagcount++;
            // get dist for each taxa
            TaxaDistribution tagTD = repGenData.getTaxaDistribution(myTag);
            if (tagTD == null) {
                numRefTagsNoTaxaDist++; // reftags have no taxadist unless they map to non-refTag
                continue;
            }
            experimentalTags++;
            int[] depths = tagTD.depths(); // gives us the depths for each taxon
            tagTaxaDistSB.append(myTag.sequence());
 
            int totalDepths = 0;
            int taxaCount = 0;
            for (int idx = 0; idx < depths.length; idx++) {
                tagTaxaDistSB.append("\t"); 
                tagTaxaDistSB.append(depths[idx]);  // add tag depth  
                totalDepths += depths[idx];
                if (depths[idx] > 0) taxaCount++;
            }
            tagTaxaDistSB.append("\n"); // end of line - start next tag
            String summarydata = myTag.sequence() + "\t" + numSamples + "\t" + taxaCount + "\t" + totalDepths + "\n";
            summarySB.append(summarydata.toString());
        }
        String tagTaxafile = outputPrefix() + "allTagTaxaDist.txt";
        String summaryFile = outputPrefix() + "tagTaxaSummaryData.txt";
        try {  
            fileWriter = Utils.getBufferedWriter(tagTaxafile);
            fileWriter.write(tagTaxaDistSB.toString());
            fileWriter.close();
            summaryWriter = Utils.getBufferedWriter(summaryFile);
            summaryWriter.write(summarySB.toString());
            summaryWriter.close();
        }
        catch(IOException e) {
            myLogger.error("Caught Exception writing to outputFile " + tagTaxafile);
            System.out.println(e);
        }
        System.out.println("Finsihed processing allTagTaxaDist - number of tags: " + tagcount + ", number without taxaDist: " + numRefTagsNoTaxaDist);
        
        // Now process tags that have been mapped to reference positions.
        // This data is reftags.  Some non-refTags have the same sequence as refTags.
        Multimap<Tag,Tuple<Position,TaxaDistribution>> tagDataMap = HashMultimap.<Tag,Tuple<Position,TaxaDistribution>>create();
        tagDataMap = repGenData.getPositionTaxaDistForTag();
        
        System.out.println("Size of tagDataMap keyset: " + tagDataMap.keySet().size() + " size of all entries: " + tagDataMap.size());
               
        String singleMapFile = outputPrefix() + "singleMappedRefTags.txt";
        String multipleMapFile = outputPrefix() + "multipleMappedRefTags.txt";
        BufferedWriter sbw = Utils.getBufferedWriter(singleMapFile);
        BufferedWriter mbw = Utils.getBufferedWriter(multipleMapFile);
        
        // Create and write header line.
        StringBuilder sb  = new StringBuilder();
        String initial = "Tag\tChromosome\tPosition";
        sb.append(initial);
        taxa.stream().forEach(item -> { // column names are the taxon names
            sb.append("\t");
            sb.append(item.getName());
        });

        sb.append("\n");
        int tagsWithMappingAndTaxa = 0;
        int dataLines1 = 0;
        int dataLinesM = 0;
        try {
            
            sbw.write(sb.toString());
            mbw.write(sb.toString());  
            for (Tag tag : tagDataMap.keySet()) {
                tagsWithMappingAndTaxa++;
                StringBuilder sb2 = new StringBuilder();
                Collection<Tuple<Position, TaxaDistribution>>  values = tagDataMap.get(tag);
                List<Tuple<Position,TaxaDistribution>> data = new ArrayList<Tuple<Position,TaxaDistribution>>(values);
                if (data.size() > 1) {
                    // add to multiple mapping list
                    for (Tuple<Position,TaxaDistribution> posTD : data) {
                        
                        String chrom = posTD.x.getChromosome().getName();
                        int pos = posTD.x.getPosition();
                        String entryString = tag.sequence() + "\t" + chrom + "\t" + pos;
                        sb2.append(entryString);
                        int[] depths = data.get(0).y.depths(); // gives us the depths for each taxon
                        
                        for (int idx = 0; idx < depths.length; idx++) {
                            sb2.append("\t"); 
                            sb2.append(depths[idx]);  // add tag depth                     
                        }
                        sb2.append("\n"); // end of line - start next tag
                        mbw.write(sb2.toString());
                        dataLinesM++;
                        sb2.setLength(0);
                    }                   
                } else {
                    // Processing for files with single mapping
                    String chrom = data.get(0).x.getChromosome().getName();
                    int pos = data.get(0).x.getPosition();
                    String entryString = tag.sequence() + "\t" + chrom + "\t" + pos;
                    sb2.append(entryString);
                    int[] depths = data.get(0).y.depths(); // gives us the depths for each taxon
                    
                    for (int idx = 0; idx < depths.length; idx++) {
                        sb2.append("\t"); 
                        sb2.append(depths[idx]);  // add tag depth                     
                    }
                    sb2.append("\n"); // end of line - start next tag
                    sbw.write(sb2.toString());
                    dataLines1++;
                    sb2.setLength(0);
                }               
            }
            sbw.close();
            mbw.close();
            ((RepGenSQLite)repGenData).close();
        } catch (Exception exc) {
            exc.printStackTrace();
        }
 
        System.out.println("\nTotal tags with physical mapping plus taxa distribution: " + tagsWithMappingAndTaxa + ", singleEntry: " + dataLines1 + ", multiple entries: " + dataLinesM);
        int numRefTags = numRefTagsNoTaxaDist + tagsWithMappingAndTaxa;
        
        System.out.println("\nNum tags in DB: " + tagcount + ", number of refTags in db: " + numRefTags + ", number experimental tags: " + experimentalTags);
        System.out.println("Num tags perfect mapping to refTag: " + tagsWithMappingAndTaxa);
        float perfectToExperimental = (float)tagsWithMappingAndTaxa/experimentalTags;
        System.out.println("Proportion perfectMapping to experimental tags: " + perfectToExperimental);
        System.out.println("\nProcess took " + (System.nanoTime() - totalTime)/1e9 + " seconds.");
 
        return null;
    }

//    public static void main (String[] args) {
// 
//    }
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
     * Input Database File
     *
     * @return Input Database File
     */
    public String inputDB() {
        return inputDB.value();
    }

    /**
     * Set Input Database File. Input Database File
     *
     * @param value Input Database File
     *
     * @return this plugin
     */
    public GetTaxaDistTableForTagsPlugin inputDB(String value) {
        inputDB = new PluginParameter<>(inputDB, value);
        return this;
    }

    /**
     * Directory and prefix for writing output.
     * There are 2 files, one for tags aligning to a single
     * position, and one for tags aigning to multiple positions.
     *
     * @return O
     */
    public String outputPrefix() {
        return outputPrefix.value();
    }

    /**
     * Set O. Directory and prefix for writing output.
     * There are 2 files, one for tags aligning to a single
     * position, and one for tags aigning to multiple positions.
     *
     * @param value O
     *
     * @return this plugin
     */
    public GetTaxaDistTableForTagsPlugin outputPrefix(String value) {
        outputPrefix = new PluginParameter<>(outputPrefix, value);
        return this;
    }
}
