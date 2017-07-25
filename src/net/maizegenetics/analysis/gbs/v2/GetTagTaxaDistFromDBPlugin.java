/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagDataSQLite;
import net.maizegenetics.dna.tag.TagDataWriter;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;

/**
 * This plugin takes a GBSv2 database as input, queries for the tags
 * and their taxa distribution, and creates a tab-delimited file of tag
 * and taxa-distribution.  IT can be used for verifying data in the db.
 * 
 * @author lcj34
 *
 */
public class GetTagTaxaDistFromDBPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(GetTagSequenceFromDBPlugin.class);

    private PluginParameter<String> myDBFile = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags").build();
    private PluginParameter<String> myOutputFile = new PluginParameter.Builder<String>("o", null, String.class).guiName("Output File").required(true).outFile()
            .description("Output txt file that can be imported to Excel").build();

    public GetTagTaxaDistFromDBPlugin() {
        super(null, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public GetTagTaxaDistFromDBPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

 
    @Override
    public DataSet processData(DataSet input) {
        TagDataWriter tdw = new TagDataSQLite(inputDB());
        try {
           //BufferedWriter fileWriter = null;
            BufferedWriter fileWriter = new BufferedWriter(new FileWriter(outputFile()));
            StringBuilder strB = new StringBuilder();
            TaxaList taxaList = tdw.getTaxaList(); // Used for printing taxon column headers
            strB.append("Tag");
            taxaList.stream().forEach(item -> { // column names are the taxon names
                strB.append("\t");
                strB.append(item.getName());
            });
            strB.append("\n");
            fileWriter.write(strB.toString());
            strB.setLength(0);
       
            Set<Tag> myTags = tdw.getTags();
            int tagcount = 0;
            for (Tag myTag: myTags) {
                tagcount++;
                // get dist for each taxa
                TaxaDistribution tagTD = tdw.getTaxaDistribution(myTag);
                if (tagTD == null) {
                    System.out.println("GetTagTaxaDist: got null tagTD at tagcount " + tagcount);
                    return null;
                }
                int[] depths = tagTD.depths(); // gives us the depths for each taxon
                strB.append(myTag.sequence());
                for (int idx = 0; idx < depths.length; idx++) {
                    strB.append("\t"); 
                    strB.append(depths[idx]);  // add tag depth                     
                }
                strB.append("\n"); // end of line - start next tag
                fileWriter.write(strB.toString());
                strB.setLength(0);
            }

            fileWriter.close();
            ((TagDataSQLite)tdw).close();  
            myLogger.info("TagsTaxaDistToTabDelim: Finished writing TaxaDistribution \n");
        } catch (Exception exc) {
            myLogger.error("GetTagTaxaDistFromDBPlugin: caught error " + exc);
            exc.printStackTrace();
        }
        return null;       
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
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
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
    public GetTagTaxaDistFromDBPlugin inputDB(String value) {
        myDBFile = new PluginParameter<>(myDBFile, value);
        return this;
    }

    /**
     * Output tab-delimited file showing tag sequences found in DB.
     *
     * @return Output File
     */
    public String outputFile() {
        return myOutputFile.value();
    }

    /**
     * Set Output File. This is tab delimited file with a list
     * of tag sequences found in the database (either a single tag
     * if the user is checking for a particular tag's presence, or
     * a list of all tags found in the db (if no user tag specified).
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public GetTagTaxaDistFromDBPlugin outputFile(String value) {
        myOutputFile = new PluginParameter<>(myOutputFile, value);
        return this;
    }

}
