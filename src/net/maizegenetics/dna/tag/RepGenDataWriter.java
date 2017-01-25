/**
 * 
 */
package net.maizegenetics.dna.tag;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.Tuple;

/**
 * @author lcj34
 *
 */
public interface RepGenDataWriter extends RepGenData {
    /**
     * Add a tag to list of known tags
     * @return true if this set did not already contain the specified element
     * @param tags
     * @param tagInstanceAverageQS map showing number of tags, and average quality score for each tag
     */
    boolean putAllTag(Set<Tag> tags,Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS);

    /**
     * Add a tag to list of known tags with their associated names (good for named contigs)
     * @return true if this set did not already contain the specified element
     * @param tagNameMap
     */
    boolean putAllNamesTag(Map<Tag,String> tagNameMap);

    /**
     * Associates a map full of the specified Tag (key) with the specified TaxaDistribution (value).
     * If there was a prior association is it replaced, as all pairs are unique.
     */
    void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap);


    /**
     * Associates the specified Reference Tag with the specified site Position (value). 
     * @param tagAnnotatedPositionMap Map of specific tag with chrom/positions specified.
     *       
     */
    void putRefTagMapping(Multimap<Tag, Position> tagAnnotatedPositionMap, String refGenome);
    
    /**
     * Stores the Smith Waterman score from2 tag alignments. 
     * tag2 chrom/pos comes from the AlignmentInfo object.  tag1 chrom/pos are separate parameters
     * @param tagAlignInfoMap Map of specific tag to tag2 alignment data
     * 
     */
    void putTagTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap);
    

    /*
    Set the specified Tag and Position combination to best, and all others were set to false.
     */
    void setTagAlignmentBest(Tag tag, Position position, boolean isBest);

    /*
    Adds a new Alignment approach name to a detailed protocol.
     */
    boolean putTagAlignmentApproach(String tagAlignmentName, String protocol);

    /*Sets the taxaList for given set of Taxa, this is the order in which the taxa distribution is recorded*/
    void putTaxaList(TaxaList taxaList);
    
    /**
     * Removes all data from the tagtaxadistribution table. 
     * This should be called from GBSSeqToTagDBPlugin when a user requests an "append"
     * option (ie, db exists, and user opted NOT to clear it).  AFter we grab the
     * existing data it is cleared.  It will be re-entered when GBSSeqToTagDBPlugin completes.
     * 
     */
    void clearTagTaxaDistributionData();
    
    /**
     * Removes all data from the DB that was added from SAMToGBSDBPlugin call.
     * The tables cleared are  CutPosition and TagCutPosition 
     */
    void clearAlignmentData();
    
    /**
     * Removes all data from the DB that was added from the DiscoverySNPCallerPluginV2
     * The tables cleared are  Allele, TagAllele and SNPPosition 
     */
    void clearDiscoveryData();
    
    /**
     * Removes all data from the snpQuality table
     */
    void clearSNPQualityData();

    /**
     * Adds a mapping approach strategy to the mappingApproach table
     * @param name
     */
    void addMappingApproach(String name);

    /**
     * Adds a mapping approach strategy to the mappingApproach table
     * @param name
     * @throws SQLException
     */
    void addReferenceGenome(String refGenome);
    
    /**
     * Adds entries to the tagCorrelations table for tag-tag correlation data 
     * @param tagCorrelationMap holds correltaions info for each tag/taxa-depth to tag/taxa-depth vector pair
     * @throws SQLException
     */
    void putTagTagCorrelationMatrix(Multimap<Tag,TagCorrelationInfo> tagCorrelationMap);
}
