/**
 * 
 */
package net.maizegenetics.dna.tag;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Tuple;

/**
 * @author lcj34
 *
 */
public interface RepGenData {
    /**
     * Returns the TaxaDistribution for a given taxa.  If distribution data is not available, null is returned.
     */
    TaxaDistribution getTaxaDistribution(Tag tag);

    /**
     * Set of all tags
     */
    Set<Tag> getTags();

    /**
     * Map of all tags with associated names
     */
    Map<Tag,String> getTagsNameMap();

    
    /**
     * Get SNPs for specified chromosomes
     * @param startChr chromosome number
     * @param endChr chromosome number
     * @return  PostionList of SNP positions for the requested chromosomes
     */
    PositionList getSNPPositionsForChromosomes(Integer startChr,Integer endChr);

    /**
     * Returns the list of tags present of a taxon.  Note this could be a very compute intensive request.
     */
    Set<Tag> getTagsForTaxon (Taxon taxon);


    /*
     * Return a map of Tags (keys) and depth (value) for a given taxon and SNP position.
     */
    Map<Tag, Integer> getTagDepth(Taxon taxon, Position position);

    /*
     * Return a map of Tags (keys) and depth (value) greater or equal to a minimum depth
     */
    Map<Tag, Integer> getTagsWithDepth(int minimumDepth);


     
    /*
    Map of alignment approaches names (key) and their protocols (value)
     */
    Map<String,String> getTagAlignmentApproaches();

    /**
     * Map of positions and with associated map of Tags and their taxa distribution and their alignment direction.
     * Warning:  This can be a very large data structure for entire chromosomes. Only the best positions are returned.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @return  Map of maps for position (key) to Map of Tag(key) to the Tuple(Direction,TaxaDistribution)(Value)
     */
    Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition);
    /**
     * For a given snp position,
     * Returns a Map of Allele with its associated Tag/TaxaDistribution 
     * Differs from getCutPositionTagTaxaMap() in that it now specifies on which
     * strand the returned tags must appear.
     * @return Map or null if not available.
     */
    Map<Position, Map<Tag, TaxaDistribution>> getCutPosForStrandTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition, boolean strand);


    /**
     * Returns the taxa list associated with taxa distribution
     * @return TaxaList or null if not available.
     */
    TaxaList getTaxaList();
        
    /**
     * Return all chromosomes stored in database.
     * Returns a list of distinct chromsome objects from the cutPosition table
     * @return List<Chromosome> or null if none found
     */
     List<Chromosome> getChromosomesFromCutPositions();

     /**
      * Return all tag/taxadistribution stored in database.
      * Returns a Hashmap of tag.taxadistributions
      * @return Map<Tag, TaxaDistribution>  - map is empty if db table is empty
      */
     Map<Tag, TaxaDistribution> getAllTagsTaxaMap();

     /**
      * Get all of the physical map positions associated with the ref tags. 
      * @return iterator of annotated positions
      */
    PositionList getPhysicalMapPositions();

    /**
     * Get the unique list of positions for ref tags for a specific chromosome within a range.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @return List of positions
     */
    PositionList getPhysicalMapPositions(Chromosome chromosome, int firstPosition, int lastPosition);
    
    /**
     * Set of all reference tags
     */
    Set<RefTagData> getRefTags();
    

    /**
     * Method returns a map of all tag (non-ref) alignments for specified non-ref tags.
     * All tags on the list must be non-ref tags
     * as the tagID is pulled from  the tagTagIDMap.
     *
     * @return
     */
    Multimap<Tag,AlignmentInfo> getTagAlignmentsForTags(List<Tag> tags, int minscore);

    /**
     * Grab correlations for specified tags
     * @param  Multimap of tags whose correlations will be pulled
     */

    Multimap<Tag,TagCorrelationInfo> getCorrelationsForTags(List<Tag> tags);

    /**
     * Method to grab the tag sequence, chrom, position and taxa distributions
     * from the db.  This returns values for all tags that have an entry in
     * the tagMapping table.
     * @return
     */

    ImmutableMultimap<Tag, Tuple<Position, TaxaDistribution>> getPositionTaxaDistForTag();

    /**
     * Method grabs data from the tag_tag_stats table for the specified method.
     * @param method_id Get data for entries with this method ID
     * @return
     */
    ImmutableMultimap<Tuple<Tag, Tag>, Tuple<Float, String>> getTag_tag_stats(int method_id);
}
