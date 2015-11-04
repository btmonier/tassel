/*
 *  FilterSiteBuilderPlugin
 * 
 *  Created on Jun 30, 2014
 */
package net.maizegenetics.analysis.filter;

import com.google.common.collect.Range;

import java.awt.Frame;
import java.lang.reflect.Field;
import java.net.URL;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.ImageIcon;

import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.FilterSite;
import net.maizegenetics.dna.snp.FilterSite.SITE_RANGE_FILTER_TYPES;
import net.maizegenetics.dna.snp.FilterSite.FILTER_SITES_ATTRIBUTES;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FilterSiteBuilderPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterSiteBuilderPlugin.class);

    private PluginParameter<String> myFilterName = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.filterName.name(), "Filter", String.class)
            .description("Filter Name")
            .build();
    private PluginParameter<Integer> mySiteMinCount = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.siteMinCount.name(), 0, Integer.class)
            .range(Range.atLeast(0))
            .description("Site Minimum Count of Alleles not Unknown")
            .build();
    private PluginParameter<Double> mySiteMinAlleleFreq = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.siteMinAlleleFreq.name(), 0.0, Double.class)
            .range(Range.closed(0.0, 1.0))
            .description("Site Minimum Minor Allele Frequency")
            .build();
    private PluginParameter<Double> mySiteMaxAlleleFreq = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.siteMaxAlleleFreq.name(), 1.0, Double.class)
            .range(Range.closed(0.0, 1.0))
            .description("Site Maximum Minor Allele Frequency")
            .build();
    private PluginParameter<Boolean> myRemoveMinorSNPStates = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.removeMinorSNPStates.name(), false, Boolean.class)
            .description("")
            .build();
    private PluginParameter<SITE_RANGE_FILTER_TYPES> mySiteFilter = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.siteRangeFilterType.name(), SITE_RANGE_FILTER_TYPES.NONE, SITE_RANGE_FILTER_TYPES.class)
            .description("True if filtering by site numbers. False if filtering by chromosome and position")
            .build();
    private PluginParameter<Integer> myStartSite = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.startSite.name(), 0, Integer.class)
            .range(Range.atLeast(0))
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.SITES)
            .description("")
            .build();
    private PluginParameter<Integer> myEndSite = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.endSite.name(), 0, Integer.class)
            .range(Range.atLeast(0))
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.SITES)
            .description("")
            .build();
    private PluginParameter<Chromosome> myStartChr = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.startChr.name(), null, Chromosome.class)
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.POSITIONS)
            .description("")
            .build();
    private PluginParameter<Integer> myStartPos = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.startPos.name(), 0, Integer.class)
            .range(Range.atLeast(0))
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.POSITIONS)
            .description("")
            .build();
    private PluginParameter<Chromosome> myEndChr = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.endChr.name(), null, Chromosome.class)
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.POSITIONS)
            .description("")
            .build();
    private PluginParameter<Integer> myEndPos = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.endPos.name(), 0, Integer.class)
            .range(Range.atLeast(0))
            .dependentOnParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.POSITIONS)
            .description("")
            .build();
    private PluginParameter<Boolean> myIncludeSites = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.includeSites.name(), true, Boolean.class)
            .description("")
            .build();
    private PluginParameter<PositionList> myPositionList = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.positionList.name(), null, PositionList.class)
            .positionList()
            .description("Filter based on position list.")
            .build();
    private PluginParameter<String> mySiteNamesList = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.siteNames.name(), null, String.class)
            .siteNameList()
            .dependentOnParameter(myPositionList, POSITION_LIST_NONE)
            .description("Filter based on site names.")
            .build();
    private PluginParameter<String> myBedFile = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.bedFile.name(), null, String.class)
            .inFile()
            .dependentOnParameter(myPositionList, POSITION_LIST_NONE)
            .description("Filter based on BED file.")
            .build();
    private PluginParameter<String> myChrPosFile = new PluginParameter.Builder<>(FILTER_SITES_ATTRIBUTES.chrPosFile.name(), null, String.class)
            .inFile()
            .description("Filter based on list of chromsome / position in file.")
            .build();

    public FilterSiteBuilderPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public FilterSiteBuilderPlugin() {
        this(null, false);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        if (input == null) {
            return;
        }

        if (isInteractive()) {

            List<Datum> filterList = input.getDataOfType(FilterSite.class);
            if (filterList.size() == 1) {
                FilterSite filter = (FilterSite) filterList.get(0).getData();
                useFilterValues(filter);
                return;
            }

            List<Datum> genotypeTableList = input.getDataOfType(GenotypeTable.class);
            if (genotypeTableList.size() == 1) {
                GenotypeTable genotypeTable = (GenotypeTable) genotypeTableList.get(0).getData();
                int lastSite = genotypeTable.numberOfSites() - 1;
                setParameter(myEndSite.cmdLineName(), lastSite);
                PositionList positions = genotypeTable.positions();
                setParameter(myStartChr, positions.chromosome(0));
                setParameter(myStartPos, positions.get(0).getPosition());
                setParameter(myEndChr, positions.chromosome(lastSite));
                setParameter(myEndPos, positions.get(lastSite).getPosition());
                return;
            }

            List<Datum> positionLists = input.getDataOfType(PositionList.class);
            if (positionLists.size() == 1) {
                PositionList positions = (PositionList) positionLists.get(0).getData();
                int lastSite = positions.numberOfSites() - 1;
                setParameter(myStartChr, positions.chromosome(0));
                setParameter(myStartPos, positions.get(0).getPosition());
                setParameter(myEndChr, positions.chromosome(lastSite));
                setParameter(myEndPos, positions.get(lastSite).getPosition());
                return;
            }

        } else {

            if ((!myStartSite.value().equals(myStartSite.defaultValue()) || (!myEndSite.value().equals(myEndSite.defaultValue())))) {
                setParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.SITES);
            }

            if ((!myStartPos.value().equals(myStartPos.defaultValue()))
                    || (myStartChr.value() != null)
                    || (!myEndPos.value().equals(myEndPos.defaultValue()))
                    || (myEndChr.value() != null)) {

                setParameter(mySiteFilter, SITE_RANGE_FILTER_TYPES.POSITIONS);

                List<Datum> genotypeTableList = input.getDataOfType(GenotypeTable.class);
                if (genotypeTableList.size() == 1) {
                    GenotypeTable genotype = (GenotypeTable) genotypeTableList.get(0).getData();
                    if (myStartPos.value().equals(myStartPos.defaultValue()) && myStartChr.value() != null) {
                        int[] firstLast = genotype.firstLastSiteOfChromosome(myStartChr.value());
                        setParameter(myStartPos, genotype.positions().get(firstLast[0]).getPosition());
                    }
                    if (myEndPos.value().equals(myEndPos.defaultValue()) && myEndChr.value() != null) {
                        int[] firstLast = genotype.firstLastSiteOfChromosome(myEndChr.value());
                        setParameter(myEndPos, genotype.positions().get(firstLast[1]).getPosition());
                    }
                }

            }

        }

    }

    @Override
    public DataSet processData(DataSet input) {

        Map<String, Object> values = new LinkedHashMap<>();
        for (Field field : getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(this);
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }
            if (current != null) {
                if (((current.value() != null)) && (!current.value().equals(current.defaultValue()))
                        || (current.cmdLineName().equals(FILTER_SITES_ATTRIBUTES.filterName.name()))) {
                    values.put(current.cmdLineName(), current.value());
                }
            }
        }

        if (siteFilter() == SITE_RANGE_FILTER_TYPES.SITES) {
            values.remove(FILTER_SITES_ATTRIBUTES.startChr.name());
            values.remove(FILTER_SITES_ATTRIBUTES.startPos.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endChr.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endPos.name());
            if (values.get(FILTER_SITES_ATTRIBUTES.startSite.name()) == null) {
                values.put(FILTER_SITES_ATTRIBUTES.startSite.name(), 0);
            }
        } else if (siteFilter() == SITE_RANGE_FILTER_TYPES.POSITIONS) {
            values.remove(FILTER_SITES_ATTRIBUTES.startSite.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endSite.name());
            if (values.get(FILTER_SITES_ATTRIBUTES.startPos.name()) == null) {
                values.put(FILTER_SITES_ATTRIBUTES.startPos.name(), 0);
            }
        } else if (siteFilter() == SITE_RANGE_FILTER_TYPES.NONE) {
            values.remove(FILTER_SITES_ATTRIBUTES.startSite.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endSite.name());
            values.remove(FILTER_SITES_ATTRIBUTES.startChr.name());
            values.remove(FILTER_SITES_ATTRIBUTES.startPos.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endChr.name());
            values.remove(FILTER_SITES_ATTRIBUTES.endPos.name());
        }

        Object sites = values.get(FILTER_SITES_ATTRIBUTES.siteNames.name());
        if (sites != null) {
            values.put(FILTER_SITES_ATTRIBUTES.siteNames.name(), getListFromString(sites.toString()));
        }

        List<Datum> result = new ArrayList<>();

        FilterSite filter = new FilterSite(values);

        List<Datum> genotypeTableList = input.getDataOfType(GenotypeTable.class);
        if (genotypeTableList.size() >= 1) {
            for (Datum datum : genotypeTableList) {
                GenotypeTable current = (GenotypeTable) datum.getData();
                GenotypeTable filteredGenotype = GenotypeTableUtils.filter(filter, current);
                if (filteredGenotype != current) {
                    Datum temp = new Datum(datum.getName() + "_" + filter.filterName(), filteredGenotype, null);
                    result.add(temp);
                    GenotypeSummaryPlugin.printSimpleSummary(temp);
                } else {
                    myLogger.info("processData: " + datum.getName() + " unchanged.");
                }
            }
        }

        result.add(new Datum(filter.filterName(), new FilterList(filter), null));

        return new DataSet(result, this);

    }

    private List<String> getListFromString(String str) {

        if ((str == null) || (str.length() == 0)) {
            return null;
        }
        String[] tokens = str.split(",");
        List<String> result = new ArrayList<>();
        for (String current : tokens) {
            current = current.trim();
            if (current.length() != 0) {
                result.add(current);
            }
        }
        return result;

    }

    private String getStringFromList(List<String> list) {

        if ((list == null) || (list.isEmpty())) {
            return null;
        }

        StringBuilder builder = new StringBuilder();
        boolean first = true;
        for (String current : list) {
            current = current.trim();
            if (current.length() != 0) {
                if (first) {
                    first = false;
                } else {
                    builder.append(",");
                }
                builder.append(current);
            }
        }
        return builder.toString();

    }

    public FilterSiteBuilderPlugin siteNamesList(List<String> value) {
        mySiteNamesList = new PluginParameter<>(mySiteNamesList, getStringFromList(value));
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = FilterSiteBuilderPlugin.class.getResource("/net/maizegenetics/analysis/images/FilterNew.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Filter Genotype Table Sites";
    }

    @Override
    public String getToolTipText() {
        return "Filter Genotype Table Sites";
    }

    public FilterSite build() {
        return (FilterSite) performFunction(null).getData(0).getData();
    }

    public FilterSiteBuilderPlugin useFilterValues(FilterSite filter) {
        setParametersToDefault();
        filter.attributes().entrySet().stream().forEach((attribute) -> {
            setParameter(attribute.getKey().name(), attribute.getValue());
        });
        return this;
    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FilterSiteBuilderPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public FilterSite runPlugin() {
        return (FilterSite) performFunction(null).getData(0).getData();
    }

    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getDataOfType(GenotypeTable.class).get(0).getData();
    }

    public GenotypeTable runPlugin(GenotypeTable input) {
        return (GenotypeTable) performFunction(DataSet.getDataSet(input)).getDataOfType(GenotypeTable.class).get(0).getData();
    }

    /**
     * Filter Name
     *
     * @return Filter Name
     */
    public String filterName() {
        return myFilterName.value();
    }

    /**
     * Set Filter Name. Filter Name
     *
     * @param value Filter Name
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin filterName(String value) {
        myFilterName = new PluginParameter<>(myFilterName, value);
        return this;
    }

    /**
     * Site Minimum Count of Alleles not Unknown
     *
     * @return Site Min Count
     */
    public Integer siteMinCount() {
        return mySiteMinCount.value();
    }

    /**
     * Set Site Min Count. Site Minimum Count of Alleles not Unknown
     *
     * @param value Site Min Count
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin siteMinCount(Integer value) {
        mySiteMinCount = new PluginParameter<>(mySiteMinCount, value);
        return this;
    }

    /**
     * Site Minimum Minor Allele Frequency
     *
     * @return Site Min Allele Freq
     */
    public Double siteMinAlleleFreq() {
        return mySiteMinAlleleFreq.value();
    }

    /**
     * Set Site Min Allele Freq. Site Minimum Minor Allele Frequency
     *
     * @param value Site Min Allele Freq
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin siteMinAlleleFreq(Double value) {
        mySiteMinAlleleFreq = new PluginParameter<>(mySiteMinAlleleFreq, value);
        return this;
    }

    /**
     * Site Maximum Minor Allele Frequency
     *
     * @return Site Max Allele Freq
     */
    public Double siteMaxAlleleFreq() {
        return mySiteMaxAlleleFreq.value();
    }

    /**
     * Set Site Max Allele Freq. Site Maximum Minor Allele Frequency
     *
     * @param value Site Max Allele Freq
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin siteMaxAlleleFreq(Double value) {
        mySiteMaxAlleleFreq = new PluginParameter<>(mySiteMaxAlleleFreq, value);
        return this;
    }

    /**
     * Remove Minor SNP States
     *
     * @return Remove Minor SNP States
     */
    public Boolean removeMinorSNPStates() {
        return myRemoveMinorSNPStates.value();
    }

    /**
     * Set Remove Minor SNP States. Remove Minor S N P States
     *
     * @param value Remove Minor S N P States
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin removeMinorSNPStates(Boolean value) {
        myRemoveMinorSNPStates = new PluginParameter<>(myRemoveMinorSNPStates, value);
        return this;
    }

    /**
     * True if filtering by site numbers. False if filtering by chromosome and
     * position
     *
     * @return Site Filter
     */
    public SITE_RANGE_FILTER_TYPES siteFilter() {
        return mySiteFilter.value();
    }

    /**
     * Set Site Filter. True if filtering by site numbers. False if filtering by
     * chromosome and position
     *
     * @param value Site Filter
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin siteFilter(SITE_RANGE_FILTER_TYPES value) {
        mySiteFilter = new PluginParameter<>(mySiteFilter, value);
        return this;
    }

    /**
     * Start Site
     *
     * @return Start Site
     */
    public Integer startSite() {
        return myStartSite.value();
    }

    /**
     * Set Start Site. Start Site
     *
     * @param value Start Site
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin startSite(Integer value) {
        myStartSite = new PluginParameter<>(myStartSite, value);
        return this;
    }

    /**
     * End Site
     *
     * @return End Site
     */
    public Integer endSite() {
        return myEndSite.value();
    }

    /**
     * Set End Site. End Site
     *
     * @param value End Site
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin endSite(Integer value) {
        myEndSite = new PluginParameter<>(myEndSite, value);
        return this;
    }

    /**
     * Start Chr
     *
     * @return Start Chr
     */
    public Chromosome startChr() {
        return myStartChr.value();
    }

    /**
     * Set Start Chr. Start Chr
     *
     * @param value Start Chr
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin startChr(Chromosome value) {
        myStartChr = new PluginParameter<>(myStartChr, value);
        return this;
    }

    /**
     * Start Pos
     *
     * @return Start Pos
     */
    public Integer startPos() {
        return myStartPos.value();
    }

    /**
     * Set Start Pos. Start Pos
     *
     * @param value Start Pos
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin startPos(Integer value) {
        myStartPos = new PluginParameter<>(myStartPos, value);
        return this;
    }

    /**
     * End Chr
     *
     * @return End Chr
     */
    public Chromosome endChr() {
        return myEndChr.value();
    }

    /**
     * Set End Chr. End Chr
     *
     * @param value End Chr
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin endChr(Chromosome value) {
        myEndChr = new PluginParameter<>(myEndChr, value);
        return this;
    }

    /**
     * End Pos
     *
     * @return End Pos
     */
    public Integer endPos() {
        return myEndPos.value();
    }

    /**
     * Set End Pos. End Pos
     *
     * @param value End Pos
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin endPos(Integer value) {
        myEndPos = new PluginParameter<>(myEndPos, value);
        return this;
    }

    /**
     * Include Site Names
     *
     * @return Include Site Names
     */
    public Boolean includeSiteNames() {
        return myIncludeSites.value();
    }

    /**
     * Set Include Site Names. Include Site Names
     *
     * @param value Include Site Names
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin includeSiteNames(Boolean value) {
        myIncludeSites = new PluginParameter<>(myIncludeSites, value);
        return this;
    }

    /**
     * Position List
     *
     * @return Position List
     */
    public PositionList positionList() {
        return myPositionList.value();
    }

    /**
     * Set Position List.
     *
     * @param value Position List
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin positionList(PositionList value) {
        myPositionList = new PluginParameter<>(myPositionList, value);
        return this;
    }

    /**
     * Site Names List
     *
     * @return Site Names List
     */
    public String siteNamesList() {
        return mySiteNamesList.value();
    }

    /**
     * Set Site Names List. Site Names List
     *
     * @param value Site Names List
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin siteNamesList(String value) {
        mySiteNamesList = new PluginParameter<>(mySiteNamesList, value);
        return this;
    }

    /**
     * Filter based on BED file.
     *
     * @return Bed File
     */
    public String bedFile() {
        return myBedFile.value();
    }

    /**
     * Set Bed File. Filter based on BED file.
     *
     * @param value Bed File
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin bedFile(String value) {
        myBedFile = new PluginParameter<>(myBedFile, value);
        return this;
    }

    /**
     * Filter based on list of chromsome / position in file.
     *
     * @return Chr Pos File
     */
    public String chrPosFile() {
        return myChrPosFile.value();
    }

    /**
     * Set Chr Pos File. Filter based on list of chromsome / position in file.
     *
     * @param value Chr Pos File
     *
     * @return this plugin
     */
    public FilterSiteBuilderPlugin chrPosFile(String value) {
        myChrPosFile = new PluginParameter<>(myChrPosFile, value);
        return this;
    }

}
