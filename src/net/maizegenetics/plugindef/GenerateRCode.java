/*
 *  GenerateRCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.analysis.association.FastMultithreadedAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.CorePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.awt.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static net.maizegenetics.dna.map.Position.STRAND_MINUS;
import static net.maizegenetics.dna.map.Position.STRAND_PLUS;

/**
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenerateRCode {

    private static final Logger myLogger = Logger.getLogger(GenerateRCode.class);

    private GenerateRCode() {

    }

    public static void main(String[] args) {
        printHeader();
        generate(FilterSiteBuilderPlugin.class, "genotypeTable", "genotypeTable");
        generate(FilterTaxaBuilderPlugin.class, "genotypeTable", "genotypeTable");
        generate(KinshipPlugin.class, "genotypeTable", "distanceMatrix");
        generate(FixedEffectLMPlugin.class, "phenotypeGenotypeTable", "tableReport");

    }

    private static void printHeader() {

        System.out.println("#!/usr/bin/env Rscript");

        System.out.println("\n#--------------------------------------------------------------------");
        System.out.println("# Script Name:   TasselPluginWrappers.R");
        System.out.println("# Description:   Generated R interface to TASSEL 5");
        System.out.println("# Author:        Brandon Monier, Ed Buckler, Terry Casstevens");
        System.out.print("# Created:       ");
        System.out.println(new Date());
        System.out.println("#--------------------------------------------------------------------");

        System.out.println("# Preamble\n");

        System.out.println("\n## Load packages");
        System.out.println("if (!requireNamespace(\"BiocManager\")) {");
        System.out.println("    install.packages(\"BiocManager\")");
        System.out.println("}");

        System.out.println("\npackages <- c(");
        System.out.println("\"rJava\"");
        System.out.println(")");
        System.out.println("BiocManager::install(packages)");
        System.out.println("library(rJava)");

        System.out.println("\n## Init JVM");
        System.out.println("rJava::.jinit()");

        System.out.println("\n## Add TASSEL 5 class path");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/lib\")");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/sTASSEL.jar\")\n");

        System.out.println("source(\"R/AllClasses.R\")\n");

    }

    public static void generate(Class currentMatch, String inputObject, String outputObject) {
        try {
            Constructor constructor = currentMatch.getConstructor(Frame.class);
            generate((AbstractPlugin) constructor.newInstance((Frame) null), inputObject, outputObject);
        } catch (Exception ex) {
            try {
                Constructor constructor = currentMatch.getConstructor(Frame.class, boolean.class);
                generate((AbstractPlugin) constructor.newInstance(null, false), inputObject, outputObject);
            } catch (NoSuchMethodException nsme) {
                myLogger.warn("Self-describing Plugins should implement this constructor: " + currentMatch.getClass().getName());
                myLogger.warn("public Plugin(Frame parentFrame, boolean isInteractive) {");
                myLogger.warn("   super(parentFrame, isInteractive);");
                myLogger.warn("}");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void generate(AbstractPlugin plugin, String inputObject, String outputObject) {

        String clazz = Utils.getBasename(plugin.getClass().getName());
        //StringBuilder sb = new StringBuilder("rTASSEL::" + stringToCamelCase(clazz));
        StringBuilder sb = new StringBuilder(stringToCamelCase(clazz));
        sb.append(" <- function(");
        sb.append(inputObject + ",\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            String methodName = current.cmdLineName();
            if ((current.defaultValue() instanceof Number) || (current.defaultValue() instanceof Boolean) || (current.defaultValue() instanceof String)) {
                sb.append("            " + methodName + "=" + current.defaultValue() + ",\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                sb.append("            " + methodName + "=\"" + current.defaultValue() + "\",\n");
            } else if ((current.defaultValue() instanceof Boolean)) {
                sb.append("            " + methodName + "=\"" + current.defaultValue() + "\",\n");
            }

        }
        sb.deleteCharAt(sb.lastIndexOf(",")); //remove last comma

        sb.append(") {\n");

        if (inputObject.equals("genotypeTable")) {
            sb.append("    if(is(genotypeTable, \"GenotypeTable\") == TRUE) {\n");
            sb.append("        genotypeTable <- genotypeTable@jtsGenotypeTable\n");
            sb.append("    }\n");
        }

        sb.append("    plugin <- new(J(\"" + plugin.getClass().getCanonicalName() + "\"), .jnull(), FALSE)\n");
        //sb.append("    plugin <- rJava::.jnew(\"" + plugin.getClass().getCanonicalName() + "\")\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            String methodName = current.cmdLineName();
            if ((current.defaultValue() instanceof Number) || (current.defaultValue() instanceof Boolean) || (current.defaultValue() instanceof String)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Boolean)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            }
        }

        if (outputObject.equals("genotypeTable")) {
            sb.append("    filteredGT <- plugin$runPlugin(" + inputObject + ")\n");
            sb.append("    new(\n");
            sb.append("        Class = \"GenotypeTable\",\n");
            sb.append("        name = paste0(\"Filtered:\"),\n");
            sb.append("        jtsGenotypeTable = filteredGT\n");
            sb.append("    )\n");
        } else {
            sb.append("    plugin$runPlugin(" + inputObject + ")\n");
        }

        sb.append("}\n");
        System.out.println(sb.toString());

    }

    private static final int DEFAULT_DESCRIPTION_LINE_LENGTH = 50;

    private static String createDescription(String description) {
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("     * ");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else if ((count > DEFAULT_DESCRIPTION_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        return builder.toString();
    }

    private static String stringToCamelCase(String str) {
        StringBuilder builder = new StringBuilder();
        builder.append(Character.toLowerCase(str.charAt(0)));
        boolean makeUpper = false;
        for (int i = 1; i < str.length(); i++) {
            char current = str.charAt(i);
            if (current == ' ') {
                makeUpper = true;
            } else if (makeUpper) {
                builder.append(Character.toUpperCase(current));
                makeUpper = false;
            } else {
                builder.append(current);
            }
        }
        return builder.toString();
    }

    private static String removeMyFromString(String str) {
        String lower = str.toLowerCase();
        if (lower.startsWith("my")) {
            str = str.substring(2);
        }
        return CaseFormat.UPPER_CAMEL.to(CaseFormat.LOWER_CAMEL, str);
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param genotype
     *
     * @return int[] in column order with NA set to R approach
     */
    public static int[] genotypeTableToDosageIntArray(GenotypeTable genotype, boolean useRef) {

        int[] result = new int[genotype.numberOfTaxa() * genotype.numberOfSites()];

        int index = 0;
        PositionList posList = genotype.positions();
        for (int site = 0; site < genotype.numberOfSites(); site++) {
            byte[] siteGenotypes = genotype.genotypeAllTaxa(site);
            int[][] alleleCounts = AlleleFreqCache.allelesSortedByFrequencyNucleotide(siteGenotypes);
            byte majorAllele = AlleleFreqCache.majorAllele(alleleCounts);

            // Use reference if available.  If not, use major allele
            Position posAtSite = posList.get(site);
            byte refAllele = Byte.MAX_VALUE;
            if (useRef) {
                myLogger.info("genotypeTableToDosageIntArray: using refAlle at sites");
                refAllele = posAtSite.getAllele(WHICH_ALLELE.Reference);
            } else {
                myLogger.info("genotypeTableToDosageIntArray: using majorAllele for ref at sites");
                refAllele = AlleleFreqCache.majorAllele(alleleCounts);
            }

            // value assigned to site / taxon is the number of alleles
            // that doesn't not match the major allele.
            for (int taxon = 0; taxon < genotype.numberOfTaxa(); taxon++) {
                int value = 0;
                byte[] alleles = GenotypeTableUtils.getDiploidValues(siteGenotypes[taxon]);
                if (alleles[0] == GenotypeTable.UNKNOWN_ALLELE || alleles[1] == GenotypeTable.UNKNOWN_ALLELE) {
                    value = Integer.MIN_VALUE;
                } else {
                    if (alleles[0] != refAllele) value++;
                    if (alleles[1] != refAllele) value++;
                }
                result[index++] = value;
            }
        }
        return result;
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param genotype
     *
     * @return int[] in column order with NA set to R approach
     */
    public static String[] genotypeTableToSampleNameArray(GenotypeTable genotype) {
        return genotypeTableToSampleNameArray(genotype.taxa());
    }

    public static String[] genotypeTableToSampleNameArray(TaxaList taxa) {
        return taxa.stream()
                .map(Taxon::getName)
                .toArray(String[]::new);
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param positions
     *
     * @return int[] in column order with NA set to R approach
     */
    public static PositionVectors genotypeTableToPositionListOfArrays(PositionList positions) {

        String[] chromosomes = new String[positions.numberOfSites()];
        int[] startPos = new int[positions.numberOfSites()];
        int[] strand = new int[positions.numberOfSites()];
        String[] refAllele = new String[positions.numberOfSites()];
        String[] altAllele = new String[positions.numberOfSites()];

        for (int site = 0; site < positions.numberOfSites(); site++) {
            Position p = positions.get(site);
            chromosomes[site] = p.getChromosome().getName();
            startPos[site] = p.getPosition();
            strand[site] = (p.getStrand() == STRAND_PLUS) ? (1) : ((p.getStrand() == STRAND_MINUS) ? (-1) : (Integer.MIN_VALUE));
            String[] variants = p.getKnownVariants();
            refAllele[site] = (variants.length > 0) ? variants[0] : "";
            altAllele[site] = (variants.length > 1) ? variants[1] : "";
        }

        return new PositionVectors(chromosomes, startPos, strand, refAllele, altAllele);
    }

    public static PositionVectors genotypeTableToPositionListOfArrays(GenotypeTable genotype) {
        return genotypeTableToPositionListOfArrays(genotype.positions());
    }

    public static class PositionVectors {
        public String[] chromosomes;
        public int[] startPos;
        public int[] strand;
        public String[] refAllele;
        public String[] altAllele;

        public PositionVectors(String[] chromosomes, int[] startPos, int[] strand, String[] refAllele, String[] altAllele) {
            this.chromosomes = chromosomes;
            this.startPos = startPos;
            this.strand = strand;
            this.refAllele = refAllele;
            this.altAllele = altAllele;
        }
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param positions
     *
     * @return int[] in column order with NA set to R approach
     */
    public static TableReportVectors tableReportToVectors(TableReport tableReport) {
        if (tableReport.getRowCount() > Integer.MAX_VALUE)
            throw new IndexOutOfBoundsException("R cannot handle more than " + Integer.MAX_VALUE + " rows");
        int rows = (int) tableReport.getRowCount();
        String[] columnNames = Stream.of(tableReport.getTableColumnNames()).map(Object::toString).toArray(String[]::new);
        List dataVector = new ArrayList();
        String[] columnType = Stream.of(tableReport.getRow(0)).map(t -> t.getClass().toString()).toArray(String[]::new);
        for (int column = 0; column < tableReport.getColumnCount(); column++) {
            Object o = tableReport.getValueAt(0, column);
            if (o instanceof Float || o instanceof Double) {
                double[] result = new double[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    result[row] = ((Number) tableReport.getValueAt(row, column)).doubleValue();
                }
                dataVector.add(result);
            } else if (o instanceof Byte || o instanceof Short || o instanceof Integer || o instanceof Long) {
                int[] result = new int[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    result[row] = ((Number) tableReport.getValueAt(row, column)).intValue();
                }
                dataVector.add(result);
            } else {
                String[] result = new String[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    result[row] = tableReport.getValueAt(row, column).toString();
                }
                dataVector.add(result);
            }
        }
        String[] annotation = new String[tableReport.getColumnCount()];
        return new TableReportVectors(columnNames, columnType, annotation, dataVector);
    }

    public static class TableReportVectors {
        public String[] columnNames;
        public String[] columnType;
        public String[] annotation;
        public List dataVector = new ArrayList();

        public TableReportVectors(String[] columnNames, String[] columnType, String[] annotation, List dataVector) {
            this.columnNames = columnNames;
            this.columnType = columnType;
            this.annotation = annotation;
            this.dataVector = dataVector;
        }
    }

    public static Phenotype createPhenotypeFromRDataFrameElements(String[] taxaArray, String[] colNames, String[] attributeType, List dataVectors) {
        List<Taxon> taxaList = Stream.of(taxaArray).map(Taxon::new).collect(Collectors.toList());

//        List<Taxon> taxaList = new ArrayList<Taxon>();
//        taxaList.add(new Taxon("t1"));
//        taxaList.add(new Taxon("t2"));
//        taxaList.add(new Taxon("t3"));
//        taxaList.add(new Taxon("t1"));
//        taxaList.add(new Taxon("t2"));
//        taxaList.add(new Taxon("t3"));


        List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
        List<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
        TaxaAttribute ta = new TaxaAttribute(taxaList);
        attributes.add(ta);
        types.add(Phenotype.ATTRIBUTE_TYPE.taxa);

        //TODO test for size of colNames
        if (colNames.length != attributeType.length || colNames.length != dataVectors.size())
            throw new IllegalArgumentException("ColNames, attributeType, and dataVectors need to be same size");
        for (int i = 0; i < colNames.length; i++) {
            Object o = dataVectors.get(i);
            if (o instanceof double[]) {
                attributes.add(new NumericAttribute(colNames[i], (double[]) o));
                System.out.println(attributeType[i]);
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                types.add(attribute_type);
            } else if (o instanceof int[]) {
                //INT min is NA
                int[] initialTraits = (int[]) o;
                float[] traitsValueWithNaN = new float[initialTraits.length];
                for (int j = 0; j < initialTraits.length; j++)
                    traitsValueWithNaN[j] = (initialTraits[j] == Integer.MIN_VALUE) ? Float.NaN : initialTraits[j];
                attributes.add(new NumericAttribute(colNames[i], traitsValueWithNaN));
                System.out.println(attributeType[i]);
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                types.add(attribute_type);
            } else if (o instanceof String[]) {
                attributes.add(new CategoricalAttribute(colNames[i], (String[]) o));
                System.out.println(attributeType[i]);
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                types.add(attribute_type);
            } else {
                throw new IllegalArgumentException("Unsupported type for phenotype table");
            }

        }
//
//        CategoricalAttribute factor = new CategoricalAttribute("factor", new String[]{"A","A","A","B","B","B"});
//        NumericAttribute trait = new NumericAttribute("trait1", new float[]{1,2,3,4,5,6});
//
//
//        attributes.add(factor);
//        attributes.add(trait);
//
//
//        types.add(Phenotype.ATTRIBUTE_TYPE.factor);
//        types.add(Phenotype.ATTRIBUTE_TYPE.data);

        return new CorePhenotype(attributes, types, "pheno1");
    }

    //    traitName traitType       traitAttribute isReponse isPredictor   newType effect
    // 1      Taxa      taxa        TaxaAttribute     FALSE       FALSE      <NA>   <NA>
    // 2  location    factor CategoricalAttribute     FALSE        TRUE    factor  fixed
    // 3     EarHT      data     NumericAttribute      TRUE       FALSE      data   <NA>
    // 4     dpoll      data     NumericAttribute     FALSE       FALSE      <NA>   <NA>
    // 5    EarDia      data     NumericAttribute      TRUE       FALSE      data   <NA>
    // 6        Q1 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 7        Q2 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 8        Q3 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 9         G  genotype             Genotype     FALSE        TRUE  genotype  fixed

    public static TableReport association(DistanceMatrix kinship, GenotypeTable genotype, Phenotype phenotype, GenotypePhenotype genoPheno) {

        if (genotype == null && phenotype == null) {
            myLogger.warn("association: genotype and phenotype are null.  Nothing calculated");
        }

        if (kinship == null) {

            myLogger.info("association: running GLM");

            FixedEffectLMPlugin plugin = new FixedEffectLMPlugin(null, false);

            DataSet input = null;

            if (genotype == null) {
                plugin.phenotypeOnly(true);
                input = DataSet.getDataSet(phenotype);
            } else {
                plugin.phenotypeOnly(false);
                input = DataSet.getDataSet(genoPheno);
            }

            DataSet output = plugin.performFunction(input);

            if (output != null) {
                return (TableReport) output.getDataOfType(TableReport.class).get(0).getData();
            }

            return null;

        } else {

            myLogger.info("association: running MLM");

            MLMPlugin plugin = new MLMPlugin(null, false);

            Datum genoDatum = new Datum("GenotypePhenotype", genoPheno, null);
            Datum kinshipDatum = new Datum("Kinship", kinship, null);
            DataSet input = new DataSet(new Datum[]{genoDatum, kinshipDatum}, null);

            DataSet output = plugin.performFunction(input);

            if (output != null) {
                return (TableReport) output.getDataOfType(TableReport.class).get(0).getData();
            }

            return null;

        }

    }

    public static TableReport fastAssociation(GenotypePhenotype genoPheno) {

        if (genoPheno == null) {
            myLogger.warn("fastAssociation: GenotypePhenotype is null.  Nothing calculated");
        }

        FastMultithreadedAssociationPlugin plugin = new FastMultithreadedAssociationPlugin(null, false);

        DataSet input = DataSet.getDataSet(genoPheno);

        DataSet output = plugin.performFunction(input);

        if (output != null) {
            return (TableReport) output.getDataOfType(TableReport.class).get(0).getData();
        }

        return null;

    }


}
