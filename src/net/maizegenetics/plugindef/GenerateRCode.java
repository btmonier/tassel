/*
 *  GenerateRCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.awt.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.Date;

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
            if ((current.defaultValue() instanceof Number)) {
                sb.append("            " + methodName + "=" + current.defaultValue() + ",\n");
            } else if ((current.defaultValue() instanceof Enum)) {
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
            if ((current.defaultValue() instanceof Number)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Enum)) {
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
    public static int[] genotypeTableToDosageIntArray(GenotypeTable genotype) {

        int[] result = new int[genotype.numberOfTaxa() * genotype.numberOfSites()];

        int index = 0;
        for (int site = 0; site < genotype.numberOfSites(); site++) {
            byte[] siteGenotypes = genotype.genotypeAllTaxa(site);
            int[][] alleleCounts = AlleleFreqCache.allelesSortedByFrequencyNucleotide(siteGenotypes);
            byte majorAllele = AlleleFreqCache.majorAllele(alleleCounts);
            // value assigned to site / taxon is the number of alleles
            // that doesn't not match the major allele.
            for (int taxon = 0; taxon < genotype.numberOfTaxa(); taxon++) {
                int value = 0;
                byte[] alleles = GenotypeTableUtils.getDiploidValues(siteGenotypes[taxon]);
                if (alleles[0] == GenotypeTable.UNKNOWN_ALLELE || alleles[1] == GenotypeTable.UNKNOWN_ALLELE) {
                    value = Integer.MIN_VALUE;
                } else {
                    if (alleles[0] != majorAllele) value++;
                    if (alleles[1] != majorAllele) value++;
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
        return genotype.taxa().stream()
                .map(Taxon::getName)
                .toArray(String[]::new);
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param genotype
     *
     * @return int[] in column order with NA set to R approach
     */
    public static PositionVectors genotypeTableToPositionListOfArrays(GenotypeTable genotype) {

        String[] chromosomes = new String[genotype.numberOfSites()];
        int[] startPos = new int[genotype.numberOfSites()];
        int[] strand = new int[genotype.numberOfSites()];
        String[] refAllele = new String[genotype.numberOfSites()];
        String[] altAllele = new String[genotype.numberOfSites()];

        for (int site = 0; site < genotype.numberOfSites(); site++) {
            Position p = genotype.positions().get(site);
            chromosomes[site] = p.getChromosome().getName();
            startPos[site] = p.getPosition();
            strand[site] = p.getStrand();
            String[] variants = p.getKnownVariants();
            refAllele[site] = (variants.length > 0) ? variants[0] : "";
            altAllele[site] = (variants.length > 1) ? variants[1] : "";
        }

        // return startPos;
        return new PositionVectors(chromosomes, startPos, strand, refAllele, altAllele);
//        //This can be replace with a custom class in the future or List.of in Java 9
//        java.util.List list=new ArrayList<>();
//        list.add(chromosomes);
//        list.add(startPos);
//        list.add(strand);
//        list.add(refAllele);
//        list.add(altAllele);
//        return list;
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

}
