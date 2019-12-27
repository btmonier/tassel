package net.maizegenetics.analysis.imputation;

import com.google.common.collect.Range;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.LoggingUtils;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

public class ProjectNamPlugin extends AbstractPlugin {
    PluginParameter<String> progenyInputDirectory = new PluginParameter.Builder<>("inputDir", "/Users/peterbradbury/temp/baoxing/USNAM_buildV4_imputed_var_founders/", String.class)
            .required(false)
            .guiName("Progeny Input Directory")
            .description("The parent calls for the NAM progeny")
            .inDir()
            .build();

    PluginParameter<String> progenyOutput = new PluginParameter.Builder<>("output", "/Users/peterbradbury/temp/baoxing/namProjectedCNS.hmp.txt", String.class)
            .required(false)
            .guiName("Progeny Output File")
            .description("The parent calls for the NAM progeny")
            .outFile()
            .build();

    GenotypeTable parentGT;
    private final String[][] namFilenames = new String[25][10];

    @Override
    public DataSet performFunction(DataSet input) {
        createNamFilenames();

        //takes the parent genotypes as input or imports
        if (input == null) parentGT = importAndFilterParents();
        else {
            List<Datum> inputGenotypes = input.getDataOfType(GenotypeTable.class);
            if (inputGenotypes.size() == 1) {
                parentGT = (GenotypeTable) inputGenotypes.get(0).getData();
            } else {
                throw new IllegalArgumentException("input must be null or contain one GenotypeTable");
            }
        }

        //which taxon is B73?
        Taxon B73gt = parentGT.taxa().stream().filter(t -> t.getName().equals("B73")).findFirst().get();
        int B73GtIndex = parentGT.taxa().indexOf(B73gt);

        //get the PositionList for positions to be imputed
        PositionList sitesToImpute = parentGT.positions();
        int numberOfSites = sitesToImpute.numberOfSites();
        int[] chrOffsets = sitesToImpute.chromosomesOffsets();

        //Make a taxa incremental GenotypeTable builder
        GenotypeTableBuilder genoBuilder = GenotypeTableBuilder.getTaxaIncremental(sitesToImpute);

        //for each taxon, project parent sites, then add a byte[] of genotypes
        for (int fam = 0; fam < 25; fam++) {
            Map<Taxon, byte[]> taxonGenotypeMap = new HashMap();

            //for each chromosome
            for (int chrndx = 0; chrndx < 10; chrndx++) {
                Chromosome myChr = Chromosome.instance(chrndx + 1);

                //   open the sorted file as a genotype table
                String filename = namFilenames[fam][chrndx];
                System.out.println("processing " + filename);
                FileLoadPlugin fileLoader = new FileLoadPlugin(null, false)
                        .sortPositions(true)
                        .keepDepth(false)
                        .fileType(FileLoadPlugin.TasselFileType.Hapmap);
                fileLoader.setOpenFiles(filename);
                GenotypeTable genoAll = (GenotypeTable) fileLoader.performFunction(null).getData(0).getData();

                //filter the nam progeny genotypes to this chromosome
                //this step is necessary to filter out sites that projected to a different chromosome
                DataSet ds = new FilterSiteBuilderPlugin(null, false)
                        .siteFilter(FilterSite.SITE_RANGE_FILTER_TYPES.POSITIONS)
                        .startChr(myChr).endChr(myChr)
                        .startPos(1).endPos(Integer.MAX_VALUE)
                        .performFunction(DataSet.getDataSet(genoAll));

                GenotypeTable geno = (GenotypeTable) ds.getDataOfType(GenotypeTable.class).get(0).getData();

                //find B73 and the other parent from the nam parent calls
                Taxon B73 = geno.taxa().stream().filter(t -> t.getName().startsWith("B73")).findFirst().orElseGet(() -> new Taxon("B73"));
                Taxon other = geno.taxa().stream().filter(t -> !t.getName().startsWith("B73") && !t.getName().startsWith("Z0")).findFirst().get();
                Optional<Taxon> optOtherGT = parentGT.taxa().stream().filter(t -> other.getName().startsWith(t.getName())).findFirst();
                if (!optOtherGT.isPresent()) continue;
                int otherGtIndex = parentGT.taxa().indexOf(optOtherGT.get());

                Map<Byte,Long> countMap = geno.streamGenotype(geno.taxa().indexOf(B73))
                        .filter(b -> b != GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
                        .collect(Collectors.groupingBy(b -> b, Collectors.counting()));
                Byte B73Byte = 0;
                long genotypeCount = 0;
                for (Map.Entry<Byte, Long> ent : countMap.entrySet()) {
                    if (countMap.get(ent.getKey()) > genotypeCount) {
                        genotypeCount = countMap.get(ent.getKey());
                        B73Byte = ent.getKey();
                    }
                }

                if (GenotypeTableUtils.isHeterozygous(B73Byte)) {
                    throw new IllegalArgumentException("B73 genotype is het for " + filename);
                }
                String B73AlleleString = NucleotideAlignmentConstants.getNucleotideIUPAC(B73Byte);
                byte B73Genotype, otherGenotype, het;

                if (B73AlleleString.equals("A")) {
                    B73Genotype = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
                    otherGenotype = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
                    het = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
                } else if (B73AlleleString.equals("C")) {
                    B73Genotype = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
                    otherGenotype = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
                    het = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
                } else {
                    throw new IllegalArgumentException("Illegal value for B73 allele : " + B73AlleleString + " in " + filename);
                }

                //   for each taxon:
                for (Taxon taxon : geno.taxa()) {
                    //      if this is the first chromosome, create a taxon -> genotype byte array map
                    if (chrndx == 0) taxonGenotypeMap.put(taxon, new byte[numberOfSites]);
                    byte[] myProjectedGenotype = taxonGenotypeMap.get(taxon);

                    //      create a range map with parent values, 0 = B73, 2 = non-B73, 1 = het
                    //TODO add ranges for start and end of chromosome unknown
                    TreeRangeMap<Integer, Integer> taxonRangeMap = TreeRangeMap.create();
                    int taxonIndex = geno.taxa().indexOf(taxon);
                    int startPos = -1;
                    int endPos = -1;
                    byte unknownGenotype = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                    byte startGenotype = unknownGenotype;
                    for (int s = 0; s < geno.numberOfSites(); s++) {
                        byte siteGenotype = geno.genotype(taxonIndex, s);
                        if (siteGenotype != unknownGenotype) {
                            if (startGenotype == unknownGenotype) {
                                startPos = geno.chromosomalPosition(s);
                                endPos = startPos;
                                startGenotype = siteGenotype;
                            } else if (siteGenotype == startGenotype) { //equal and not unknown
                                endPos = geno.chromosomalPosition(s);
                            } else {  //not equal and not unknown
                                if (startGenotype == B73Genotype) taxonRangeMap.put(Range.closed(startPos, endPos), 0);
                                else if (startGenotype == otherGenotype) taxonRangeMap.put(Range.closed(startPos, endPos), 2);
                                else if (startGenotype == het) taxonRangeMap.put(Range.closed(startPos, endPos), 1);
                                startPos = geno.chromosomalPosition(s);
                                endPos = startPos;
                                startGenotype = siteGenotype;
                            }
                        }
                    }

                    //add the final interval
                    if (startGenotype == B73Genotype) taxonRangeMap.put(Range.closed(startPos, endPos), 0);
                    else if (startGenotype == otherGenotype) taxonRangeMap.put(Range.closed(startPos, endPos), 2);
                    else if (startGenotype == het) taxonRangeMap.put(Range.closed(startPos, endPos), 1);

                    //      for each position in chromosome look up the value in the range map
                    //      set value in byte array
                    int startSite = chrOffsets[chrndx];
                    int endSite = numberOfSites;
                    if (chrndx < 9) endSite = chrOffsets[chrndx + 1];
                    for (int s = startSite; s < endSite; s++) {
                        Integer genoval = taxonRangeMap.get(sitesToImpute.get(s).getPosition());
                        if (genoval == null) {
                            myProjectedGenotype[s] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                        } else if (genoval == 0) {
                            myProjectedGenotype[s] = parentGT.genotype(B73GtIndex, s);
                        } else if (genoval == 2) {
                            myProjectedGenotype[s] = parentGT.genotype(otherGtIndex, s);
                        } else if (genoval == 1) {
                            myProjectedGenotype[s] = heterozygote(parentGT.genotype(B73GtIndex, s), parentGT.genotype(otherGtIndex, s));
                        } else {
                            myProjectedGenotype[s] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                        }
                    }
                }

            }
            //add the taxa to the genotype builder from the taxonGenotypeMap
            taxonGenotypeMap.entrySet().stream().filter(ent -> !ent.getKey().getName().startsWith("B73")).forEach(ent -> {
                genoBuilder.addTaxon(ent.getKey(), ent.getValue());
            });
        }

        //export GenotypeTable to a file
        LoggingUtils.setupDebugLogging();
        ExportUtils.writeToHapmap(genoBuilder.build(), progenyOutput.value());
        return null;
    }

    private byte heterozygote(byte g1, byte g2) {
        if (g1 == g2) return g1;
        if (g1 == GenotypeTable.UNKNOWN_DIPLOID_ALLELE || g2 == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        if (GenotypeTableUtils.isHeterozygous(g1)) return g1;
        if (GenotypeTableUtils.isHeterozygous(g1)) return g2;
        else return GenotypeTableUtils.getUnphasedDiploidValue(GenotypeTableUtils.getDiploidValues(g1)[0], GenotypeTableUtils.getDiploidValues(g2)[0]);
    }

    private void createNamFilenames() {
        String inputDir = progenyInputDirectory.value();
        if (!(inputDir.endsWith("/"))) inputDir +="/";

        for (int fam = 0; fam < 25; fam++) {
            String famstr;
            if (fam < 16) {
                famstr = ("0" + (fam + 1));
            } else {
                famstr = ("0" + (fam + 2));
            }
            if (fam < 9) famstr = "0" + (fam + 1);
            else if (fam < 16) famstr = Integer.toString(fam + 1);
            else famstr = Integer.toString(fam + 2);

            for (int chr = 0; chr < 10; chr++) {
                //familySuffix[i] = String.format(".family.Z0%sparents.hmp.txt", famstr);
                namFilenames[fam][chr] = String.format( inputDir + "USNAM_buildV4_imputed_var_founders_chr%d.family.Z0%sparents.hmp.txt", chr + 1, famstr);
            }
        }
    }

    @Override
    public String pluginDescription() {
        return "Project parent genotypes on to NAM progeny";
    }

    private void reportParentAlleleFrequencies() {
        String outfile = "/Users/peterbradbury/temp/parent_allele_counts.txt";
        try(PrintWriter pw = new PrintWriter(outfile)) {
            for (int fam = 0; fam < 25; fam++) {
                for (int chr = 0; chr < 10; chr++) {
                    Chromosome myChr = Chromosome.instance(chr + 1);
                    String filename = namFilenames[fam][chr];
                    FileLoadPlugin fileLoader = new FileLoadPlugin(null, false)
                            .sortPositions(true)
                            .keepDepth(false)
                            .fileType(FileLoadPlugin.TasselFileType.Hapmap);
                    fileLoader.setOpenFiles(filename);


                    GenotypeTable allGenos = (GenotypeTable) fileLoader.performFunction(null).getData(0).getData();
                    DataSet ds = new FilterSiteBuilderPlugin(null, false)
                            .siteFilter(FilterSite.SITE_RANGE_FILTER_TYPES.POSITIONS).startChr(myChr).endChr(myChr)
                            .startPos(1).endPos(Integer.MAX_VALUE)
                            .processData(DataSet.getDataSet(allGenos));
                    GenotypeTable geno = (GenotypeTable) ds.getDataOfType(GenotypeTable.class).get(0).getData();

                    geno.taxa().stream().filter(t -> !t.getName().startsWith("Z0")).forEach(t -> {
                        Map<String,Long> countMap = geno.streamGenotype(geno.taxa().indexOf(t))
                                .collect(Collectors.groupingBy(b -> NucleotideAlignmentConstants.getNucleotideIUPAC(b), Collectors.counting()));
                        pw.printf("%s: %s, %s%n", filename, t.getName(), countMap.toString());
                    });
                }
                pw.println("----------------------");
            }

        } catch(IOException e) {
            throw new IllegalArgumentException("file write error");
        }
        System.out.println("Finished.");
    }

    private GenotypeTable importAndFilterParents() {
        String parentFilename = "/Users/peterbradbury/temp/baoxing/NAM.parents.hapmap";
        FileLoadPlugin loader = new FileLoadPlugin(null, false)
                .fileType(FileLoadPlugin.TasselFileType.Hapmap)
                .sortPositions(true)
                .keepDepth(false);
        loader.setOpenFiles(parentFilename);
        GenotypeTable importedGeno = (GenotypeTable) loader.performFunction(null).getData(0).getData();

        return (GenotypeTable) new FilterSiteBuilderPlugin(null, false)
                .siteMinCount(10)
                .siteMinAlleleFreq(0.01)
                .performFunction(DataSet.getDataSet(importedGeno))
                .getDataOfType(GenotypeTable.class)
                .get(0).getData();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    public static void main(String[] args) {
//        LoggingUtils.setupDebugLogging();
        ProjectNamPlugin pnp = new ProjectNamPlugin();

//        pnp.createNamFilenames();
//        pnp.reportParentAlleleFrequencies();

        projectData();
    }

    public static void projectData() {
        ProjectNamPlugin pnp = new ProjectNamPlugin();
        pnp.performFunction(null);
    }

}
