/*
 *  GenotypeTableBuilder
 */
package net.maizegenetics.dna.snp;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.snp.depth.AlleleDepthBuilder;
import net.maizegenetics.dna.snp.depth.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeMergeRule;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Builder for GenotypeTables.  New genotypeTables are built from a minimum of TaxaList, PositionList, and
 * GenotypeCallTable.  Depth and Scores are optional features of GenotypeTables.
 * <p></p>
 If you know the taxa,position, and genotypes are known from the beginning use:
 GenotypeTable a=GenotypeTableBuilder.getInstance(genotype, positionList, taxaList);

 In many situations only GenotypeTables are built incrementally, either by Taxa or Site.
 <p></p>
 For taxa building:
 <pre>
 *{@code
 *    GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(gbs.positions(),outFile);
 *    for (int i=0; i<hm2.numberOfTaxa(); i++) {
 *        Taxon taxon=hm2.taxa().get(i);
 *        byte[] geno=hm2.genotypeAllSites(i);
 *        gtb.addTaxon(taxon,geno);
 *        }
 *    GenotypeTable gt=gtb.build();
 *}
 </pre>
 <p></p>
 In many cases, genotype want to add taxa to an existing genotypeTable.  Direct addition is not possible, as GenotypeTables
 are immutable, but the GenotypeTableBuilder.getTaxaIncremental provides a strategy for creating and merging taxa together.
 Key to the process is that GenotypeMergeRule defines how the taxa with identical names will be merged.<br></br>
 Merging is possible with HDF5 files, but only if the closeUnfinished() method was used with the previous building.
 <pre>{@code
GenotypeTable existingGenotypeTable1, existingGenotypeTable2;
GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(existingGenotypeTable1,
    new BasicGenotypeMergeRule(0.01));
for (int i=0; i<existingGenotypeTable2.numberOfTaxa(); i++) {
    gtb.addTaxon(existingGenotypeTable2.taxa().get(i), existingGenotypeTable2.genotypeAllSites(i)
        existingGenotypeTable2.depth().depthAllSitesByte(i));
}
 }</pre>

 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenotypeTableBuilder {

    private GenotypeCallTable genotype=null;

    //Fields for incremental taxa
    private PositionList positionList=null;
    private TaxaListBuilder taxaListBuilder=null;
    private ArrayList<byte[]> incGeno=null;
    private ArrayList<byte[][]> incDepth=null;
    private HashMap<Taxon,Integer> incTaxonIndex=null;
    private boolean sortAlphabetically=false;

    //Fields for incremental sites
    private TaxaList taxaList=null;
    private PositionListBuilder posListBuilder=null;
    private boolean isTaxaMerge=false; //if in taxa merge mode, this only works with TAXA_INC build type;//, GENO_EDIT}; //GENO_EDIT is not
    private GenotypeMergeRule mergeRule=null;
    private boolean isHDF5=false;
    private IHDF5Writer writer=null;
    private BuildType myBuildType;
    /*
    Builder for in memory taxa incremental
     */
    private GenotypeTableBuilder(PositionList positionList, GenotypeMergeRule mergeRule) {
        this.positionList=positionList;
        this.myBuildType=BuildType.TAXA_INC;
        this.mergeRule=mergeRule;
        if(mergeRule!=null) {
            this.isTaxaMerge=true;

        }
        incGeno=new ArrayList<>();
        incDepth=new ArrayList<>();
        incTaxonIndex=new HashMap<>();
        taxaListBuilder=new TaxaListBuilder();
    }

    /*
        Builder for in memory site incremental
     */
    private GenotypeTableBuilder(TaxaList taxaList) {
        this.taxaList=taxaList;
        this.myBuildType=BuildType.SITE_INC;
        incGeno=new ArrayList<>();
        posListBuilder=new PositionListBuilder();
    }

    /**
     * Creates a new HDF5 file if positionList is not null.  Opens an existing HDF5 File if positionList is null.
     * Merging is allowed depending on whether a mergeRule is included
     * that can used with TaxaIncremental addition.
     * @param hdf5File
     * @param positionList
     */
    private GenotypeTableBuilder(String hdf5File, PositionList positionList, GenotypeMergeRule mergeRule) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        config.dontUseExtendableDataTypes();
        writer = config.writer();
        if(HDF5Utils.doesGenotypeModuleExist(writer) && HDF5Utils.isHDF5GenotypeLocked(writer)) {
            writer.close();
            throw new UnsupportedOperationException("This file is locked for genotypic additions");
        }
        if(positionList!=null) {
            this.positionList=new PositionListBuilder(writer,positionList).build();  //create a new position list
            setupGenotypeTaxaInHDF5(writer);
        } else {
            this.positionList=PositionListBuilder.getInstance(writer);

        }
        this.mergeRule=mergeRule;
        if(mergeRule!=null) {
            this.isTaxaMerge=true;
        }
        this.myBuildType=BuildType.TAXA_INC;
        isHDF5=true;
        //taxaListBuilder=new TaxaListBuilder();
    }

    /**
     * Creates an in memory builder for addition by taxon.  Each taxon can only be added once,
     * i.e. merging is not possible
     * @param positionList The positions used for the builder
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList) {
        return new GenotypeTableBuilder(positionList,(GenotypeMergeRule)null);
    }

    /**
     * Creates an in memory builder for addition by taxon, which permits the merging of taxa.
     * @param positionList The positions used for the builder
     * @param mergeRule rules for merging identically named taxa
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList, GenotypeMergeRule mergeRule) {
        return new GenotypeTableBuilder(positionList,(GenotypeMergeRule)mergeRule);
    }

    /**
     * Creates a builder initialized with the Genotypes in a existing GenotypeTable.  The position list and initial taxa list is
     * derived from the positions, taxa, and genotypes already in the GenotypeTable.  The initial GenotypeTable is
     * not changed as it is immutable.
     * @param genotypeTable input genotype table
     * @param mergeRule rules for merging identically named taxa
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(GenotypeTable genotypeTable, GenotypeMergeRule mergeRule) {
        PositionList positionList=genotypeTable.positions();
        GenotypeTableBuilder gtb= new GenotypeTableBuilder(positionList, mergeRule);
        boolean hasDepth=genotypeTable.hasDepth();
        for (int i=0; i<genotypeTable.numberOfTaxa(); i++) {
            if(hasDepth) {
                gtb.addTaxon(genotypeTable.taxa().get(i),genotypeTable.genotypeAllSites(i), genotypeTable.depth().depthAllSitesByte(i));
            } else {
                gtb.addTaxon(genotypeTable.taxa().get(i),genotypeTable.genotypeAllSites(i));
            }
        }
        return gtb;
    }

    /**
     * Create a new taxa incremental HDF5 GenotypeTableBuilder
     * @param positionList the defined list of positions
     * @param newHDF5File  hdf5 file to be created
     * @return the builder to add taxa to
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList, String newHDF5File) {
        return new GenotypeTableBuilder(newHDF5File, positionList, null);
    }

    /**
     * Merges taxa to an existing HDF5 file.  The position list is derived from the positions already in the
     * existing HDF5 file.
     * @param existingHDFFile
     * @param mergeRule
     * @return
     */
    public static GenotypeTableBuilder mergeTaxaIncremental(String existingHDFFile, GenotypeMergeRule mergeRule) {
        return new GenotypeTableBuilder(existingHDFFile, null, mergeRule);
    }




    public static GenotypeTableBuilder getSiteIncremental(TaxaList taxaList) {
        return new GenotypeTableBuilder(taxaList);
    }

    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        if (genotype.numberOfSites() != positionList.numberOfSites()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + genotype.numberOfSites() + " doesn't equal number of sites in position list: " + positionList.numberOfSites());
        }
        if (genotype.numberOfTaxa() != taxaList.numberOfTaxa()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + genotype.numberOfTaxa() + " doesn't equal number of taxa in taaxa list: " + taxaList.numberOfTaxa());
        }
        return new CoreGenotypeTable(genotype, positionList, taxaList, siteScore, alleleDepth);
    }

    /**
     * Standard approach for creating a new Alignment
     * @param genotype
     * @param positionList
     * @param taxaList
     * @return new alignment
     */
    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList) {
        if (genotype.numberOfSites() != positionList.numberOfSites()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + genotype.numberOfSites() + " doesn't equal number of sites in position list: " + positionList.numberOfSites());
        }
        if (genotype.numberOfTaxa() != taxaList.numberOfTaxa()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + genotype.numberOfTaxa() + " doesn't equal number of taxa in taaxa list: " + taxaList.numberOfTaxa());
        }
        return new CoreGenotypeTable(genotype, positionList, taxaList);
    }

    /**
     * Creates a new HDF5 file alignment based on existing Genotype, PositionList, and TaxaList.
     * @param genotype
     * @param positionList
     * @param taxaList
     * @param hdf5File name of the file
     * @return alignment backed by new HDF5 file
     */
    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, String hdf5File) {
        if (genotype.numberOfSites() != positionList.numberOfSites()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + genotype.numberOfSites() + " doesn't equal number of sites in position list: " + positionList.numberOfSites());
        }
        if (genotype.numberOfTaxa() != taxaList.numberOfTaxa()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + genotype.numberOfTaxa() + " doesn't equal number of taxa in taaxa list: " + taxaList.numberOfTaxa());
        }
        GenotypeTableBuilder aB=GenotypeTableBuilder.getTaxaIncremental(positionList,hdf5File);
        for (int i=0; i<taxaList.numberOfTaxa(); i++) {
            aB.addTaxon(taxaList.get(i),genotype.genotypeAllSites(i));
        }
        return aB.build();
    }

    /**
     * Creates a new HDF5 file alignment based on an existing alignment.
     * @param a existing alignment
     * @param hdf5File name of the file
     * @return alignment backed by new HDF5 file
     */
    public static GenotypeTable getInstance(GenotypeTable a, String hdf5File) {
        return getInstance(a.genotypeMatrix(), a.positions(), a.taxa(), hdf5File);
    }

    public static GenotypeTable getInstance(String hdf5File) {
        IHDF5Reader reader=HDF5Factory.openForReading(hdf5File);
        TaxaList tL=new TaxaListBuilder().buildFromHDF5Genotypes(reader);
        PositionList pL=PositionListBuilder.getInstance(reader);
        GenotypeCallTable geno=GenotypeCallTableBuilder.buildHDF5(reader);
        AlleleDepth depth=AlleleDepthBuilder.getExistingHDF5Instance(reader);
        return GenotypeTableBuilder.getInstance(geno, pL, tL,null,depth);
    }

    public static GenotypeTable getInstanceOnlyMajorMinor(GenotypeTable alignment) {
        int numTaxa = alignment.numberOfTaxa();
        int numSites = alignment.numberOfSites();
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstance(numTaxa, numSites);
        byte[] majorAllele = new byte[64];
        byte[] minorAllele = new byte[64];
        for (int bigS = 0; bigS < numSites; bigS += 64) {
            int blockSize = Math.min(64, numSites - bigS);

            for (int s = 0; s < blockSize; s++) {
                majorAllele[s] = alignment.majorAllele(s + bigS);
                minorAllele[s] = alignment.minorAllele(s + bigS);
            }

            for (int t = 0; t < numTaxa; t++) {
                for (int s = 0; s < blockSize; s++) {
                    byte[] currentAlleles = alignment.genotypeArray(t, s + bigS);
                    if ((currentAlleles[0] != majorAllele[s]) && (currentAlleles[0] != minorAllele[s])) {
                        currentAlleles[0] = GenotypeTable.UNKNOWN_ALLELE;
                    }
                    if ((currentAlleles[1] != majorAllele[s]) && (currentAlleles[1] != minorAllele[s])) {
                        currentAlleles[1] = GenotypeTable.UNKNOWN_ALLELE;
                    }
                    builder.setBase(t, s, GenotypeTableUtils.getDiploidValue(currentAlleles[0], currentAlleles[1]));
                }
            }
        }
        return new CoreGenotypeTable(builder.build(), alignment.positions(), alignment.taxa());
    }

    public static GenotypeTable getHomozygousInstance(GenotypeTable alignment) {
        int numTaxa = alignment.numberOfTaxa();
        int numSites = alignment.numberOfSites();
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstance(numTaxa, numSites);
        //TODO this would be even faster to work through the SuperByteMatrix, as knowledge of site or taxa is not needed.
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                byte currGeno=alignment.genotype(t, s);
                if(GenotypeTableUtils.isHeterozygous(currGeno)) {
                    builder.setBase(t, s, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                } else {
                    builder.setBase(t, s, currGeno);
                }
            }

        }
        return new CoreGenotypeTable(builder.build(), alignment.positions(), alignment.taxa());
    }

    /**
     * Returns a taxa optimized version of a filtered alignment.  Only needed in performance critical situations
     * like imputation.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    public static GenotypeTable getGenotypeCopyInstance(FilterGenotypeTable alignment) {
        return copyGenotypeInstance(alignment);
    }

    /**
     * Returns a taxa optimized version of a combine alignment.  Only needed in performance critical situations
     * like imputation.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    public static GenotypeTable getGenotypeCopyInstance(CombineGenotypeTable alignment) {
        return copyGenotypeInstance(alignment);
    }

    /**
     * This is private as the method used by combine and filter alignment, as other data structures
     * are immutable and optimized that is unneeded except for these alignment types.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    private static GenotypeTable copyGenotypeInstance(GenotypeTable alignment) {
        int numTaxa = alignment.numberOfTaxa();
        int numSites = alignment.numberOfSites();
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstance(numTaxa, numSites);
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) { builder.setBase(t, s, alignment.genotype(t, s));}
        }
        return new CoreGenotypeTable(builder.build(), alignment.positions(), alignment.taxa());
    }

    public GenotypeTableBuilder addSite(Position pos, byte[] genos) {
        if((myBuildType!=BuildType.SITE_INC)||isHDF5) throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental and without HDF5");
        if(genos.length!=taxaList.numberOfTaxa()) throw new IndexOutOfBoundsException("Number of taxa and genotypes do not agree");
        posListBuilder.add(pos);
        incGeno.add(genos);
        return this;
    }

    public GenotypeTableBuilder addTaxon(Taxon taxon, byte[] genos) {
        return addTaxon(taxon, genos, null);
    }

    public GenotypeTableBuilder addTaxon(Taxon taxon, byte[] genos, byte[][] depth) {
        if(myBuildType!=BuildType.TAXA_INC) throw new IllegalArgumentException("addTaxon only be used with AlignmentBuilder.getTaxaIncremental");
        if(genos.length!=positionList.numberOfSites()) throw new IndexOutOfBoundsException("Number of sites and genotypes do not agree");
        if(isHDF5) {
            if(isTaxaMerge && HDF5Utils.doTaxonCallsExist(writer,taxon)) {
                mergeTaxonInHDF5(writer,taxon,genos,depth);
            } else {
                addTaxon(writer, taxon, genos, depth);
            }
        } else {
            if(isTaxaMerge && incTaxonIndex.containsKey(taxon)) {
                mergeTaxonInMemory(taxon,genos,depth);
            }  else {
                taxaListBuilder.add(taxon);
                incGeno.add(genos);
                incDepth.add(depth);
                incTaxonIndex.put(taxon,incGeno.size()-1);
            }
        }

        return this;
    }
    
    public GenotypeTableBuilder addTaxon(Taxon taxon, int[][] depths, byte[] genos) {  // reversed signature so it does clash with above statement: "return addTaxon(taxon, genos, null);"
        byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(depths);
        return addTaxon(taxon, genos, byteDepths);
    }

    private void mergeTaxonInMemory(Taxon taxon, byte[] genos, byte[][] depth) {
        int taxonIndex=incTaxonIndex.get(taxon);
        byte[] combGenos=new byte[genos.length];
        if(depth!=null) {
            byte[][] existingDepth=incDepth.get(taxonIndex);
//            System.out.println("ExistingDepth");
//            System.out.println(Arrays.deepToString(existingDepth));
//            System.out.println(Arrays.deepToString(depth));
            byte[][] combDepth=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][genos.length];
            byte[] currDepths=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (int site=0; site<combDepth[0].length; site++) {
                for (int allele=0; allele<combDepth.length; allele++) {
                    currDepths[allele]=combDepth[allele][site]=AlleleDepthUtil.addByteDepths(depth[allele][site],existingDepth[allele][site]);
                }
                combGenos[site]=mergeRule.callBasedOnDepth(currDepths);
            }
            incGeno.set(taxonIndex,combGenos);
            incDepth.set(taxonIndex,combDepth);
        } else {
            byte[] existingGenos=incGeno.get(taxonIndex);
            for (int site=0; site<combGenos.length; site++) {
                combGenos[site]=mergeRule.mergeCalls(genos[site],existingGenos[site]);
            }
            incGeno.set(taxonIndex,combGenos);
        }
    }

    public boolean isHDF5() {
       return isHDF5;
    }

    /*
    Set the builder so that when built it will sort the taxa
     */
    public GenotypeTableBuilder sortTaxa() {
        if(myBuildType!=BuildType.TAXA_INC) throw new IllegalArgumentException("sortTaxa can only be used with AlignmentBuilder.getTaxaIncremental");
        sortAlphabetically=true;
        return this;
    }

    /**
     * Finishes building the GenotypeTable.  For HDF5 files it locks the taxa and genotype modules so that cannot be
     * modified again.
     * @return a genotype table
     */
    public GenotypeTable build(){
        if(isHDF5) {
            String name=writer.getFile().getAbsolutePath();
            annotateHDF5File(writer);
            HDF5Utils.lockHDF5GenotypeModule(writer);
            HDF5Utils.lockHDF5TaxaModule(writer);
            writer.close();
            return getInstance(name);
        }
        switch (myBuildType) {
            case TAXA_INC: {
                TaxaList tl=(sortAlphabetically)?taxaListBuilder.sortTaxaAlphabetically().build():taxaListBuilder.build();
                GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(tl.numberOfTaxa(),positionList.numberOfSites());
                boolean hasDepth=(incDepth.size()==tl.numberOfTaxa() && incDepth.get(0)!=null);
                AlleleDepthBuilder adb=null;
                if(hasDepth) {adb=AlleleDepthBuilder.getNucleotideInstance(tl.numberOfTaxa(),positionList.numberOfSites());}
                for (int i=0; i<incGeno.size(); i++) {
                    gB.setBaseRangeForTaxon(i, 0, incGeno.get(incTaxonIndex.get(tl.get(i))));
                    if(hasDepth) {
                        adb.setDepth(i,incDepth.get(incTaxonIndex.get(tl.get(i))));
                    }
                }
                AlleleDepth ad=(hasDepth)?adb.build():null;
                return new CoreGenotypeTable(gB.build(), positionList, tl,null,ad);
            }
            case SITE_INC: {
                //TODO validate sort order, sort if needed
                PositionList pl=posListBuilder.build();
                GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(taxaList.numberOfTaxa(),pl.numberOfSites());
                for (int s=0; s<pl.numberOfSites(); s++) {
                    byte[] b=incGeno.get(s);
                    for (int t=0; t<b.length; t++) {
                        gB.setBase(t,s,b[t]);
                    }
                }
                return new CoreGenotypeTable(gB.build(), pl, taxaList);
            }
        }
        return null;
    }

    /**
     * Used to close an HDF5 GenotypeTableBuilder, when it will be reopened later and appended.  This file cannot be
     * used for other purposes in this unfinished state.
     */
    public void closeUnfinished(){
        if(isHDF5==false) throw new UnsupportedOperationException("Only a HDF5 GenotypeTableBuilder can be closed");
        taxaListBuilder=null;
        writer.close();
    }

    /*
    HDF5 Alignment section.
     */

    private synchronized void setupGenotypeTaxaInHDF5(IHDF5Writer writer) {
        HDF5Utils.createHDF5TaxaModule(writer);
        HDF5Utils.createHDF5GenotypeModule(writer);
        HDF5Utils.writeHDF5GenotypesMaxNumAlleles(writer,NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
        HDF5Utils.writeHDF5GenotypesRetainRareAlleles(writer,false);
        HDF5Utils.writeHDF5GenotypesNumTaxa(writer,0);
        HDF5Utils.writeHDF5GenotypesAlleleStates(writer,NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    /**
     * Code needed to add a Taxon to HDF5
     */
    private synchronized void addTaxon(IHDF5Writer myWriter, Taxon id, byte[] genotype, byte[][] depth) {
        boolean goodAdd=HDF5Utils.addTaxon(myWriter,id);
        if(goodAdd==false) throw new IllegalStateException("Taxon ["+id.getName()+"] already exists in the HDF5 file.  Duplicated taxa not allowed.");
        HDF5Utils.writeHDF5GenotypesCalls(myWriter,id.getName(),genotype);
        if(depth!=null) {
            if(depth.length!=6) throw new IllegalStateException("Just set A, C, G, T, -, + all at once");
            if(depth[0].length!=positionList.numberOfSites()) throw new IllegalStateException("Setting all depth in addTaxon.  Wrong number of sites");
            HDF5Utils.writeHDF5GenotypesDepth(myWriter,id.getName(),depth);
        }
    }

    private synchronized void mergeTaxonInHDF5(IHDF5Writer myWriter, Taxon id, byte[] genotype, byte[][] depth) {
        String[] existingFlowCellLanes = HDF5Utils.getTaxon(writer, id.getName()).getTextAnnotation("Flowcell_Lane");
        String[] newFlowCellLanes = id.getTextAnnotation("Flowcell_Lane");
        if (newFlowCellLanes.length > 0) {
            for (String existingFL : existingFlowCellLanes) {
                if (existingFL.equals(newFlowCellLanes[0])) 
                    throw new IllegalStateException("mergeTaxonInHDF5: Reads from flowcell_lane "+id.getTextAnnotation("Flowcell_Lane")[0]
                            +" previously added to taxon "+id.getName());
            }
            Taxon modifiedTaxon = new Taxon.Builder(HDF5Utils.getTaxon(writer, id.getName())).addAnno("Flowcell_Lane",newFlowCellLanes[0]).build();
            HDF5Utils.replaceTaxonAnnotations(myWriter, modifiedTaxon);
        }
        byte[] combGenos=new byte[genotype.length];
        if(depth!=null) {
            byte[][] existingDepth=HDF5Utils.getHDF5GenotypesDepth(myWriter,id.getName());
            byte[][] combDepth=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][genotype.length];
            byte[] currDepths=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (int site=0; site<combDepth[0].length; site++) {
                for (int allele=0; allele<combDepth.length; allele++) {
                    currDepths[allele]=combDepth[allele][site]=AlleleDepthUtil.addByteDepths(depth[allele][site],existingDepth[allele][site]);
                }
                combGenos[site]=mergeRule.callBasedOnDepth(currDepths);
            }
            HDF5Utils.replaceHDF5GenotypesCalls(myWriter,id.getName(),combGenos);
            HDF5Utils.replaceHDF5GenotypesDepth(myWriter,id.getName(),combDepth);
        } else {
            byte[] existingGenos=HDF5Utils.getHDF5GenotypesCalls(myWriter,id.getName());
            for (int site=0; site<combGenos.length; site++) {
                combGenos[site]=mergeRule.mergeCalls(genotype[site],existingGenos[site]);
            }
            HDF5Utils.replaceHDF5GenotypesCalls(myWriter,id.getName(),combGenos);
        }
    }

    /**
     * Annotates the HDF5 Genotype file with allele frequency information.  Can only be called on unlocked HDF5 files.
     * Currently, placed in the GenotypeTableBuilder as it still above genotypes, taxa, and sites.
     * @param writer
     */
    public static void annotateHDF5File(IHDF5Writer writer) {
       // int hdf5GenoBlock=writer.getIntAttribute(Tassel5HDF5Constants.DEFAULT_ATTRIBUTES_PATH, Tassel5HDF5Constants.BLOCK_SIZE);
        if(HDF5Utils.isHDF5GenotypeLocked(writer)) throw new UnsupportedOperationException("This is a locked HDF5 file");
        int hdf5GenoBlock=Tassel5HDF5Constants.BLOCK_SIZE;
        int sites=writer.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        TaxaList tL=new TaxaListBuilder().buildFromHDF5Genotypes(writer);
        int taxa=tL.numberOfTaxa();
        writer.setIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA,taxa);
        int[][] af=new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        byte[][] afOrder=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        float[] coverage=new float[taxa];
        float[] hets=new float[taxa];
        for (int taxon = 0; taxon < taxa; taxon++) {
            String basesPath = Tassel5HDF5Constants.getGenotypesCallsPath(tL.taxaName(taxon));
            byte[] genotype=writer.readByteArray(basesPath);
            int covSum=0;  //coverage of the taxon
            int hetSum=0;
            for (int s = 0; s < sites; s++) {
                byte[] b = GenotypeTableUtils.getDiploidValues(genotype[s]);
                if(b[0]<6) af[b[0]][s]++;
                if(b[1]<6) af[b[1]][s]++;
                if(GenotypeTableUtils.isHeterozygous(genotype[s])) hetSum++;
                if(genotype[s]!=GenotypeTable.UNKNOWN_DIPLOID_ALLELE) covSum++;
            }
            coverage[taxon]=(float)covSum/(float)sites;
            hets[taxon]=(float)hetSum/(float)covSum;
        }
        float[] maf=new float[sites];
        float[] paf=new float[sites];
        int baseMask=0xF;
        for (int s = 0; s < sites; s++) {
            int sum=0;
            int[] cntAndAllele=new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i]=(af[i][s]<<4)|(5-i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                sum+=af[i][s];
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort, there are faster ways
            //http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
            for (byte i = 0; i < 6; i++) {
                afOrder[5-i][s]=(cntAndAllele[i]>0xF)?((byte)(5-(baseMask&cntAndAllele[i]))):GenotypeTable.UNKNOWN_ALLELE;
            }
            if(afOrder[1][s]!=GenotypeTable.UNKNOWN_ALLELE) maf[s]=(float)af[afOrder[1][s]][s]/(float)sum;
            paf[s]=(float)sum/(float)(2*taxa);
        }
        writer.createGroup(Tassel5HDF5Constants.GENO_DESC);
        int chunk=(sites<hdf5GenoBlock)?sites:hdf5GenoBlock;
        writer.createIntMatrix(Tassel5HDF5Constants.ALLELE_CNT, 6, sites, 1, chunk, Tassel5HDF5Constants.intDeflation);
        writer.createByteMatrix(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 6, sites, 1, chunk, Tassel5HDF5Constants.intDeflation);
        writer.createFloatArray(Tassel5HDF5Constants.MAF, sites, chunk, Tassel5HDF5Constants.floatDeflation);
        writer.createFloatArray(Tassel5HDF5Constants.SITECOV, sites, chunk, Tassel5HDF5Constants.floatDeflation);
      //  writer.createGroup(Tassel5HDF5Constants.TAXA_DESC);
        chunk=(tL.numberOfTaxa()<hdf5GenoBlock)?tL.numberOfTaxa():hdf5GenoBlock;
        writer.createFloatArray(Tassel5HDF5Constants.TAXACOV, tL.numberOfTaxa(), chunk, Tassel5HDF5Constants.floatDeflation);
        writer.createFloatArray(Tassel5HDF5Constants.TAXAHET, tL.numberOfTaxa(), chunk, Tassel5HDF5Constants.floatDeflation);
        if(af[0].length>0) HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.ALLELE_CNT, writer, af[0].length, 1<<16, af);
        if(afOrder[0].length>0) HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.ALLELE_FREQ_ORD, writer, afOrder[0].length, 1<<16, afOrder);
        if(maf.length>0) HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.MAF, writer, maf.length, 1<<16, maf);
        if(paf.length>0) HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.SITECOV, writer, paf.length, 1<<16, paf);

        if(coverage.length>0) HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.TAXACOV, writer, coverage.length, 1<<16, coverage);
        if(hets.length>0) {
            System.out.println("Hets Length:"+hets.length);
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.TAXAHET, writer, hets.length, 1<<16, hets);
        }
    }

private enum BuildType{TAXA_INC, SITE_INC}
}
