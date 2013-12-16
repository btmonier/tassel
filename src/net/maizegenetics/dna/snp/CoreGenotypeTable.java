/*
 * CoreGenotypeTable
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.bit.DynamicBitStorage;
import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Map;

/**
 * Basic implementation of a {@link GenotypeTable}.  Use the GenotypeTableBuilder to construct.
 *
 * @see GenotypeTable
 * @see GenotypeTableBuilder
 *
 * @author Terry Casstevens
 */
public class CoreGenotypeTable implements GenotypeTable {

    private static final Logger myLogger = Logger.getLogger(CoreGenotypeTable.class);
    private final GenotypeCallTable myGenotype;
    private final Map<ALLELE_SORT_TYPE, BitStorage> myBitStorage = new HashMap<ALLELE_SORT_TYPE, BitStorage>();
    private final PositionList myPositionList;
    private final TaxaList myTaxaList;
    private final SiteScore mySiteScore;
    private final AlleleDepth myAlleleDepth;
    private final int mySiteCount;
    private final int myTaxaCount;

    CoreGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        //todo need check dimensions
        myGenotype = genotype;
        myPositionList = positionList;
        myTaxaList = taxaList;
        mySiteScore = siteScore;
        myAlleleDepth = alleleDepth;
        mySiteCount = myPositionList.numberOfSites();
        myTaxaCount = myTaxaList.numberOfTaxa();
    }

    CoreGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList) {
        this(genotype, positionList, taxaList, null, null);
    }

    @Override
    public GenotypeCallTable genotypeMatrix() {
        return myGenotype;
    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.genotype(taxon, site);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return myGenotype.genotypeArray(taxon, site);
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.genotype(taxon, myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        return myGenotype.genotypeForAllTaxa(site);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllSites(taxon, alleleNumber);
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForSitesBlock(taxon, alleleNumber, startBlock, endBlock);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllSites(taxon, firstParent, alleleNumber);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllTaxa(site, firstParent, alleleNumber);
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForSitesBlock(taxon, firstParent, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myGenotype.genotypeAsString(taxon, site);
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return myGenotype.genotypeAsStringRow(taxon);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        return myGenotype.genotypeAsStringArray(taxon, site);
    }

    @Override
    public byte referenceGenotype(int site) {
        return myPositionList.referenceGenotype(site);
    }

    @Override
    public byte[] referenceGenotypes(int startSite, int endSite) {
        return myPositionList.referenceGenotypes(startSite, endSite);
    }

    @Override
    public byte[] referenceGenotypeForAllSites() {
        return myPositionList.referenceGenotypeForAllSites();
    }

    @Override
    public boolean hasReference() {
        return myPositionList.hasReference();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int heterozygousCount(int site) {
        return myGenotype.heterozygousCount(site);
    }

    @Override
    public PositionList positions() {
        return myPositionList;
    }

    @Override
    public String siteName(int site) {
        return myPositionList.siteName(site);
    }

    @Override
    public int numberOfSites() {
        return mySiteCount;
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        return myPositionList.chromosomeSiteCount(chromosome);
    }

    @Override
    public int[] startAndEndOfChromosome(Chromosome chromosome) {
        return myPositionList.startAndEndOfChromosome(chromosome);
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaCount;
    }

    @Override
    public int chromosomalPosition(int site) {
        return myPositionList.chromosomalPosition(site);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
    }

    @Override
    public int[] physicalPositions() {
        return myPositionList.physicalPositions();
    }

    @Override
    public String chromosomeName(int site) {
        return myPositionList.chromosomeName(site);
    }

    @Override
    public Chromosome chromosome(int site) {
        return myPositionList.chromosome(site);
    }

    @Override
    public Chromosome chromosome(String name) {
        return myPositionList.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myPositionList.chromosomes();
    }

    @Override
    public int numChromosomes() {
        return myPositionList.numChromosomes();
    }

    @Override
    public int[] chromosomesOffsets() {
        return myPositionList.chromosomesOffsets();
    }

    @Override
    public float siteScore(int taxon, int site) {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScore: This Alignment has no Site Scores.");
        }
        return mySiteScore.siteScore(taxon, site);
    }

    @Override
    public float[][] siteScores() {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScores: This Alignment has no Site Scores.");
        }
        return mySiteScore.siteScores();
    }

    @Override
    public boolean hasSiteScores() {
        if (mySiteScore == null) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public SITE_SCORE_TYPE siteScoreType() {
        return mySiteScore.siteScoreType();
    }

    @Override
    public int indelSize(int site) {
        return myPositionList.indelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return myPositionList.isIndel(site);
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public byte majorAllele(int site) {
        return myGenotype.majorAllele(site);
    }

    @Override
    public String majorAlleleAsString(int site) {
        return myGenotype.majorAlleleAsString(site);
    }

    @Override
    public byte minorAllele(int site) {
        return myGenotype.minorAllele(site);
    }

    @Override
    public String minorAlleleAsString(int site) {
        return myGenotype.minorAlleleAsString(site);
    }

    @Override
    public byte[] minorAlleles(int site) {
        return myGenotype.minorAlleles(site);
    }

    @Override
    public byte[] alleles(int site) {
        return myGenotype.alleles(site);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        return myGenotype.minorAlleleFrequency(site);
    }

    @Override
    public double majorAlleleFrequency(int site) {
        return myGenotype.majorAlleleFrequency(site);
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.taxaName(index);
    }

    @Override
    public String genomeVersion() {
        return myPositionList.genomeVersion();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myPositionList.isPositiveStrand(site);
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        return new GenotypeTable[]{this};
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        return myGenotype.genosSortedByFrequency(site);
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] alleleDefinitions() {
        return myGenotype.alleleDefinitions();
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myGenotype.alleleDefinitions(site);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return myGenotype.genotypeAsString(site, value);
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myGenotype.diploidAsString(site, value);
    }

    @Override
    public int maxNumAlleles() {
        return myGenotype.maxNumAlleles();
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        return myGenotype.totalGametesNonMissingForSite(site);
    }

    @Override
    public int totalNonMissingForSite(int site) {
        return myGenotype.totalNonMissingForSite(site);
    }

    @Override
    public int minorAlleleCount(int site) {
        return myGenotype.minorAlleleCount(site);
    }

    @Override
    public int majorAlleleCount(int site) {
        return myGenotype.majorAlleleCount(site);
    }

    @Override
    public Object[][] genoCounts() {
        return myGenotype.genoCounts();
    }

    @Override
    public Object[][] majorMinorCounts() {
        return myGenotype.majorMinorCounts();
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        return myGenotype.totalGametesNonMissingForTaxon(taxon);
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        return myGenotype.heterozygousCountForTaxon(taxon);
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        return myGenotype.totalNonMissingForTaxon(taxon);
    }

    @Override
    public byte[] depthForAlleles(int taxon, int site) {
        return myAlleleDepth.getDepthForAlleles(taxon, site);
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        switch (scope) {
            case Frequency:
                return alleles(site);
            default:
                myLogger.warn("getAllelesByScope: Unsupported type: " + scope);
                return null;
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxaBySortType(ALLELE_SORT_TYPE type, int site, int alleleNumber) {
        return bitStorage(type).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public BitStorage bitStorage(ALLELE_SORT_TYPE scopeType) {

        BitStorage result = myBitStorage.get(scopeType);
        if (result != null) {
            return result;
        }

        switch (scopeType) {
            case Frequency:
                result = new DynamicBitStorage(myGenotype, scopeType, myGenotype.majorAlleleForAllSites(), myGenotype.minorAlleleForAllSites());
                break;
            default:
                myLogger.warn("getBitStorage: Unsupported type: " + scopeType);
                return null;
        }

        myBitStorage.put(scopeType, result);
        return result;
    }
}
