package net.maizegenetics.util;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.taxa.Taxon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public final class HDF5Utils {

    private HDF5Utils() {
    }

    public static boolean isTASSEL4HDF5Format(IHDF5Reader reader) {
        return reader.exists(HapMapHDF5Constants.LOCI);
    }

 //TAXA Module
     public static void createHDF5TaxaModule(IHDF5Writer h5w) {
         h5w.createGroup(Tassel5HDF5Constants.TAXA_MODULE);
         h5w.setBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED,false);
     }

    public static void lockHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED,true);
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, String taxonName){
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxonName));
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, Taxon taxon){
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxon.getName()));
    }

    public static boolean isTaxaLocked(IHDF5Reader reader){
        return reader.getBooleanAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH,Tassel5HDF5Constants.TAXA_LOCKED);
    }

    public static boolean addTaxon(IHDF5Writer h5w, Taxon taxon) {
        if(isTaxaLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String path=Tassel5HDF5Constants.getTaxonPath(taxon.getName());
        if(h5w.exists(path)) return false;
        h5w.createGroup(path);
        SetMultimap<String, String> annoMap=taxon.getAnnotationAsMap();
        for (String keys : annoMap.keys()) {
            String s=Joiner.on(",").join(annoMap.get(keys));
            h5w.setStringAttribute(path,keys,s);
        }
        return true;
    }

    public static Taxon getTaxon(IHDF5Reader reader, String taxonName) {
        String taxonPath=Tassel5HDF5Constants.getTaxonPath(taxonName);
        if(!reader.exists(taxonPath)) return null;
        Taxon.Builder tb=new Taxon.Builder(taxonName);
        for (String a : reader.getAllAttributeNames(taxonPath)) {
            for(String s: Splitter.on(",").split(reader.getStringAttribute(taxonPath,a))) {
                tb.addAnno(a,s);
            }
        }
        return tb.build();
    }

    public static List<String> getAllTaxaNames(IHDF5Reader reader) {
        List<String> taxaNames=new ArrayList<>();
        List<HDF5LinkInformation> fields = reader.getAllGroupMemberInformation(Tassel5HDF5Constants.TAXA_MODULE, true);
        for (HDF5LinkInformation is : fields) {
            taxaNames.add(is.getName());
        }
        return taxaNames;
    }


    public static void writeHDF5TaxaNumTaxa(IHDF5Writer h5w, int numTaxa) {
        h5w.setIntAttribute(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_NUM_TAXA, numTaxa);
    }


 //GENOTYPE Module

    public static int getHDF5GenotypeTaxaNumber(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_NUM_TAXA);
    }

    public static void createHDF5GenotypeModule(IHDF5Writer h5w) {
        if(h5w.exists(Tassel5HDF5Constants.GENOTYPES_MODULE)) throw new UnsupportedOperationException("Genotypes module already exists in HDF5 file");
        h5w.createGroup(Tassel5HDF5Constants.GENOTYPES_MODULE);
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED,false);
    }

    public static void lockHDF5GenotypeModule(IHDF5Writer h5w) {
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED,true);
    }

    public static boolean isHDF5GenotypeLocked(IHDF5Reader reader){
        return reader.getBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH,Tassel5HDF5Constants.GENOTYPES_LOCKED);
    }

    public static void writeHDF5GenotypesMaxNumAlleles(IHDF5Writer h5w, int maxNumAlleles) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, maxNumAlleles);
    }

    public static void writeHDF5GenotypesRetainRareAlleles(IHDF5Writer h5w, boolean retainRareAlleles) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setBooleanAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_RETAIN_RARE_ALLELES, retainRareAlleles);
    }

    public static void writeHDF5GenotypesNumTaxa(IHDF5Writer h5w, int numTaxa) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setIntAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA, numTaxa);
    }

    public static void writeHDF5GenotypesScoreType(IHDF5Writer h5w, String scoreType) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        h5w.setStringAttribute(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, scoreType);
    }

    public static void writeHDF5GenotypesAlleleStates(IHDF5Writer h5w, String[][] aEncodings) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        int numEncodings = aEncodings.length;
        int numStates = aEncodings[0].length;
        MDArray<String> alleleEncodings = new MDArray<>(String.class, new int[]{numEncodings, numStates});
        for (int s = 0; s < numEncodings; s++) {
            for (int x = 0; x < numStates; x++) {
                alleleEncodings.set(aEncodings[s][x], s, x);
            }
        }
        h5w.createStringMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, 100, new int[]{numEncodings, numStates});
        h5w.writeStringMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, alleleEncodings);
    }

    public static byte[] getHDF5GenotypesCalls(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        return reader.readAsByteArray(callsPath);
    }

    public static void writeHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if(h5w.exists(callsPath)) throw new IllegalStateException("Taxa Calls Already Exists:"+taxon);
        h5w.createByteArray(callsPath, calls.length, Math.min(Tassel5HDF5Constants.BLOCK_SIZE,calls.length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static void replaceHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if(!h5w.exists(callsPath)) throw new IllegalStateException("Taxa Calls Do Not Already Exists to replace");
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static byte[][] getHDF5GenotypesDepth(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        return reader.readByteMatrix(callsPath);
    }

    public static void writeHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if(h5w.exists(callsPath)) throw new IllegalStateException("Taxa Depth Already Exists:"+taxon);
        h5w.createByteMatrix(callsPath, depth.length, depth[0].length, 6, Math.min(Tassel5HDF5Constants.BLOCK_SIZE,depth.length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }

    public static void replaceHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if(isHDF5GenotypeLocked(h5w)==true) throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if(!h5w.exists(callsPath)) throw new IllegalStateException("Taxa Depth Does Not Already Exists to Replace");
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }




//    Positions/numSites
    public static int getHDF5PositionNumber(IHDF5Reader reader){
        return reader.getIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH,Tassel5HDF5Constants.POSITION_NUM_SITES);
    }

    public static void createHDF5PositionModule(IHDF5Writer h5w) {
        h5w.createGroup(Tassel5HDF5Constants.POSITION_MODULE);
    }

    public static void writeHDF5PositionNumSite(IHDF5Writer h5w, int numSites) {
        h5w.setIntAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES, numSites);
    }

  //Writers for these should also be implemented, but there are some data scale issue that they are written in blocks.
    //see public PositionListBuilder(IHDF5Writer h5w, PositionList a)
    //    Positions/Positions
//    Positions/Chromosome
//    Positions/ChromosomeIndices
//    Positions/SnpIds



    /**
     *
     * @param objectPath
     * @param myWriter
     * @param objMaxLength
     * @param blockSize
     * @param val
     */
    public static void writeHDF5EntireArray(String objectPath, IHDF5Writer myWriter, int objMaxLength, int blockSize, Object val) {
        int blocks = ((objMaxLength - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = ((objMaxLength - startPos) > blockSize) ? blockSize : objMaxLength - startPos;
            if (val instanceof byte[][]) {
                byte[][] oval = (byte[][]) val;
                byte[][] sval = new byte[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof int[][]) {
                int[][] oval = (int[][]) val;
                int[][] sval = new int[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof byte[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((byte[]) val, startPos, startPos + length));
            } else if (val instanceof float[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((float[]) val, startPos, startPos + length));
            } else if (val instanceof int[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((int[]) val, startPos, startPos + length));
            } else if (val instanceof String[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((String[]) val, startPos, startPos + length));
            }
        }
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param blockSize
     * @param block
     * @param val
     */
    public static void writeHDF5Block(String objectPath, IHDF5Writer myWriter, int blockSize, int block, Object val) {
        int startPos = block * blockSize;
        if (val instanceof byte[][]) {
            byte[][] bval = (byte[][]) val;
            myWriter.writeByteMatrixBlockWithOffset(objectPath, bval, bval.length, bval[0].length, 0l, (long) startPos);
        } else if (val instanceof byte[]) {
            byte[] fval = (byte[]) val;
            myWriter.writeByteArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof float[]) {
            float[] fval = (float[]) val;
            myWriter.writeFloatArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[]) {
            int[] fval = (int[]) val;
            myWriter.writeIntArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[][]) {
            int[][] ival = (int[][]) val;
            myWriter.writeIntMatrixBlockWithOffset(objectPath, ival, ival.length, ival[0].length, 0l, (long) startPos);
        } else  if (val instanceof String[]) {
            String[] sval = (String[]) val;
            myWriter.writeStringArrayBlockWithOffset(objectPath, sval, sval.length, (long) startPos);
        }
    }

}
