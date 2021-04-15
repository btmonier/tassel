package net.maizegenetics.analysis.gbs.v2;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Set;
import org.ini4j.Ini;
import org.ini4j.Profile;

/**
 *
 * @author dpavelec
 */
public class EnzymeList {
    
    /**
     * Looks for "enzymes.ini" config file in the lib path. 
     * Loads default Enzymes if not found 
     */
    public EnzymeList() {
        String parentPath = getJarPath(Ini.class).getPath();
        File iniFile = new File(parentPath + "/enzymes.ini");
        if (!loadFromFile(iniFile)) {
            loadDefaults();
        }
    }

    /**
     * Custom Config file
     * @param iniFile 
     */
    public EnzymeList(File iniFile) {
        if (!loadFromFile(iniFile)) {
            loadDefaults();
        }
    }

    private synchronized boolean loadFromFile(File iniFile) {
        if(map==null){
            map = new LinkedHashMap();
        }
        try {
            Ini ini = new Ini(iniFile);
            if (!iniFile.exists()) {
                return false;
            }
            ini.load();
            
            //List<Profile.Section> all = ini.getAll(ini);
            Set<String> keySet = ini.keySet();
            String empty = null;
            for(String s : keySet){
                Profile.Section e = (Profile.Section) ini.get(s);
                String name = e.get("name", empty);
                String initialCutSiteRemnant = e.get("initialCutSiteRemnant", empty);
                String likelyReadEnd = e.get("likelyReadEnd", empty);
                String readLenString = e.get("readEndCutSiteRemnantLength", "-1");
                int readEndCutSiteRemnantLength = Integer.parseInt(readLenString);
                if(name==null|initialCutSiteRemnant==null|likelyReadEnd==null|readEndCutSiteRemnantLength<0){
                    System.err.println("WARNING! Cannot load Enzyme " + s);
                    continue;
                }else{
                    map.put(s.trim().toUpperCase(), new Enzyme(name, initialCutSiteRemnant.split(","), likelyReadEnd.split(","), readEndCutSiteRemnantLength));
                }
            }
        } catch (IOException ex) {
            System.err.println("ERROR! Cannot load Enzyme List -- " + iniFile.getPath());
            return false;
        }
        return true;
    }
    
    private void loadDefaults() {
        if(map==null){
            map = new LinkedHashMap();
        }
        map.put("APEKI", new Enzyme("ApeKI",
                new String[]{"CAGC", "CTGC"},
                new String[]{"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                4
        ));
        map.put("PSTI", new Enzyme("PstI",
                new String[]{"TGCAG"},
                new String[]{"CTGCAG", "CTGCAAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                5));
        
        map.put("NSII", new Enzyme("NsiI",
                new String[]{"TGCAA"},
                new String[]{"CTGCAA", "CTGCAAGAA"}, // full cut site (from partial digest or chimera) or common adapter start
                5));

        map.put("ECOT22I", new Enzyme("EcoT22I",
                new String[]{"TGCAT"},
                new String[]{"ATGCAT", "ATGCAAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                5));

        map.put("PASI", new Enzyme("PasI",
                new String[]{"CAGGG", "CTGGG"},
                new String[]{"CCCAGGG", "CCCTGGG", "CCCTGAGAT", "CCCAGAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                5));

        map.put("HPAII", new Enzyme("HpaII",
                new String[]{"CGG"},
                new String[]{"CCGG", "CCGAGATCGG"}, // full cut site (from partial digest or chimera) or common adapter start
                3));

        map.put("MSPI", new Enzyme("MspI", // MspI and HpaII are isoschizomers (same recognition seq and overhang)
                new String[]{"CGG"},
                new String[]{"CCGG", "CCGAGATCGG"}, // full cut site (from partial digest or chimera) or common adapter start
                3));

        map.put("PSTI-ECOT22I", new Enzyme("PstI-EcoT22I",
                new String[]{"TGCAG", "TGCAT"},
                new String[]{"ATGCAT", "CTGCAG", "CTGCAAGAT", "ATGCAAGAT"}, // look for EcoT22I site, PstI site, or common adapter for PstI/EcoT22I
                5));

        map.put("PSTI-MSPI", new Enzyme("PstI-MspI",
                new String[]{"TGCAG"},
                // corrected, change from  CCGAGATC to CCGCTCAGG, as Y adapter was used for MspI  -QS
                new String[]{"CCGG", "CTGCAG", "CCGCTCAGG"}, // look for MspI site, PstI site, or common adapter for MspI
                3));
        map.put("NSII-MSPI", new Enzyme("NsiI-MspI",
                new String[]{"TGCAA"},
                // corrected, change from  CCGAGATC to CCGCTCAGG, as Y adapter was used for MspI  -QS
                new String[]{"CCGG", "CTGCAA", "CCGCTCAGG"}, // look for MspI site, PstI site, or common adapter for MspI
                3));

        map.put("PSTI-TAQI", new Enzyme("PstI-TaqI",
                // corrected, change from  TCGAGATC to TCGCTCAGG, as Y adapter was used for TaqI  -QS
                new String[]{"TGCAG"},
                new String[]{"TCGA", "CTGCAG", "TCGCTCAGG"}, // look for TaqI site, PstI site, or common adapter for TaqI
                3));

        map.put("PAER7I-HHAI", new Enzyme("PaeR7I-HhaI",
                // Requested by Ye, Songqing, use same Y adapter as Polland paper  -QS
                new String[]{"TCGAG"},
                new String[]{"GCGC", "CTCGAG", "GCGAGATC"}, // look for HhaI site, PaeR7I site, or common adapter for HhaI
                3));

        map.put("SBFI-MSPI", new Enzyme("SbfI-MspI",
                new String[]{"TGCAGG"},
                // corrected, change from  CCGAGATC to CCGCTCAGG, as Y adapter was used for MspI  -QS
                new String[]{"CCGG", "CCTGCAGG", "CCGCTCAGG"}, // look for MspI site, SbfI site, or common adapter for MspI
                3));

        map.put("ASISI-MSPI", new Enzyme("AsiSI-MspI",
                new String[]{"ATCGC"},
                // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
                new String[]{"CCGG", "GCGATCGC", "CCGCTCAGG"}, // look for MspI site, AsiSI site, or common adapter for MspI
                3));

        map.put("BSSHII-MSPI", new Enzyme("BssHII-MspI",
                new String[]{"CGCGC"},
                // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
                new String[]{"CCGG", "GCGCGC", "CCGCTCAGG"}, // look for MspI site, BssHII site, or common adapter for MspI
                3));

        map.put("FSEI-MSPI", new Enzyme("FseI-MspI",
                new String[]{"CCGGCC"},
                // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
                new String[]{"CCGG", "GGCCGGCC", "CCGCTCAGG"}, // look for MspI site, FseI site, or common adapter for MspI
                3));

        map.put("SALI-MSPI", new Enzyme("SalI-MspI",
                new String[]{"TCGAC"},
                // likelyReadEnd for common adapter is CCGCTCAGG, as the Poland et al.(2012) Y adapter was used for MspI
                new String[]{"CCGG", "GTCGAC", "CCGCTCAGG"}, // look for MspI site, SalI site, or common adapter for MspI
                3));

        map.put("APOI", new Enzyme("ApoI",
                new String[]{"AATTC", "AATTT"}, // corrected from {"AATTG","AATTC"} by Jeff Glaubitz on 2012/09/12
                new String[]{"AAATTC", "AAATTT", "GAATTC", "GAATTT", "AAATTAGAT", "GAATTAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                5));

        map.put("BAMHI", new Enzyme("BamHI",
                new String[]{"GATCC"},
                new String[]{"GGATCC", "GGATCAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                //            new String[]{"GGATCC", "AGATCGGAA", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"}, // <-- corrected from this by Jeff Glaubitz on 2012/09/12
                5));

        map.put("MSEI", new Enzyme("MseI",
                new String[]{"TAA"},
                new String[]{"TTAA", "TTAAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                3));

        map.put("SAU3AI", new Enzyme("Sau3AI",
                new String[]{"GATC"},
                new String[]{"GATC", "GATCAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                4));

        map.put("NDEI", new Enzyme("NdeI",
                new String[]{"TATG"},
                new String[]{"CATATG", "CATAAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                4));

        map.put("HINPLI", new Enzyme("HinP1I",
                new String[]{"CGC"},
                new String[]{"GCGC", "GCGAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                3));

        map.put("SBFI", new Enzyme("SbfI",
                new String[]{"TGCAGG"},
                new String[]{"CCTGCAGG", "CCTGCAAGAT"}, // full cut site (from partial digest or chimera) or common adapter start
                6));

        map.put("RBSTA", new Enzyme("RBSTA",
                new String[]{"TA"},
                new String[]{"TTAA", "GTAC", "CTAG", "TTAAGAT", "GTAAGAT", "CTAAGAT"}, // full cut site (from partial digest or chimera) of MseI, CVIQi, XspI or common adapter start
                3));

        map.put("RBSCG", new Enzyme("RBSCG",
                new String[]{"CG"},
                new String[]{"CCGC", "TCGA", "GCGC", "CCGG", "ACGT", "CCGAGAT", "TCGAGAT", "GCGAGAT", "ACGAGAT"}, // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
                3));
    }

    /**
     * Gets the file location of the jar file containing the specified class.
     * <br>NOTE: may not work with remote files.
     *
     * @param clazz class to find
     * @return the location of the jar
     */
        static File getJarPath(Class clazz) {
        try {
            URL location = clazz.getProtectionDomain().getCodeSource().getLocation();
            File f = new File(location.toURI().getPath());
            if (!f.exists()) {
                File nf = new File(location.getPath());
                if (nf.exists()) {
                    f = nf;
                }
            }
            if (f.getName().toLowerCase().endsWith(".jar")) {
                return f.getParentFile();
            }
            return f;
        } catch (URISyntaxException ex) {
            return new File(clazz.getProtectionDomain().getCodeSource().getLocation().getPath());
        }
    }

    /**
     * Replaces the previous new GBSEnzyme(String enzyme)
     * @param name name of the enzyme
     * @return 
     */
    public synchronized Enzyme getEnzyme(String name) {
        if(map==null){
            loadDefaults();
        }
        return map.get(name.trim().toUpperCase());
    }
    
    /**
     * Add enzyme to a a configuration file
     * @param ini the config file 
     * @param e The enzyme to add
     */
    static void add(Ini ini,Enzyme e) {
        String key = e.name.trim().toUpperCase();
        Profile.Section section = ini.containsKey(key) ? (Profile.Section) ini.get(key) : ini.add(key);
        section.put("name", e.name);
        //initialCutSiteRemnant
        StringBuilder builder = new StringBuilder(e.initialCutSiteRemnant[0]);
        for (int i = 1; i < e.initialCutSiteRemnant.length; ++i) {
            builder.append(",");
            builder.append(e.initialCutSiteRemnant[i]);
        }
        section.put("initialCutSiteRemnant", builder.toString());
        //likelyReadEnd
        builder = new StringBuilder(e.likelyReadEnd[0]);
        for (int i = 1; i < e.likelyReadEnd.length; ++i) {
            builder.append(",");
            builder.append(e.likelyReadEnd[i]);
        }
        section.put("likelyReadEnd", builder.toString());
        section.put("readEndCutSiteRemnantLength", e.readEndCutSiteRemnantLength);
    }
    
    /**
     * Prints the list of Enzymes loaded in this object
     * @param print 
     * @throws IOException 
     */
    public void printEnzymes(PrintStream print) throws IOException{
        Ini ini = new Ini();
        ini.getConfig().setLineSeparator("\n");
        if(this.map==null){
            loadDefaults();
        }
        for(String key : this.map.keySet()){
            Enzyme e = this.map.get(key);
            add(ini,e);
        }
        ini.store(print);
    }
    
    public static void main(String[] args) throws IOException{
        new EnzymeList().printEnzymes(System.out);
    }

    private LinkedHashMap<String, Enzyme> map;
    /**
     * Replaces the previous GBSEnzyme class
     */
    public static class Enzyme {
        private Enzyme(
                String name,
                String[] initialCutSiteRemnant,
                String[] likelyReadEnd,
                int readEndCutSiteRemnantLength) {
            this.name = name;
            this.initialCutSiteRemnant = initialCutSiteRemnant;
            this.likelyReadEnd = likelyReadEnd;
            this.readEndCutSiteRemnantLength = readEndCutSiteRemnantLength;
        }
        public final String name;
        public final String[] initialCutSiteRemnant;
        public final String[] likelyReadEnd;
        public final int readEndCutSiteRemnantLength;
        
        /** 
         * This methods is for convenience based on old code. Use Enzyme.name
         * @deprecated 
         */
        public String enzyme() {return name;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.initialCutSiteRemnant
         * @deprecated 
         */
        public String[] initialCutSiteRemnant() {return initialCutSiteRemnant;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.likelyReadEnd
         * @deprecated 
         */
        public String[] likelyReadEnd() {return likelyReadEnd;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.readEndCutSiteRemnantLength
         * @deprecated 
         */
        public int readEndCutSiteRemnantLength() {return readEndCutSiteRemnantLength;}
    }
}
