/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import java.awt.Frame;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;

import com.google.common.collect.*;
import com.google.common.primitives.Ints;
import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;
import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.Tuple;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.log4j.Logger;

import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.tag.RepGenDataWriter;
import net.maizegenetics.dna.tag.RepGenSQLite;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TaxaDistribution;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;

/**
 * This class takes a rAmpSeq (formerly RepGen) database and for
 * each tag in the tag table, performs the following tag-tag
 * correlations based on the taxa distribution for each tag
 *   (1) tag-tag Pearson's correlation
 *   (2) tag-tag Spearman's correlation
 *   (3) tag-tag presence/absence correlations
 *   (4) r-squared
 *   
 *   The vectors presented to the analysis methods represent
 *   a list of taxa and the number of times the tag was seen
 *   in that taxa.  The presence/absence vectors have a 1 or 0
 *   as values in each slot.
 * @author lcj34
 *
 */
public class RepGenLDAnalysisPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(RepGenLDAnalysisPlugin.class);
    
    private PluginParameter<String> myInputDB = new PluginParameter.Builder<String>("db", null, String.class).guiName("Input DB").required(true).inFile()
            .description("Input database file with tags and taxa distribution").build();
    private PluginParameter<Integer> minTaxa = new PluginParameter.Builder<>("minTaxa", 50, Integer.class).guiName("Min Taxa for RSquared")
            .description("Minimum number of taxa that must be present for R-squared to be calculated.").build();
    private PluginParameter<Integer> minTaxaDepth = new PluginParameter.Builder<>("minTaxaDepth", 5000, Integer.class).guiName("Min Taxa for RSquared")
            .description("Minimum number of taxa that must be present for R-squared to be calculated.").build();
    private PluginParameter<TaxaList> myTaxaList = new PluginParameter.Builder<>("taxaList", null, TaxaList.class).guiName("Taxa list to test")
            .description("Minimum number of taxa that must be present for R-squared to be calculated.").build();
    private PluginParameter<Integer> minTagCount = new PluginParameter.Builder<>("minTagCount", 100, Integer.class).guiName("Min tag count after filtering")
            .description("Minimum number of taxa that must be present for R-squared to be calculated.").build();
    private PluginParameter<String> myProcessName = new PluginParameter.Builder<>("methodName", null, String.class).required(true).guiName("Name for process")
            .description("Name of the process ").build();
    public RepGenLDAnalysisPlugin() {
        super(null, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public RepGenLDAnalysisPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    // THis Map defines segregation expectations as punnett square values.
    // Each defined string gets 2 values in its double[][] array: One for the expected
    // probability of the upper row of the square, one containing the expected
    // probabilities for the lower row in the square.  These are used when
    // calculating F2 deviations below.
    private static final Map<String, double[][]> segExpectations;
    static {
        segExpectations=new HashMap<>();
        segExpectations.put("RIL",new double[][]{{0,0.495},{0.495,0.01}});
        segExpectations.put("F2",new double[][]{{0,0.25},{0.25,0.5}});
    }

    
    @Override
    public void postProcessParameters() {

        if (myInputDB.isEmpty() || !Files.exists(Paths.get(inputDB()))) {
            throw new IllegalArgumentException("RepGenLDAnalysisPlugin: postProcessParameters: Input DB not set or found");
        }
    }
    
    @Override
    public DataSet processData(DataSet input) {
        long totalTime = System.nanoTime();
        long time=System.nanoTime();
 
        try {           
            System.out.println("RepGenLDAnalysis:processData begin, get all tags/taxadist from db");
            RepGenDataWriter repGenData=new RepGenSQLite(inputDB());

            Map <Tag, TaxaDistribution> tagTaxaMap = repGenData.getAllTagsTaxaMap();
            Map <Tag, double[]> tagDepthMap = new HashMap<>();
            Map <Tag, Queue<Tuple<Double,Tag>>> tagQueueMap = new ConcurrentHashMap<>();

            System.out.println("Time to get all tags with taxa from db: " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
            time = System.nanoTime();
            System.out.println("\nStart processing tag correlations.  Number of tags in db: " + tagTaxaMap.keySet().size());
            Set<Tag> tagSet = tagTaxaMap.keySet();
            List<Tag> tagList = new ArrayList<>(tagSet);
            //MinMaxPriorityQueue<Tuple<Double,Integer>>[] corrQueList = new MinMaxPriorityQueue[tagList.size()];

            double[][] depthsTagTaxa=new double[tagList.size()][];
//            Ordering<Tuple<Double,Integer>> byCorrOrdering = new Ordering<Tuple<Double,Integer>>() {
//                public int compare(Tuple<Double,Integer>left, Tuple<Double,Integer> right) {
//                    return Double.compare(left.x,right.x);
//                }
//            };
            Ordering<Tuple<Double,Tag>> byCorrOrderingTag = new Ordering<Tuple<Double,Tag>>() {
                public int compare(Tuple<Double,Tag>left, Tuple<Double,Tag> right) {
                    return Double.compare(left.x,right.x);
                }
            };

            TaxaList taxons=repGenData.getTaxaList();
            final int[] depthByTaxon=new int[taxons.size()];
            tagTaxaMap.values().stream()
                .map(TaxaDistribution::depths)
                .forEach(d -> {
                    for (int t = 0; t < depthByTaxon.length; t++) {
                        depthByTaxon[t]+=d[t];
                    }
                });

            int minCombDept=5;
            System.out.println(Arrays.toString(depthByTaxon));

            //Filter taxa by depths and taxa in list.
            // This creates a list of the indices (indices, not depth, is stored).
            // Maximum size of entries for each tag is set to 10
            int[] taxaWithSufficientDepth=IntStream.range(0,depthByTaxon.length)
                    .filter(ti -> taxaList()==null || taxaList().contains(taxons.get(ti)))
                    .filter(ti -> depthByTaxon[ti]>minTaxaDepth())
                    .toArray();
            System.out.println(Arrays.toString(taxaWithSufficientDepth));

            for (Tag tag : tagTaxaMap.keySet()) {
                int[] depthTags=tagTaxaMap.get(tag).depths();
                double[] depthsForNewTL=IntStream.of(taxaWithSufficientDepth).mapToDouble(ti -> depthTags[ti]).toArray();
                double[] stdDepthsForNewTL=IntStream.of(taxaWithSufficientDepth).mapToDouble(ti -> (double)(depthTags[ti]+0.01)/(double)depthByTaxon[ti]).toArray();
                double sum= DoubleStream.of(depthsForNewTL).sum();
                if(sum<minTagCount()) continue;
                tagDepthMap.put(tag,stdDepthsForNewTL);
                tagQueueMap.put(tag,Queues.synchronizedQueue(MinMaxPriorityQueue.orderedBy(byCorrOrderingTag).maximumSize(10).create()));
            }

//            for (int i = 0; i < 20; i++) {
//                System.out.println(i);
//                Tag testTag=tagList.get(i);
//                System.out.println(Arrays.toString(tagTaxaMap.get(testTag).depths()));
//                System.out.println(Arrays.toString(tagDepthMap.get(testTag)));
//            }


            for (int t = 0; t < taxaWithSufficientDepth.length; t++) {
                System.out.printf("%s\t%d%n",taxons.get(taxaWithSufficientDepth[t]).getName(), depthByTaxon[taxaWithSufficientDepth[t]]);
            }


            // For each tag in the tagDepth Map:
            //  1.  loop through the depth map (for tag/tag correlation)
            //       - skip if entry.hashCode() >- entry2.hashCode().  this
            //          is ensuring just the lover triangle is populated
            //       = calculate the F2 Deviations using defined values in seExpectations map defined above
            //       - add the calcF2Deviations() value to tagQueueMap for both tag1 and tag2
            tagDepthMap.entrySet().stream().parallel().forEach(entry -> {
                for (Map.Entry<Tag, double[]> entry2 : tagDepthMap.entrySet()) {
                    if(entry.hashCode()>=entry2.hashCode()) continue;
                    double gF2 = calcF2Deviations(entry.getValue(),entry2.getValue(),segExpectations.get("RIL"),0.00001);
                    tagQueueMap.get(entry.getKey()).add(new Tuple<Double,Tag>(gF2,entry2.getKey()));
                    tagQueueMap.get(entry2.getKey()).add(new Tuple<Double,Tag>(gF2,entry.getKey()));
                }
            });

            System.out.println("Done");

            // For each entry in the tagQueueMap (populated above):
            //    1.  get all values for the tag
            //    2.  create an "evenness" array of doubles, store calculations from "evenRatios()" method
            //    3.  Map<Tuple<Tag,Tag>,Tuple<Double,String>> results  gets:
            //         - the tags as the key in the Tuple<tag,tag>,
            //         - the value is a tuple of <CalcF2Deviation double value for this pair, the evenness array as a String> 
            Map<Tuple<Tag,Tag>,Tuple<Double,String>> results=new HashMap<>();
            tagQueueMap.entrySet().stream().forEach(tagQueueEntry -> {
                tagQueueEntry.getValue().forEach(x2tagTuple -> {
                    double[] evenness=evenRatios(tagDepthMap.get(tagQueueEntry.getKey()),tagDepthMap.get(x2tagTuple.getY()), 0.00001);
                    results.put(new Tuple<>(tagQueueEntry.getKey(),x2tagTuple.getY()),
                            new Tuple<>(x2tagTuple.getX(),Arrays.toString(evenness)));
                });
            });

            System.out.println(repGenData.getTags().size());
 //           System.out.println(results.toString());
            repGenData.addMappingApproach(processName());
            repGenData.putTagTagStats(results, processName());

            
            ((RepGenSQLite)repGenData).close();
            
        } catch (Exception exc) {
            System.out.println("RepGenLDAnalysis:process_data:  processing error");
            exc.printStackTrace();
        }
        System.out.println("Process took " + (System.nanoTime() - totalTime)/1e9 + " seconds.\n");
        return null;
    }

    // THis method is similar to "odds ratio" calculation.  IT is used to show
    // association between tag1/tag2 as determined by their presence in the same taxa.
    private double[] evenRatios(double[] tag1cnt, double[] tag2cnt, double threshold) {
        double depth11=0, depth12of1=0, depth12of2=0, depth22=0;
        double count11=0, count12=0, count22=0;
        for (int i = 0; i < tag1cnt.length; i++) {
            if ((tag1cnt[i]>threshold)&&(tag2cnt[i]>threshold)) {
                depth12of1+=tag1cnt[i];
                depth12of2+=tag2cnt[i];
                count12++;
            } else if((tag1cnt[i]>threshold)&&(tag2cnt[i]<threshold)){
                depth11+=tag1cnt[i];
                count11++;
            } else if((tag1cnt[i]<threshold)&&(tag2cnt[i]>threshold)){
                depth22+=tag2cnt[i];
                count22++;
            }
        }
        return new double[]{depth11/(2*count11),depth12of1/count12, depth12of2/count12, depth22/(2*count11)};
    }

    private double calcF2Deviations(double[] tag1cnt, double[] tag2cnt, double[][] expFreq, double threshold) {
        double[][] cnts=new double[2][2];
        for (int i = 0; i < tag1cnt.length; i++) {
            int tag1present=(tag1cnt[i]>threshold)?1:0;
            int tag2present=(tag2cnt[i]>threshold)?1:0;
            cnts[tag1present][tag2present]+=1.0;
        }
        double xchi=0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if(i==0 && j==0) continue;  //consider whether to test the zero class
                double exp=tag1cnt.length*expFreq[i][j];
                xchi+=(cnts[i][j]-exp)*(cnts[i][j]-exp)/exp;
            }
        }
//        if(Math.abs(xchi)<3) {
//            System.out.println(Arrays.toString(tag1cnt));
//            System.out.println(Arrays.toString(tag2cnt));
//            System.out.println(xchi+":"+Arrays.deepToString(cnts));
//        }
        return xchi;
    }


    private Optional<double[][]> filterCombinedDepth(int minimumCombinedDepth, int minObservations, double[] a1, double[] a2) {
        double[][] result=new double[2][];
        result[0]=IntStream.range(0,a1.length).filter(i -> a1[i]+a2[i]>minimumCombinedDepth).mapToDouble(i -> a1[i]).toArray();
        result[1]=IntStream.range(0,a1.length).filter(i -> a1[i]+a2[i]>minimumCombinedDepth).mapToDouble(i -> a2[i]).toArray();
        if(result[0].length<minObservations) return Optional.empty();
        return Optional.of(result);
    }



    
    public void calculateCorrelations(Multimap<Tag,TagCorrelationInfo> tagTagCorrelations, Map <Tag, TaxaDistribution> tagTaxaMap, 
            Tag tag1, Tag tag2, double[] ddepths1,double[] depthsPrime1) {
        TaxaDistribution tag2TD = tagTaxaMap.get(tag2);
        if (tag2TD == null) {                       
            System.out.println("GetTagTaxaDist: got null tagTD for sequence " + tag2.sequence());
            return ;
        }
        
        // I need doubles below !!
        int[] depths2 = tag2TD.depths(); // gives us the depths for each taxon
        double[] ddepths2 = new double[depths2.length];
        // Apparently no shorter method of casting int to double 
        for (int idx = 0; idx < depths2.length; idx++) {
            ddepths2[idx] = (double)depths2[idx];
        }
        double[] depthsPrime2 = new double[depths2.length];
        for (int idx = 0; idx < depthsPrime1.length; idx++) {
            // boolean - presence/absence
            if (ddepths2[idx] > 0) depthsPrime2[idx] = 1;
            else depthsPrime2[idx] = 0;
        }
        // From analysis.association.GenomicSelectionPlugin.java:
        PearsonsCorrelation Pearsons = new PearsonsCorrelation();
        double p1 = Pearsons.correlation(ddepths1,ddepths2);
        
        SpearmansCorrelation Spearmans = new SpearmansCorrelation();
        double spearval = Spearmans.correlation(ddepths1, ddepths2);
        
        double p2 = Pearsons.correlation(depthsPrime1, depthsPrime2);
        
        // Count number of times both tags appeared in a taxa, number
        // of times neither tag appeared in a taxa, number of times
        // tag1 appeared but not tag2, and number of times tag2 appeared by not tag1
        int t1Nott2 = 0;
        int t2Nott1 = 0;
        int neither = 0;
        int both = 0;
        
        for (int didx = 0; didx < depthsPrime2.length; didx++) {
            if (depthsPrime1[didx] > 0 && depthsPrime2[didx] > 0) both++;
            if (depthsPrime1[didx] == 0 && depthsPrime2[didx] == 0) neither++;
            if (depthsPrime1[didx] > 0 && depthsPrime2[didx] == 0) t1Nott2++;
            if (depthsPrime1[didx] == 0 && depthsPrime2[didx] > 0) t2Nott1++;
        }
        // Calculate r-squared based on presence/absence of tags at each taxa.
        double r2 = LinkageDisequilibrium.calculateRSqr(neither, t1Nott2, t2Nott1, both, minTaxa());
        TagCorrelationInfo tci = new TagCorrelationInfo(tag2,p1,spearval,p2,r2);
        tagTagCorrelations.put(tag1, tci);       
    }
    
    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getToolTipText() {
        // TODO Auto-generated method stub
        return null;
    }


    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
     public static void main(String[] args) {
         GeneratePluginCode.generate(RepGenLDAnalysisPlugin.class);
     }

    /**
     * Convenience method to run plugin with one return object.
     */
//    // TODO: Replace <Type> with specific type.
//    public <Type> runPlugin(DataSet input) {
//        return (<Type>) performFunction(input).getData(0).getData();
//    }

    /**
     * Input database file with tags and taxa distribution
     *
     * @return Input DB
     */
    public String inputDB() {
        return myInputDB.value();
    }

    /**
     * Set Input DB. Input database file with tags and taxa
     * distribution
     *
     * @param value Input DB
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin inputDB(String value) {
        myInputDB = new PluginParameter<>(myInputDB, value);
        return this;
    }

    /**
     * Minimum number of taxa that must be present for R-squared
     * to be calculated.
     *
     * @return Min Taxa for RSquared
     */
    public Integer minTaxa() {
        return minTaxa.value();
    }

    /**
     * Set Min Taxa for RSquared. Minimum number of taxa that
     * must be present for R-squared to be calculated.
     *
     * @param value Min Taxa for RSquared
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin minTaxa(Integer value) {
        minTaxa = new PluginParameter<>(minTaxa, value);
        return this;
    }

    /**
     * Minimum number of taxa that must be present for R-squared
     * to be calculated.
     *
     * @return Min Taxa for RSquared
     */
    public Integer minTaxaDepth() {
        return minTaxaDepth.value();
    }

    /**
     * Set Min Taxa for RSquared. Minimum number of taxa that
     * must be present for R-squared to be calculated.
     *
     * @param value Min Taxa for RSquared
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin minTaxaDepth(Integer value) {
        minTaxaDepth = new PluginParameter<>(minTaxaDepth, value);
        return this;
    }

    /**
     * Minimum number of taxa that must be present for R-squared
     * to be calculated.
     *
     * @return Taxa list to test
     */
    public TaxaList taxaList() {
        return myTaxaList.value();
    }

    /**
     * Set Taxa list to test. Minimum number of taxa that
     * must be present for R-squared to be calculated.
     *
     * @param value Taxa list to test
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin taxaList(TaxaList value) {
        myTaxaList = new PluginParameter<>(myTaxaList, value);
        return this;
    }

    /**
     * Minimum number of taxa that must be present for R-squared
     * to be calculated.
     *
     * @return Min tag count after filtering
     */
    public Integer minTagCount() {
        return minTagCount.value();
    }

    /**
     * Set Min tag count after filtering. Minimum number of
     * taxa that must be present for R-squared to be calculated.
     *
     * @param value Min tag count after filtering
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin minTagCount(Integer value) {
        minTagCount = new PluginParameter<>(minTagCount, value);
        return this;
    }

    /**
     * Name of the process
     *
     * @return Name for process
     */
    public String processName() {
        return myProcessName.value();
    }

    /**
     * Set Name for process. Name of the process
     *
     * @param value Name for process
     *
     * @return this plugin
     */
    public RepGenLDAnalysisPlugin processName(String value) {
        myProcessName = new PluginParameter<>(myProcessName, value);
        return this;
    }
}
