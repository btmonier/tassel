
/**
 * @author lcj34
 *
 */

package net.maizegenetics.dna.tag;

import com.google.common.collect.*;
import com.google.common.io.CharStreams;
import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;
import org.sqlite.SQLiteConfig;

import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.*;
import java.util.*;


public class RepGenSQLite implements RepGenDataWriter, AutoCloseable {
    private Connection connection = null;

    /*These maps contain  objects that are most queried by users.  This is a not the simplest way to do this, which
        would probably be done cleaner just with queries against the databases.  However, there are large performance
        gains in the case of SQLite (or at least with my ability to optimize).

        The logic behind this most of the datasets are relatively small (except tagTagIDMap), and this prevents creation
        of these objects over and over again.
     */
    private BiMap<Tag,Integer> tagTagIDMap; // Tag is AbstractTag.java
    
    // RefTagData is the key to BiMap refTagRefTagIDMAP because:
    // The reference tag sequence can show up in different places.
    // refTag is unique based on sequence, chromosome, seqlen, position, refGenomeID.  
    // Sequence can be duplicated on same chrom and on multiple chroms.
    private Map<String,Integer> mappingApproachToIDMap;
    private BiMap<String,Integer> referenceGenomeToIDMap;
    private SortedMap<Position,Integer> physicalMapPositionToIDMap;


    private TaxaList myTaxaList;

    PreparedStatement tagTaxaDistPS;
    PreparedStatement posTagMappingInsertPS;
    PreparedStatement taxaDistWhereTagMappingIDPS;
    PreparedStatement tagtagAlignmentInsertPS;
    PreparedStatement tagAlignForNonRefTagPS;
    PreparedStatement tagTagCorrelationInsertPS;
    PreparedStatement tagCorrelationsForTag1PS;
    PreparedStatement tagCorrelationsForTag2PS;

    public RepGenSQLite(String filename) {
        try{
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        // create a database connection

        try {
            boolean doesDBExist= Files.exists(Paths.get(filename));

            SQLiteConfig config=new SQLiteConfig();

            connection = DriverManager.getConnection("jdbc:sqlite:"+filename,config.toProperties());
            connection.setAutoCommit(true);  //This has massive performance effects
            Statement statement = connection.createStatement();
            statement.setQueryTimeout(30);  // set timeout to 30 sec.
            //                System.out.println(schema);
            if(doesDBExist==false) {
                String schema = CharStreams.toString(new InputStreamReader(RepGenSQLite.class.getResourceAsStream("RepGenSchema.sql")));
                statement.executeUpdate(schema);
            }
            initPreparedStatements();
            loadReferenceGenomeHash();
            loadTagHash();
            loadMappingApproachHash();
            loadTaxaList();
        }
        catch(Exception e)
        {
            // if the error message is "out of memory",
            // it probably means no database file is found
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }

    @Override
    public void close() throws Exception {
        System.out.println("Closing SQLDB");
        connection.close();
    }

    private void initPreparedStatements() {
        try{
            posTagMappingInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into TagMapping (tagid, posid, method_id, bp_error, cm_error)" +
                    " values(?,?,?,?,?)");
            
            tagTaxaDistPS=connection.prepareStatement("select depthsRLE from tagtaxadistribution where tagid=?");

            tagtagAlignmentInsertPS=connection.prepareStatement(
                    "INSERT into tag_tag_Alignments (tag1id, tag2id, score )" +
                    " values(?,?,?)");
            tagTagCorrelationInsertPS=connection.prepareStatement(
                    "INSERT into tagCorrelations (tag1id, tag2id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 )" +
                    " values(?,?,?,?,?,?)");
            // because there can be a tagID X in both tag and refTag table, you must
            // specify that this query only wants the values where tag1 is NOT a ref
            // (IE tag1ID comes from the tagTagIDMap)
            // THis gets both the non-ref and ref tag alignments for a particular non-ref tag.
            // NOTE: SQLite has no real boolean field.  The values are stored as 0 and 1
            tagAlignForNonRefTagPS= connection.prepareStatement(
                    "select tag2id,  score " +
                    "from tag_tag_Alignments where tag1id=? and score >= ?");
            tagCorrelationsForTag1PS = connection.prepareStatement(
                    "select tag2id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 " +
                    "from tagCorrelations where tag1id=?");
            tagCorrelationsForTag2PS = connection.prepareStatement(
                    "select tag1id, t1t2_pearson, t1t2_spearman, pres_abs_pearson, r2 " +
                    "from tagCorrelations where tag2id=?");
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTagHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from tag");
            int size=rs.getInt(1);
            System.out.println("size of all tags in tag table=" + size);
            if(tagTagIDMap==null || size/(tagTagIDMap.size()+1)>3) tagTagIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from tag");
            boolean hasName;
            try{rs.findColumn("tagName");hasName=true;}catch (SQLException e) {hasName=false;}
            while(rs.next()) {
                TagBuilder tagBuilder=TagBuilder.instance(rs.getBytes("sequence"),rs.getShort("seqlen"));
                if(hasName) tagBuilder.name(rs.getString("tagName"));
                tagTagIDMap.putIfAbsent(tagBuilder.build(),rs.getInt("tagid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }


    private void loadPhysicalMapPositionHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from physicalmapposition");
            int size=rs.getInt(1);
            System.out.println("Before loading new positions, size of all positions in physicalMapPosiiton table="+size);
            if(physicalMapPositionToIDMap==null) {physicalMapPositionToIDMap=new TreeMap<>();}
            else if(size==physicalMapPositionToIDMap.size()) return;
            rs=connection.createStatement().executeQuery("select * from physicalMapPosition");
            while(rs.next()) {
                Position p=new GeneralPosition
                        .Builder(new Chromosome(rs.getString("chromosome")),rs.getInt("physical_position"))
                        .strand(rs.getByte("strand"))
                        .build();
                physicalMapPositionToIDMap.putIfAbsent(p, rs.getInt("posid"));
            }
            rs=connection.createStatement().executeQuery("select count(*) from physicalmapposition");
            size=rs.getInt(1);
            System.out.println("After loading new positions, size of all positions in physicalMapPosiiton table="
            +size + ", size of physicalMapPositionToIDMAP: " + physicalMapPositionToIDMap.size());
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadReferenceGenomeHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from referenceGenome");
            int size=rs.getInt(1);
            System.out.println("size of all references in referenceGenome table="+size);
            if(size==0) {
                connection.createStatement().executeUpdate("insert into referenceGenome (refname) " +
                        "values('unknown')");
                size=1;
            }
            referenceGenomeToIDMap=HashBiMap.create(size);
            rs=connection.createStatement().executeQuery("select * from referenceGenome");
            while(rs.next()) {
                referenceGenomeToIDMap.put(rs.getString("refname"), rs.getInt("refid"));
                System.out.println("refence name from referenceGenome: " + rs.getString("refname"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    private void loadMappingApproachHash() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from mappingApproach");
            int size=rs.getInt(1);
            System.out.println("size of all approaches in mappingApproach table="+size);
            if(size==0) {
                connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                        "values('unknown','unknown','unknown')");
                size=1;
            }
            mappingApproachToIDMap=new HashMap<>(size);
            rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                mappingApproachToIDMap.put(rs.getString("approach"), rs.getInt("mapappid"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private void loadTaxaList() {
        try{
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from taxa");
            int size=rs.getInt(1);
            System.out.println("size of all taxa in taxa table="+size);
            TaxaListBuilder tlb=new TaxaListBuilder();
            rs=connection.createStatement().executeQuery("select * from taxa");
            while(rs.next()) {
                tlb.add(new Taxon(rs.getString("name")));
            }
            myTaxaList=tlb.build();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public TaxaList getTaxaList() {
        if(myTaxaList==null) loadTaxaList();
        return myTaxaList;
    }

    @Override
    public Map<Tag, String> getTagsNameMap() {
        try{
            Map<Tag, String> tagNameMap=new HashMap<>(tagTagIDMap.size()+1);
            ResultSet rs=connection.createStatement().executeQuery("select * from tag");
            while(rs.next()) {
                tagNameMap.put(TagBuilder.instance(rs.getBytes("sequence"), rs.getShort("seqlen")).build(), rs.getString("tagName"));
            }
            return tagNameMap;
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }

    // Only when called from RepGenLoadSeqToDBPLugin do we have the number of instances
    // and the quality score.  Other places (from within this file) don't have that data.
    // Populate table based on information sent in.
    @Override
    public boolean putAllTag(Set<Tag>tags,Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS) {
    //public boolean putAllTag(Set<Tag> tags) {
        int batchCount=0, totalCount=0;
        if (tagInstanceAverageQS != null) {
            try {
                connection.setAutoCommit(false);
                PreparedStatement tagInsertPS=
                        connection.prepareStatement("insert into tag (sequence, seqlen,isReference,qualityScore,numTagInstances,sequencetext) values(?,?,?,?,?,?)");
                for (Map.Entry<Tag, Tuple<Integer,String>> entry : tagInstanceAverageQS.entrySet()) {
                    Tag tag = entry.getKey();
                    if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                    int numInstances = entry.getValue().x;
                    String qscore = entry.getValue().y;
                    tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                    tagInsertPS.setShort(2, tag.seqLength());
                    tagInsertPS.setBoolean(3, tag.isReference());
                    tagInsertPS.setString(4, qscore);
                    tagInsertPS.setInt(5, numInstances);
                    tagInsertPS.setString(6, tag.sequence());
                    tagInsertPS.addBatch();
                    batchCount++;
                    totalCount++;
                    if(batchCount>10000) {
                       // System.out.println("tagInsertPS.executeBatch() "+batchCount);
                        tagInsertPS.executeBatch();
                        //connection.commit();
                        batchCount=0;
                    }
                }
                tagInsertPS.executeBatch();
                connection.setAutoCommit(true);
            } catch (SQLException e) {
                e.printStackTrace();
                return false;
            }
        } else {
            try {
                connection.setAutoCommit(false);
                PreparedStatement tagInsertPS=
                        connection.prepareStatement("insert into tag (sequence, seqlen,isReference,sequencetext) values(?,?,?,?)");
                for (Tag tag: tags) {                
                    if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip                    
                    tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                    tagInsertPS.setShort(2, tag.seqLength());
                    tagInsertPS.setBoolean(3, tag.isReference());
                    tagInsertPS.setString(4, tag.sequence());
                    tagInsertPS.addBatch();
                    batchCount++;
                    totalCount++;
                    if(batchCount>100000) {
                       // System.out.println("tagInsertPS.executeBatch() "+batchCount);
                        tagInsertPS.executeBatch();
                        //connection.commit();
                        batchCount=0;
                    }
                }
                tagInsertPS.executeBatch();
                connection.setAutoCommit(true);
            } catch (SQLException e) {
                e.printStackTrace();
                return false;
            }
        }

        if(totalCount>0) {
            System.out.println("RepGenSQLite:putAllTag, totalCount=" + totalCount + ",loadingHash");
            loadTagHash();
        } 
        return true;
    }


    @Override
    public boolean putAllNamesTag(Map<Tag, String> tagNameMap) {
        int batchCount=0, totalCount=0;
        try {
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tag (sequence, seqlen, tagName) values(?,?,?)");
            for (Map.Entry<Tag, String> entry : tagNameMap.entrySet()) {
                Tag tag=entry.getKey();
                if(tagTagIDMap.containsKey(tag)) continue;  //it is already in the DB skip
                tagInsertPS.setBytes(1, tag.seq2BitAsBytes());
                tagInsertPS.setShort(2, tag.seqLength());
                tagInsertPS.setString(3, entry.getValue());
                tagInsertPS.addBatch();
                batchCount++;
                totalCount++;
                if(batchCount>100000) {
                    System.out.println("tagInsertPS.executeBatch() "+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
            return false;
        }
        if(totalCount>0) loadTagHash();
        return true;
    }

    @Override
    public void putTaxaList(TaxaList taxaList) {
        try {
            connection.createStatement().execute("delete from taxa");
            connection.setAutoCommit(false);
            PreparedStatement taxaInsertPS=connection.prepareStatement("insert into taxa (taxonid, name) values(?,?)");
            for (int i = 0; i < taxaList.size(); i++) {
                taxaInsertPS.setInt(1, i);
                taxaInsertPS.setString(2, taxaList.get(i).getName());
                taxaInsertPS.addBatch();

            }
            taxaInsertPS.executeBatch();
            connection.setAutoCommit(true);
            loadTaxaList();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }


    @Override
    public void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap) {
        int batchCount=0;
        try {
            int numTaxa=myTaxaList.numberOfTaxa();
            connection.setAutoCommit(false);
            PreparedStatement tagInsertPS=connection.prepareStatement("insert into tagtaxadistribution (tagid, depthsRLE, totalDepth) values(?,?,?)");
            for (Map.Entry<Tag, TaxaDistribution> entry : tagTaxaDistributionMap.entrySet()) {
                int tagID=tagTagIDMap.get(entry.getKey());
                tagInsertPS.setInt(1,tagID);
                if(entry.getValue().maxTaxa()!=numTaxa) throw new IllegalStateException("Number of taxa does not agree with taxa distribution");
                tagInsertPS.setBytes(2, entry.getValue().encodeTaxaDepth());
                tagInsertPS.setInt(3, entry.getValue().totalDepth());
                tagInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                    System.out.println("putTaxaDistribution next"+batchCount);
                    tagInsertPS.executeBatch();
                    //connection.commit();
                    batchCount=0;
                }
            }
            tagInsertPS.executeBatch();
            connection.setAutoCommit(true);  
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    //This method is called from RepGenAlignerPlugin to add  tags
    // to the  tagAlignments table.
    @Override
    public void putTagTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap) {
        int batchCount=0;
        try {
            connection.setAutoCommit(false);
            
            for (Map.Entry<Tag, AlignmentInfo> entry : tagAlignInfoMap.entries()) {
                // Put tag alignments into the tagAlignments table
                AlignmentInfo ai=entry.getValue();
                int ind=1;
 
                tagtagAlignmentInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                tagtagAlignmentInsertPS.setInt(ind++, tagTagIDMap.get(ai.tag2()));
                tagtagAlignmentInsertPS.setInt(ind++, ai.score());  // alignment score

                tagtagAlignmentInsertPS.addBatch();
                batchCount++;
                if(batchCount>100000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    tagtagAlignmentInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            tagtagAlignmentInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numAlignments from tag_tag_Alignments");
            if (rs.next()) {
                System.out.println("Total alignments in tag_tag_Alignments table: " + rs.getInt("numAlignments"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }


    //THis method is called from RepGenAlignerPlugin to add reference kmers
    // to the db.  This method adds data to the tag, physicalMapPosition, and
    // tagMapping tables.
    //KEEP THIS ED
    @Override
    public void putTagPositionMapping(Multimap<Tag, Position> refTagPositionMap, String refGenome, String mappingMethod) {
        int batchCount=0;
        loadReferenceGenomeHash();
        try {
            //putAllRefTag(refTagPositionMap,refGenome);
            putAllTag(refTagPositionMap.keySet(),null);
            System.out.println("putREfTagMaping: size of map: " + refTagPositionMap.size()
              + ", keyset size: " + refTagPositionMap.keySet().size() + ", values size: " 
                    + refTagPositionMap.values().size());
            putPhysicalMapPositionsIfAbsent(refTagPositionMap.values(),refGenome);
            int mappingApproachID=mappingApproachToIDMap.get(mappingMethod);
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, Position> entry : refTagPositionMap.entries()) {
                // add reference tags to tagMapping table - map refTag to physicalMapPosition
                Position entrypos=entry.getValue();
                int ind=1;
               // Defined above:  posTagInsertPS=connection.prepareStatement(
                //"INSERT OR IGNORE into TagMapping (tagid, posid, method_id, bp_error, cm_error)" +
                // " values(?,?,?,?,?)");
                posTagMappingInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey())); // refTagID
                posTagMappingInsertPS.setInt(ind++, physicalMapPositionToIDMap.get(entrypos)); // position
                posTagMappingInsertPS.setInt(ind++, mappingApproachID); // method_id
                posTagMappingInsertPS.setInt(ind++, 0);  //todo this needs to be input data (bp_error)
                posTagMappingInsertPS.setFloat(ind++, 0);  //todo this needs to be input data (cm_error)

                posTagMappingInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) { // LCJ - changed from 100000 - writes are REALLY slow when db is big, 79921 tags
                   // System.out.println("putTagAlignments next"+batchCount);
                    posTagMappingInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posTagMappingInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (DISTINCT physical_position) as numPhysicalSites from physicalMapPosition");
            if (rs.next()) {
                System.out.println("Total number of distinct physical position sites: " + rs.getInt("numPhysicalSites"));
            }
            rs = connection.createStatement().executeQuery("select count (*) as numPhysicalSites from physicalMapPosition");
            if (rs.next()) {
                System.out.println("Total number of physical position sites: " + rs.getInt("numPhysicalSites"));
            }
            PreparedStatement physMapNumFromTCPMP = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,physical_position from physicalMapPosition " +
                    "GROUP BY physical_position) where tgcnt=?");
            physMapNumFromTCPMP.setInt(1, 1);// having 1 tag
            rs = physMapNumFromTCPMP.executeQuery();

            if (rs.next()) {
                System.out.println("Number of physical position sites with 1 tag: " + rs.getInt("numSites"));
            }
            physMapNumFromTCPMP.setInt(1, 2);// having 2 tag
            rs = physMapNumFromTCPMP.executeQuery();
            if (rs.next()) {
                System.out.println("Number of physical position sites with 2 tags: " + rs.getInt("numSites"));
            }
            physMapNumFromTCPMP.setInt(1, 3);// having 3 tags
            rs = physMapNumFromTCPMP.executeQuery();
            if (rs.next()) {
                System.out.println("Number of physical position sites with 3 tags: " + rs.getInt("numSites"));
            }

            PreparedStatement cutSiteGreaterThanPS = connection.prepareStatement(
                    "select count(*) as numSites from (select count(*) as tgcnt,physical_position from physicalMapPosition " +
                    "GROUP BY physical_position) where tgcnt>?");
            cutSiteGreaterThanPS.setInt(1, 3);// having > 3 tags
            rs = cutSiteGreaterThanPS.executeQuery();
            if (rs.next()) {
                System.out.println("Number of cut sites with more than 3 tags: " + rs.getInt("numSites"));
            }           
        } catch (SQLException e) {
            e.printStackTrace();
        }

    }

    private int getMappingApproachID(Position p) throws SQLException{
        String mapApp=p.getAnnotation().getTextAnnotation("mappingapproach")[0];
        if(mapApp==null) return mappingApproachToIDMap.get("unknown");
        Integer val=mappingApproachToIDMap.get(mapApp);
        if(val==null) {
            connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                    "values('"+mapApp+"','unknown','unknown')");
            loadMappingApproachHash();
            return mappingApproachToIDMap.get(mapApp);
        } else return val;
    }

    private int getReferenceGenomeID(String refGenome) throws SQLException{                
        Integer val=referenceGenomeToIDMap.get(refGenome);
        if(val==null) {
            connection.createStatement().executeUpdate("insert into referenceGenome (refname) " +
                    "values('"+refGenome+"')");
            loadReferenceGenomeHash();
            return referenceGenomeToIDMap.get(refGenome);
        } else return val;
    }
    
    @Override
    public void addReferenceGenome(String name) {
        Integer val=mappingApproachToIDMap.get(name);
        if(val==null) {
            try {
                connection.createStatement().executeUpdate("insert or ignore into referenceGenome (refname) " +
                        "values('"+name+"')");
            } catch (SQLException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            loadReferenceGenomeHash();
        } 
    }

    @Override
    public void addMappingApproach(String name) {
        Integer val=mappingApproachToIDMap.get(name);
        if(val==null) {
            try {
                connection.createStatement().executeUpdate("insert into mappingApproach (approach, software, protocols) " +
                        "values('"+name+"','unknown','unknown')");
            } catch (SQLException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            loadMappingApproachHash();
        } 
    }
    
    @Override
    public void setTagAlignmentBest(Tag tag, Position position, boolean isBest) {

    }

    @Override
    public boolean putTagAlignmentApproach(String tagAlignmentName, String protocol) {
        return false;
    }

    @Override
    public TaxaDistribution getTaxaDistribution(Tag tag) {
        int tagid=tagTagIDMap.get(tag);
        try {
            tagTaxaDistPS.setInt(1,tagid);
            ResultSet rs=tagTaxaDistPS.executeQuery();
            while (rs.next()) {
                // only non-ref tags have taxa distribution
                return TaxaDistBuilder.create(rs.getBytes(1));
            }
            return null; // no taxa dist for this tag
            
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public Set<Tag> getTags() {
        return tagTagIDMap.keySet();
    }
    
    @Override
    public Set<RefTagData> getRefTags() {
        return null;
    }
    
    @Override
    public PositionList getPhysicalMapPositions() {
        if(physicalMapPositionToIDMap == null) loadPhysicalMapPositionHash();
        PositionListBuilder plb=new PositionListBuilder();
        physicalMapPositionToIDMap.keySet().stream()
                .forEach(p -> plb.add(p));
        plb.sortPositions();
        return plb.build();
    }

    @Override
    public PositionList getPhysicalMapPositions(Chromosome chromosome, int firstPosition, int lastPosition) {
        PositionListBuilder plb=new PositionListBuilder();
        plb.addAll(getPositionSubMap(chromosome,firstPosition,lastPosition).keySet());
        return plb.build();
    }

    @Override
    public PositionList getSNPPositionsForChromosomes(Integer startChr, Integer endChr) {
        return null;
    }

    @Override
    public Set<Tag> getTagsForTaxon(Taxon taxon) {
        ImmutableSet.Builder<Tag> tagBuilder=new ImmutableSet.Builder<>();
        int taxonIndex=myTaxaList.indexOf(taxon);
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from tagtaxadistribution");
            while(rs.next()) {
                if(TaxaDistBuilder.create(rs.getBytes("depthsRLE")).depths()[taxonIndex]>0) {
                    tagBuilder.add(tagTagIDMap.inverse().get(rs.getInt("tagid")));
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }

    @Override
    public Map<Tag, Integer> getTagDepth(Taxon taxon, Position position) {
        return null;
    }

    @Override
    public Map<Tag, Integer> getTagsWithDepth(int minimumDepth) {
        ImmutableMap.Builder<Tag, Integer> tagBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery(
                    "select tagid, totalDepth from tagtaxadistribution where totalDepth >= "+minimumDepth);
            while(rs.next()) {
                tagBuilder.put(tagTagIDMap.inverse().get(rs.getInt(1)),rs.getInt(2));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tagBuilder.build();
    }



    private Map<Position,Integer> getPositionSubMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        if(physicalMapPositionToIDMap==null) loadPhysicalMapPositionHash();
        Position startPos=new GeneralPosition.Builder(chromosome,firstPosition).build();
        if(lastPosition<0) lastPosition=Integer.MAX_VALUE;
        Position lastPos=new GeneralPosition.Builder(chromosome,lastPosition).build();
        return physicalMapPositionToIDMap.subMap(startPos,lastPos);
    }

    @Override
    public Map<String, String> getTagAlignmentApproaches() {
        ImmutableMap.Builder<String,String> appBuilder=new ImmutableMap.Builder<>();
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from mappingApproach");
            while(rs.next()) {
                appBuilder.put(rs.getString("approach"), rs.getString("software") + ":" + rs.getString("approach"));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return appBuilder.build();
    }

    @Override
    public Map<Position, Map<Tag, Tuple<Boolean, TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition) {
        return null;
    }

    @Override
    public Map<Position, Map<Tag, TaxaDistribution>> getCutPosForStrandTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition, boolean strand) {
        return null;
    }

    private void putPhysicalMapPositionsIfAbsent(Collection<Position> positions, String refGenome) {
        try {
            int batchCount=0;
            int positionExists = 0;
            int newPosition = 0;
            if(physicalMapPositionToIDMap==null) loadPhysicalMapPositionHash();
            connection.setAutoCommit(false);
            System.out.println("putPhysicalMapPositionsIfAbsent: size of positions: " + positions.size());
            PreparedStatement posInsertPS=connection.prepareStatement(
                    "INSERT OR IGNORE into physicalMapPosition (reference_genome_id, chromosome, physical_position, strand) values(?,?,?,?)");
            for (Position p : positions) {
                if(physicalMapPositionToIDMap.containsKey(p)) {
                    positionExists++;
                    continue;
                }
                newPosition++;
                posInsertPS.setInt(1, getReferenceGenomeID(refGenome)); 
                posInsertPS.setString(2, p.getChromosome().getName());
                posInsertPS.setInt(3, p.getPosition());
                posInsertPS.setByte(4, p.getStrand());
                posInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putPhysicalMapPositionsIfAbsent next"+batchCount);
                    posInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            posInsertPS.executeBatch();
            if(batchCount>0) loadPhysicalMapPositionHash();
            connection.setAutoCommit(true);
            System.out.println("putPhysicalMapPositionsIfAbsent: end, positionExists: " + positionExists
                    + ", newPositions: " + newPosition);
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public List<Chromosome> getChromosomesFromCutPositions() {
        // Get a list of chromosomes from the cutPosition table
        List<Chromosome> chromList = new ArrayList<>();
        try {
            ResultSet rs = connection.createStatement().executeQuery("select DISTINCT chromosome from cutPosition");
            while(rs.next()) {
                Chromosome chrom=new Chromosome(rs.getString("chromosome"));
                chromList.add(chrom);
            }
        } catch (SQLException exc) {
            exc.printStackTrace();
            return null;
        }
        return chromList;
    }

    @Override
    public void clearTagTaxaDistributionData() {
        // Clear tags and tagtaxadist tables populated from GBSSeqToTagDBPlugin
        // TODO add refTag clearing when this is revisited fro REpGEn
        try {
            boolean rs = connection.createStatement().execute("delete FROM tagtaxadistribution"); 
            rs = connection.createStatement().execute("delete FROM tag");
            tagTagIDMap = null;
            loadTagHash();
            
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting tagtaxadistribution data");
            exc.printStackTrace();
        }              
    }

    @Override
    public void clearAlignmentData() {
        // Clear tables populated from SAMToGBSdbPluging
        try {           
            boolean rs = connection.createStatement().execute("delete FROM tagCutPosition");
            rs = connection.createStatement().execute("delete FROM cutPosition");
            rs = connection.createStatement().execute("delete FROM mappingApproach");
            physicalMapPositionToIDMap = null;
            mappingApproachToIDMap = null;
            loadMappingApproachHash(); // this adds "unknown" to mappingApproachToIDMap
            loadReferenceGenomeHash();
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting alignment data");
            exc.printStackTrace();
        }              
    }

    @Override
    public void clearDiscoveryData() {
        // Clear all entries from tables populated from the DiscoverySNPCallerPluginV2 
        // The "delete" removes data, but keeps the table size.  The "vacuum" command
        // is NOT called as it rebuilds the entire data base which can be quite time intensive.
        // The rows needed for this table will be needed again in the subsequent run.  Vacuum also creates
        // a new file, so disk space requirements could double while vacuuming.
        try {         
            boolean rs = connection.createStatement().execute("delete FROM tagallele");            
            rs = connection.createStatement().execute("delete FROM snpposition"); 
            rs = connection.createStatement().execute("delete FROM allele");
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting discovery data");
            exc.printStackTrace();
        }       
    }

    @Override
    public void clearSNPQualityData() {
        // Clear table populated via SNPQualityProfilerPlugin
        try {
            boolean rs = connection.createStatement().execute("delete FROM snpQuality");         
        } catch (SQLException exc) {
            System.out.println("ERROR - problem deleting snpQuality data");
            exc.printStackTrace();
        }       
    }

    // This method is called  when a user wishes to append data to an existing db from
    // RegGenLoadSeqToDBPlugin;  also called from RepGenLDAnalysis for use in matrix
    // correlation calculations.
    // The map created must NOT be a fixed taxaDist map.  We intend to increment the values
    @Override
    public Map<Tag, TaxaDistribution> getAllTagsTaxaMap() {
        ImmutableMap.Builder<Tag,TaxaDistribution> tagTDBuilder = ImmutableMap.builder();
        loadTagHash(); // get updated tag map
        try {
            ResultSet rs = connection.createStatement().executeQuery("select tagid, depthsRLE from tagtaxadistribution");
            while(rs.next()) {
                Tag myTag = tagTagIDMap.inverse().get(rs.getInt("tagid"));
                TaxaDistribution myTD = TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                TaxaDistribution myTDIncr = TaxaDistBuilder.create(myTD); // convert to incrementable version
                tagTDBuilder.put(myTag, myTDIncr);
            }
            return tagTDBuilder.build();
        } catch (SQLException exc) {
            System.out.println("getAllTaxaMap: caught SQLException attempting to grab taxa Distribution ");
            exc.printStackTrace();
        }
        return tagTDBuilder.build();
    }

    @Override
    public Multimap<Tag, AlignmentInfo> getTagAlignmentsForTags(List<Tag> tags, int minscore) {
        // This method gets all tag-tag alignments for a list of tags (all data is for non-reference tags)
        ImmutableMultimap.Builder<Tag,AlignmentInfo> tagAIBuilder = ImmutableMultimap.builder();
        loadTagHash(); // get updated tag map
 
        try {
            for (Tag tag: tags){
                // For each tag on the list, get all its alignments
                Integer tagID = tagTagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    // should this return a NULL to indicate an error with the tag?
                    System.out.println("getAlignmentsForTag: no tagID in alignments table for tag: " + tag.sequence());
                    continue;
                }
                // The query specifies that tag1 is not a reference tag
                tagAlignForNonRefTagPS.setInt(1,tagID);
                tagAlignForNonRefTagPS.setInt(2, minscore);
                ResultSet rs = tagAlignForNonRefTagPS.executeQuery();
                while (rs.next()) {
                    int score = rs.getInt("score");
                    Tag tag2 = tagTagIDMap.inverse().get(rs.getInt("tag2id"));
                    AlignmentInfo ai = new AlignmentInfo(tag2,null,-1,-1, -1, null,score);
                    tagAIBuilder.put(tag,ai);
                }              
            }
        } catch (SQLException exc) {
            System.out.println("getTagAlignmentsForTag: caught SQLException attempting to grab alignment data ");
            exc.printStackTrace();
        }
        return tagAIBuilder.build();       
    }

    @Override
    public void  putTagTagCorrelationMatrix(Multimap<Tag,TagCorrelationInfo> tagCorrelationMap){
        int batchCount=0;
        loadTagHash(); // get updated list of tags
        try {
            connection.setAutoCommit(false);
            for (Map.Entry<Tag, TagCorrelationInfo> entry : tagCorrelationMap.entries()) {
                // Put tag alignments into the tagAlignments table
                TagCorrelationInfo tci=entry.getValue();
                int ind=1;
 
                tagTagCorrelationInsertPS.setInt(ind++, tagTagIDMap.get(entry.getKey()));
                tagTagCorrelationInsertPS.setInt(ind++, tagTagIDMap.get(tci.tag2()));
                tagTagCorrelationInsertPS.setDouble(ind++, tci.t1t2_pearson());
                tagTagCorrelationInsertPS.setDouble(ind++, tci.t1t2_spearman()); 
                tagTagCorrelationInsertPS.setDouble(ind++, tci.pres_abs_pearson());  // presence/absence vector matrix Pearson result
                tagTagCorrelationInsertPS.setDouble(ind++, tci.r2()); // presence/absence vector matrix r-squared results
                
                tagTagCorrelationInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                   // System.out.println("putTagAlignments next"+batchCount);
                    tagTagCorrelationInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            tagTagCorrelationInsertPS.executeBatch();
            connection.setAutoCommit(true);
            // print some metrics for debugging
            ResultSet rs = connection.createStatement().executeQuery("select count (*) as numCorrelations from tagCorrelations");
            if (rs.next()) {
                System.out.println("Total tag-tag correlations in tagCorrelations table: " + rs.getInt("numCorrelations"));
            }
 
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    // A tag may show up as tag1 or tag2 in the tagCorrelations table.
    // for a tagX/tagY correlation, the table will contain either an entry for tag1=x, tag2=Y, or
    // an entry for tag1=y/tag2=x.  The info is the same - only 1 entry is present.
    @Override
    public Multimap<Tag, TagCorrelationInfo> getCorrelationsForTags(List<Tag> tags) {
        ImmutableMultimap.Builder<Tag,TagCorrelationInfo> tagCorBuilder = ImmutableMultimap.builder();
        
        //loadTagHash(); // get updated tag map
        //loadRefTagHash();
 
        try {
//            ResultSet rs = connection.createStatement().executeQuery("select count(*) from tagCorrelations");
//            int size=rs.getInt(1);
//            System.out.println("size of all tagCorrelations in tagCorrelations table=" + size);
            for (Tag tag: tags){
                // For each tag on the list, get all its alignments
                Integer tagID = tagTagIDMap.get(tag);
                if (tagID == null) {
                    // what is best to print out here?
                    // should this return a NULL to indicate an error with the tag?
                    System.out.println("getCorrelationsForTag: no tagID in tagCorrelations table for tag: " + tag.sequence());
                    continue;
                }
                // The query specifies that tag1 is not a reference tag
                tagCorrelationsForTag1PS.setInt(1,tagID);
                ResultSet rs = tagCorrelationsForTag1PS.executeQuery();
                int numCorrelations = 0;
                while (rs.next()) {
                    numCorrelations++;                  
                    int tag2id = rs.getInt("tag2id");
                    double t1t2_p = rs.getFloat("t1t2_pearson");
                    double t1t2_s = rs.getFloat("t1t2_spearman");
                    double pa_pearson = rs.getFloat("pres_abs_pearson");
                    double r2 = rs.getFloat("r2");
                    
                    Tag tag2 = tagTagIDMap.inverse().get(tag2id);
                    TagCorrelationInfo tci = new TagCorrelationInfo(tag2, t1t2_p, t1t2_s, pa_pearson,r2);
 
                    tagCorBuilder.put(tag,tci);
                } 
                //System.out.println("LCJ - RepGenSQLite:getCorrelationsForTag - num correlations found for tag as tag1" + numCorrelations);
                // grab the correlations when tag2 matches the tagID
                tagCorrelationsForTag2PS.setInt(1,tagID);
                rs = tagCorrelationsForTag2PS.executeQuery();
                int numCorrelations2 = 0;
                while (rs.next()) {
                    numCorrelations2++;
                    
                    int tag1id = rs.getInt("tag1id");
                    double t1t2_p = rs.getFloat("t1t2_pearson");
                    double t1t2_s = rs.getFloat("t1t2_spearman");
                    double pa_pearson = rs.getFloat("pres_abs_pearson");
                    double r2 = rs.getFloat("r2");
                    
                    Tag tag2 = tagTagIDMap.inverse().get(tag1id);
                    TagCorrelationInfo tci = new TagCorrelationInfo(tag2, t1t2_p, t1t2_s, pa_pearson,r2);
 
                    tagCorBuilder.put(tag,tci);
                } 
//                System.out.println("LCJ - RepGenSQLite:getCorrelationsForTag - num correlations found for tag as tag1 " + numCorrelations
//                        + ", num correlations found for tag as tag2 " + numCorrelations2);
            }
        } catch (SQLException exc) {
            System.out.println("getAllTaxaMap: caught SQLException attempting to grab taxa Distribution ");
            exc.printStackTrace();
        }
        return tagCorBuilder.build();       
    }

    @Override
    public void putTagTagStats(Map<Tuple<Tag, Tag>, Tuple<Double, String>> tagTagStats, String method) {
        System.out.println("putTagTagStats");
        try {

            int batchCount=0;
            if(physicalMapPositionToIDMap==null) loadPhysicalMapPositionHash();
            connection.setAutoCommit(false);
            System.out.println("putTagTagStats: size of stats: " + tagTagStats.size());
            PreparedStatement statInsertPS=connection.prepareStatement(
                    "INSERT into tag_tag_stats (tag1id, tag2id, method_id, stat_value, stats) values(?,?,?,?,?)");
            int method_id=mappingApproachToIDMap.get(method);
            for (Map.Entry<Tuple<Tag, Tag>, Tuple<Double, String>> entry : tagTagStats.entrySet()) {
                statInsertPS.setInt(1, tagTagIDMap.get(entry.getKey().getX()));
                statInsertPS.setInt(2, tagTagIDMap.get(entry.getKey().getY()));
                statInsertPS.setInt(3, method_id);
                statInsertPS.setDouble(4, entry.getValue().getX());
                statInsertPS.setString(5, entry.getValue().getY());
                System.out.printf("%d\t%d\t%d\t%g\n",tagTagIDMap.get(entry.getKey().getX()),tagTagIDMap.get(entry.getKey().getY()),method_id, entry.getValue().getX());
                statInsertPS.addBatch();
                batchCount++;
                if(batchCount>10000) {
                    System.out.println("putTagTagStats next"+batchCount);
                    statInsertPS.executeBatch();
                    batchCount=0;
                }
            }
            statInsertPS.executeBatch();
            connection.setAutoCommit(true);
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    @Override
    public ImmutableMultimap<Tuple<Tag, Tag>, Tuple<Float, String>> getTag_tag_stats(int method_id) {
        ImmutableMultimap.Builder<Tuple<Tag,Tag>, Tuple<Float,String>> tagStatData = ImmutableMultimap.builder();
        try{
            // Get all entries with specified method id.
            String query = "select * from tag_tag_stats where method_id = " + method_id;
            ResultSet rs=connection.createStatement().executeQuery(query);
            while (rs.next()) {
                Tag tag1 = tagTagIDMap.inverse().get(rs.getInt("tag1id"));
                Tag tag2 = tagTagIDMap.inverse().get(rs.getInt("tag2id"));
                Tuple<Tag,Tag> tagTagTuple = new Tuple<Tag,Tag> (tag1,tag2);
                Float value = rs.getFloat("stat_value");
                String otherStats = rs.getString("stats");
                Tuple<Float,String> data = new Tuple<Float,String> (value,otherStats);
                tagStatData.put(tagTagTuple, data);
            }
            return tagStatData.build();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    // This method grabs the chrom, position and taxa distribution data for all tags
    // that aligned perfectly to a reference tag.
    @Override
    public ImmutableMultimap<Tag,Tuple<Position,TaxaDistribution>> getPositionTaxaDistForTag() {
        ImmutableMultimap.Builder<Tag,Tuple<Position,TaxaDistribution>> ptdBuilder = ImmutableMultimap.builder();
        String query = "select tag.sequence, tag.seqlen,physicalMapPosition.chromosome, physicalMapPosition.physical_position, physicalMapPosition.strand, tagtaxadistribution.* " +
                "from tag, physicalMapPosition, tagtaxadistribution, tagMapping where tagMapping.tagid = tag.tagid and " +
                "physicalMapPosition.posid = tagMapping.posid and tag.tagid=tagtaxadistribution.tagid ";

        // The code below grabs the requested data, creates the Tag, POsition and TaxaDistribution
        // objects, and returns the data in a map.
        int numValues = 0;
        try {
            ResultSet rs=connection.createStatement().executeQuery("select count(*) from tagMapping");
            int size=rs.getInt(1);
            rs=connection.createStatement().executeQuery("select count(*) from physicalMapPosition");
            int pmtsize = rs.getInt(1);
            System.out.println("size of all entries in tagMapping table=" + size + ", size of physicalMapPosition " + pmtsize);
            rs = connection.createStatement().executeQuery(query);

            while (rs.next()){
                numValues++;
                Tag tag = TagBuilder.instance(rs.getBytes("sequence"),rs.getShort("seqlen")).build();
                Chromosome chrom = new Chromosome(rs.getString("chromosome"));
                int posInt = rs.getInt("physical_position");
                int strand = rs.getInt("strand");
                Position pos=new GeneralPosition
                        .Builder(chrom,posInt)
                        .strand((byte)strand)
                        .build();
                TaxaDistribution td = TaxaDistBuilder.create(rs.getBytes("depthsRLE"));
                Tuple<Position,TaxaDistribution> posTD = new Tuple<Position,TaxaDistribution>(pos,td);
                ptdBuilder.put(tag,posTD);
            }
            System.out.println("\nGetPositionTaxaDistForTag:  number of query values returned: " + numValues );

        } catch (Exception exc) {
            exc.printStackTrace();
        }
        return ptdBuilder.build();
    }

    @Override
    public Tag getTagFromIndexTemp(int tagIndexInDB) {
        return tagTagIDMap.inverse().get(tagIndexInDB);
    }
}
