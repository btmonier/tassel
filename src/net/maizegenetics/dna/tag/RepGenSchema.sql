-- Table: tag
CREATE TABLE tag (
    tagid    INTEGER PRIMARY KEY,
    tagName VARCHAR,
    --convert sequence to VARCHAR if possible
    sequence BLOB NOT NULL,
    sequencetext TEXT NOT NULL,
    seqlen INTEGER NOT NULL,
    --remove if reference if possible
    isReference BOOLEAN,
    qualityScore TEXT,
    numTagInstances INTEGER,
    UNIQUE (sequence,seqlen)
);

-- Table: tagTagAlignments
CREATE TABLE tag_tag_stats (
    tagtagstatId INTEGER   PRIMARY KEY,
    tag1id	INTEGER NOT NULL,
    tag2id  INTEGER NOT NULL,
    method_id INTEGER NOT NULL,
    stat_value FLOAT,
    stats VARCHAR,
    UNIQUE (tag1id, tag2id, method_id)
);

-- Table: tagMapping
-- Junction (link) table between reftag, and physicalMapPosition
CREATE TABLE tagMapping (
    reftagid       INTEGER NOT NULL,
    position_id  INTEGER NOT NULL,
    method_id    INTEGER NOT NULL,
    bp_error	INTEGER,
    cm_error	FLOAT(2),
    PRIMARY KEY (reftagid, position_id)
);

-- Table: referenceGenome
-- Identifies reference genomes for physicalMapPosition table
CREATE TABLE referenceGenome (
	refid	INTEGER PRIMARY KEY,
	refname	TEXT,
	UNIQUE (refname)
);
 
-- Table: physicalMapPosition
CREATE TABLE physicalMapPosition (
    posid INTEGER   PRIMARY KEY,
    reference_genome_id	INTEGER,
    chromosome TEXT      NOT NULL,
    physical_position   INTEGER   NOT NULL,
    strand     INTEGER(1)  NOT NULL
);
CREATE UNIQUE INDEX physpos_idx ON physicalMapPosition(chromosome,physical_position,strand);
CREATE INDEX phychrpos_idx ON physicalMapPosition(chromosome);

-- Table: tagTagAlignments
CREATE TABLE tag_tag_Alignments (
    tagAlignId INTEGER   PRIMARY KEY,
    tag1id	INTEGER,
    tag2id  INTEGER,
    score INTEGER
);
CREATE INDEX tag_tag_Alignments_idx1 on tag_tag_Alignments(tag1id,tag2id);

-- Table: tagCorrelations
-- t1t2_pearson is the pearson correlation of tag1 depths vector by tag2 depths vector
-- t1t2_spearman is the spearman correlation of tag1 depths vector by tag2 depths vector
-- pres_abs_pearson is the pearson correltaion of tag1 prime vec by tag2 prime vec (presence/absence)
-- r2 is the t1' x t2' r2 value from analysis.popgen.LinkageDisequilibrium.calculateRSqr()
CREATE TABLE tagCorrelations (
    tagcorrelationsId INTEGER   PRIMARY KEY,
    tag1id	INTEGER,
    tag2id  INTEGER,
    t1t2_pearson  real,
    t1t2_spearman  real,
    pres_abs_pearson real,
    r2  real
);
CREATE INDEX tagCorrelationsions_idx1 on tagCorrelations(tag1id,tag2id);


-- Table: mappingApproach
CREATE TABLE mappingApproach (
    mapappid INTEGER   PRIMARY KEY,
    approach TEXT NOT NULL UNIQUE,
    software   TEXT NOT NULL,
    protocols   TEXT NOT NULL
);


-- Table: tagtaxadistribution
CREATE TABLE tagtaxadistribution (
    tagtxdstid  INTEGER   PRIMARY KEY,
    tagid      INTEGER NOT NULL,
    depthsRLE  BLOB,
    totalDepth INTEGER
);
CREATE INDEX tagid_idx ON tagtaxadistribution(tagid);


-- Table: taxa
CREATE TABLE taxa (
    taxonid INTEGER PRIMARY KEY,
    name    TEXT NOT NULL
);
