DROP TABLE interaction_input;

CREATE TABLE interaction_input (
    id_cp_interaction        REAL,
    partner_a                TEXT,
    partner_b                TEXT,
    protein_name_a           TEXT,
    protein_name_b           TEXT,
    annotation_strategy      TEXT,
    source                   TEXT,
    is_ppi                   INTEGER,
    curator                  TEXT,
    reactome_complex         TEXT,
    reactome_reaction        TEXT,
    reactome_pathway         TEXT,
    complexportal_complex    REAL,
    comments                 TEXT,
    version                  TEXT,
    interactors              TEXT,
    classification           TEXT,
    directionality           TEXT,
    modulatory_effect        TEXT
);

INSERT INTO interaction_input
SELECT * FROM interaction_input_old;

UPDATE interaction_input
SET id_cp_interaction = (
    SELECT it.id_cp_interaction
    FROM interaction_table it
    JOIN multidata_table md1 ON md1.id_multidata = it.multidata_1_id
    JOIN multidata_table md2 ON md2.id_multidata = it.multidata_2_id
    WHERE
        (
            md1.name = interaction_input.partner_a AND
            md2.name = interaction_input.partner_b
        )
        OR
        (
            md1.name = interaction_input.partner_b AND
            md2.name = interaction_input.partner_a
        )
    LIMIT 1
);

SELECT COUNT(*) AS filled
FROM interaction_input
WHERE id_cp_interaction IS NOT NULL;

SELECT COUNT(*) AS missing
FROM interaction_input
WHERE id_cp_interaction IS NULL;

ALTER TABLE interaction_input RENAME TO interaction_input_old2;

CREATE TABLE interaction_input (
    id_cp_interaction        REAL PRIMARY KEY,
    partner_a                TEXT,
    partner_b                TEXT,
    protein_name_a           TEXT,
    protein_name_b           TEXT,
    annotation_strategy      TEXT,
    source                   TEXT,
    is_ppi                   INTEGER,
    curator                  TEXT,
    reactome_complex         TEXT,
    reactome_reaction        TEXT,
    reactome_pathway         TEXT,
    complexportal_complex    REAL,
    comments                 TEXT,
    version                  TEXT,
    interactors              TEXT,
    classification           TEXT,
    directionality           TEXT,
    modulatory_effect        TEXT
);

INSERT INTO interaction_input
SELECT * FROM interaction_input_old2;

SELECT id_cp_interaction, COUNT(*)
FROM interaction_input_old2
GROUP BY id_cp_interaction
HAVING COUNT(*) > 1;

DROP TABLE interaction_input_old2;

PRAGMA table_info(protein_table);

ALTER TABLE protein_table RENAME TO protein_table_old;

CREATE TABLE protein_table (
    id_protein             INTEGER PRIMARY KEY,
    protein_name           TEXT,
    tags                   TEXT,
    tags_reason            TEXT,
    tags_description       TEXT,
    protein_multidata_id   INTEGER,
    FOREIGN KEY (protein_multidata_id) REFERENCES multidata_table(id_multidata)
);

INSERT INTO protein_table
SELECT * FROM protein_table_old;

DROP TABLE protein_table_old;

ALTER TABLE protein_table RENAME TO protein_table_old;


PRAGMA table_info(protein_input);

ALTER TABLE protein_input RENAME TO protein_input_old;

CREATE TABLE protein_table (
    id_protein             INTEGER PRIMARY KEY,
    protein_name           TEXT UNIQUE,  -- ? now enforced as unique
    tags                   TEXT,
    tags_reason            TEXT,
    tags_description       TEXT,
    protein_multidata_id   INTEGER,
    FOREIGN KEY (protein_multidata_id) REFERENCES multidata_table(id_multidata)
);

INSERT INTO protein_table
SELECT * FROM protein_table_old;

CREATE TABLE protein_input (
    uniprot                    TEXT,
    protein_name               TEXT,
    transmembrane              INTEGER,
    peripheral                 INTEGER,
    secreted                   INTEGER,
    secreted_desc              TEXT,
    secreted_highlight         INTEGER,
    receptor                   INTEGER,
    receptor_desc              TEXT,
    integrin                   INTEGER,
    biosynthetic_pathway_desc  TEXT,
    biosynthetic_pathway       INTEGER,
    other                      INTEGER,
    other_desc                 TEXT,
    tags                       TEXT,
    tags_reason                TEXT,
    tags_description           TEXT,
    pfam                       TEXT,
    version                    TEXT,
    FOREIGN KEY (protein_name) REFERENCES protein_table(protein_name)
);

INSERT INTO protein_input
SELECT * FROM protein_input_old;

DROP TABLE protein_table_old;
DROP TABLE protein_input_old;

PRAGMA table_info(complex_table_old);

ALTER TABLE complex_table RENAME TO complex_table_old;

CREATE TABLE complex_table (
    id_complex               INTEGER PRIMARY KEY,
    complex_multidata_id     INTEGER,
    pdb_structure            TEXT,
    pdb_id                   TEXT,
    stoichiometry            TEXT,
    comments_complex         TEXT,
    reactome_reaction        TEXT,
    reactome_complex         TEXT,
    complexportal_complex    REAL,
    rhea_reaction            TEXT,
    FOREIGN KEY (complex_multidata_id) REFERENCES multidata_table(id_multidata),
    UNIQUE (id_complex, complex_multidata_id)  -- Composite constraint
);


INSERT INTO complex_table
SELECT * FROM complex_table_old;

DROP TABLE complex_table_old;

PRAGMA table_info(complex_composition_table_old);

ALTER TABLE complex_composition_table RENAME TO complex_composition_table_old;

CREATE TABLE complex_composition_table (
	id_complex_composition INTEGER PRIMARY KEY,
    complex_multidata_id   INTEGER,
    protein_multidata_id   INTEGER,
    total_protein          INTEGER,    -- keep any other fields you use
    FOREIGN KEY (complex_multidata_id) REFERENCES multidata_table(id_multidata),
    FOREIGN KEY (protein_multidata_id) REFERENCES multidata_table(id_multidata),
	UNIQUE (complex_multidata_id, protein_multidata_id)
);

INSERT INTO complex_composition_table
SELECT * FROM complex_composition_table_old;

DROP TABLE complex_composition_table_old;

PRAGMA table_info(complex_input);

ALTER TABLE complex_input RENAME TO complex_input_old;

CREATE TABLE complex_input (
    complex_multidata_id     INTEGER,
    complex_name             TEXT,
    uniprot_1                TEXT,
    uniprot_2                TEXT,
    uniprot_3                TEXT,
    uniprot_4                TEXT,
    uniprot_5                TEXT,
    transmembrane            INTEGER,
    peripheral               INTEGER,
    secreted                 INTEGER,
    secreted_desc            TEXT,
    secreted_highlight       INTEGER,
    receptor                 INTEGER,
    receptor_desc            TEXT,
    integrin                 INTEGER,
    other                    INTEGER,
    other_desc               TEXT,
    pdb_id                   TEXT,
    pdb_structure            TEXT,
    stoichiometry            TEXT,
    comments_complex         TEXT,
    reactome_reaction        TEXT,
    reactome_complex         TEXT,
    complexportal_complex    REAL,
    rhea_reaction            TEXT,
    curator                  TEXT,
    version                  TEXT,
    FOREIGN KEY (complex_multidata_id) REFERENCES multidata_table(id_multidata)
);

INSERT INTO complex_input (
    complex_multidata_id,
    complex_name,
    uniprot_1, uniprot_2, uniprot_3, uniprot_4, uniprot_5,
    transmembrane, peripheral, secreted, secreted_desc, secreted_highlight,
    receptor, receptor_desc, integrin, other, other_desc,
    pdb_id, pdb_structure, stoichiometry, comments_complex,
    reactome_reaction, reactome_complex, complexportal_complex,
    rhea_reaction, curator, version
)
SELECT
    mt.id_multidata,
    ci.complex_name,
    ci.uniprot_1, ci.uniprot_2, ci.uniprot_3, ci.uniprot_4, ci.uniprot_5,
    ci.transmembrane, ci.peripheral, ci.secreted, ci.secreted_desc, ci.secreted_highlight,
    ci.receptor, ci.receptor_desc, ci.integrin, ci.other, ci.other_desc,
    ci.pdb_id, ci.pdb_structure, ci.stoichiometry, ci.comments_complex,
    ci.reactome_reaction, ci.reactome_complex, ci.complexportal_complex,
    ci.rhea_reaction, ci.curator, ci.version
FROM complex_input_old ci
JOIN multidata_table mt ON mt.name = ci.complex_name
WHERE mt.is_complex = 1;

DROP TABLE complex_input_old;

PRAGMA table_info(gene_table);

ALTER TABLE gene_table RENAME TO gene_table_old;

CREATE TABLE gene_table (
    id_gene        INTEGER PRIMARY KEY,
    ensembl        TEXT,  
    gene_name      TEXT,
    hgnc_symbol    TEXT,
    protein_name   TEXT,
    protein_id     INTEGER,
    FOREIGN KEY (protein_id) REFERENCES protein_table(id_protein)
);

INSERT INTO gene_table
SELECT * FROM gene_table_old;

DROP TABLE gene_table_old;
