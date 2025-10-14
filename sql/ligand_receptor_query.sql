WITH non_complex_genes AS (
  SELECT DISTINCT
         md.id_multidata,
         TRIM(go.mgi_symbol)  AS mgi_symbol,
         TRIM(go.hgnc_symbol) AS hgnc_symbol
  FROM multidata_table md
  JOIN protein_table pt   ON pt.protein_multidata_id = md.id_multidata
  JOIN gene_table gt      ON gt.protein_id = pt.id_protein
  JOIN gene_orthologue go ON go.hgnc_symbol = gt.hgnc_symbol
  WHERE md.is_complex = 0
    AND go.mgi_symbol IS NOT NULL
    AND TRIM(go.mgi_symbol) <> ''
),
complex_genes AS (
  SELECT DISTINCT
         cc.complex_multidata_id AS id_multidata,
         TRIM(go.mgi_symbol)  AS mgi_symbol,
         TRIM(go.hgnc_symbol) AS hgnc_symbol
  FROM complex_composition_table cc
  JOIN protein_table pt   ON pt.protein_multidata_id = cc.protein_multidata_id
  JOIN gene_table gt      ON gt.protein_id = pt.id_protein
  JOIN gene_orthologue go ON go.hgnc_symbol = gt.hgnc_symbol
  WHERE go.mgi_symbol IS NOT NULL
    AND TRIM(go.mgi_symbol) <> ''
),
gene_map AS (
  SELECT DISTINCT id_multidata, mgi_symbol, hgnc_symbol
  FROM (
    SELECT * FROM non_complex_genes
    UNION ALL
    SELECT * FROM complex_genes
  ) u
),
base_interactions AS (
  SELECT DISTINCT id_cp_interaction, multidata_1_id, multidata_2_id
  FROM interaction_table
  WHERE directionality = 'Ligand-Receptor'
),
all_components AS (
  SELECT DISTINCT
         bi.id_cp_interaction,
         'ligand'       AS component_type,
         gm.mgi_symbol,
         gm.hgnc_symbol
  FROM base_interactions bi
  JOIN gene_map gm ON gm.id_multidata = bi.multidata_1_id

  UNION

  SELECT DISTINCT
         bi.id_cp_interaction,
         'receptor'     AS component_type,
         gm.mgi_symbol,
         gm.hgnc_symbol
  FROM base_interactions bi
  JOIN gene_map gm ON gm.id_multidata = bi.multidata_2_id
)
SELECT *
FROM all_components
ORDER BY id_cp_interaction, component_type, mgi_symbol;

