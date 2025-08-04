WITH non_complex_genes AS (
  SELECT
    md.id_multidata,
    go.mgi_symbol
  FROM multidata_table md
  JOIN protein_table pt ON pt.protein_multidata_id = md.id_multidata
  JOIN gene_table gt ON gt.protein_id = pt.id_protein
  JOIN gene_orthologue go ON go.hgnc_symbol = gt.hgnc_symbol
  WHERE md.is_complex = 0
),
complex_genes AS (
  SELECT
    cc.complex_multidata_id AS id_multidata,
    go.mgi_symbol
  FROM complex_composition_table cc
  JOIN protein_table pt ON pt.protein_multidata_id = cc.protein_multidata_id
  JOIN gene_table gt ON gt.protein_id = pt.id_protein
  JOIN gene_orthologue go ON go.hgnc_symbol = gt.hgnc_symbol
),
all_genes AS (
  SELECT * FROM non_complex_genes
  UNION ALL
  SELECT * FROM complex_genes
),
ligands AS (
  SELECT
    it.id_cp_interaction,
    ag.mgi_symbol AS ligand_mgi
  FROM interaction_table it
  JOIN all_genes ag ON ag.id_multidata = it.multidata_1_id
  WHERE it.directionality = 'Ligand-Receptor'
),
receptors AS (
  SELECT
    it.id_cp_interaction,
    ag.mgi_symbol AS receptor_mgi
  FROM interaction_table it
  JOIN all_genes ag ON ag.id_multidata = it.multidata_2_id
  WHERE it.directionality = 'Ligand-Receptor'
)
SELECT
  l.id_cp_interaction,
  l.ligand_mgi,
  r.receptor_mgi
FROM ligands l
JOIN receptors r ON l.id_cp_interaction = r.id_cp_interaction
WHERE l.ligand_mgi IS NOT NULL AND r.receptor_mgi IS NOT NULL
ORDER BY l.id_cp_interaction;
