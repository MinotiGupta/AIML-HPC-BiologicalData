# Neem–Alternaria public omics dataset (preliminary)

This workspace combines public genomic, transcriptomic, protein and compound resources for *Alternaria alternata* and *Azadirachta indica*. Raw source files are retained under `data/raw/`; the `data/processed/*.csv` files record source identifiers, local paths and selection notes. Retrieval date for newly obtained items: 2026-10-15.

## Current dataset contents

- **Alternaria genome and annotation:** NCBI RefSeq `GCF_001642055.1` (Altal1, strain SRC1lrK2f), scaffold assembly: 79 genomic FASTA records, GFF, 13,466 CDS records and 13,466 predicted protein records. Matching GenBank/RefSeq archive supplied by the user is retained.
- **Neem genome:** NCBI GenBank `GCA_022749755.1` / WGS `JAGQDM000000000.1`, BioProject `PRJNA645650`: 70 sequences and 281,703,417 bases, chromosome-level assembly. Its package includes no GFF/CDS/protein FASTA; the WGS flatfile reports no proteins.
- **Transcriptome metadata:** Run-level metadata exports for fungal `PRJNA813997` and neem `PRJNA645650`, `PRJNA176672`, and `PRJNA344442`; the processed table has one row per transcriptomic-source run (34 records). Reads have not been downloaded and no expression matrix has been generated. The 12 RNA-Seq runs under PRJNA176672 map to one experiment/BioSample, so they are not 12 independent biological replicates. WGA and CLONE strategy entries are explicitly flagged and should not be treated as conventional RNA-Seq without checking the study methods. The article's legacy label `SRA053330` is not silently equated with a modern project accession.
- **Proteins:** 13,466 predicted fungal RefSeq proteins, plus 39 species-filtered UniProtKB neem records. The UniProt subset is not a complete neem proteome.
- **Neem compounds:** Five limonoids with PubChem CIDs, formulas, molecular weights, canonical and isomeric SMILES, and InChIKeys in `compound_metadata.csv`. The companion PubChem SDF has five records. These are a focused seed set, not an exhaustive neem metabolome.
- **Broader natural products:** COCONUT October 2026 2D-lite SDF transfer was attempted, but the resulting partial file failed ZIP integrity validation after concurrent writes. It is marked invalid in the manifest and is not part of the usable dataset. The verified five-compound PubChem set is the current chemical dataset.

## Tables

- `genomic_metadata.csv`: assembly/accession, source, paths, counts and genomic-resource caveats.
- `transcriptomic_metadata.csv`: run, experiment, sample/BioSample, assay, layout and platform metadata; source exports remain in raw folders.
- `proteomic_metadata.csv`: protein source and completeness status.
- `compound_metadata.csv`: PubChem chemical identifiers and structures.
- `dataset_manifest.csv`: file-level inventory of raw and processed files.

## Sources

NCBI Datasets/RefSeq/GenBank and SRA; UniProtKB REST; PubChem PUG REST; COCONUT download portal. URLs and accessions are provided in the corresponding processed table or raw export. COCONUT download page: https://coconut.naturalproducts.net/download (bulk transfer was not validated; see manifest). PubChem compound records are linked by CID in `compound_metadata.csv`.

## Readiness and limitations

This is a preliminary, traceable collection for the next project stage. Before cross-omics integration, verify fungal and neem strain/tissue relevance, resolve the biological meaning of library labels from primary study methods, download selected SRA runs if expression analysis is needed, and acquire an annotated neem gene/protein set if the selected assembly can be annotated or a matching annotated assembly is found. Preserve accession-level provenance and avoid treating repeated runs from one BioSample as biological replication.

## Detailed workflow

See [PROJECT_WORKFLOW.md](PROJECT_WORKFLOW.md) for the full data-collection and AI/ML workflow, methods, outputs, results, limitations, and reproduction instructions.

