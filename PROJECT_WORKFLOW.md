# Project workflow: neem–Alternaria public data and exploratory protein ML

## 1. Project aim and current scope

This project assembles public biological data for the fungal pathogen *Alternaria alternata* and the neem plant *Azadirachta indica*. The collected data types are genomic, transcriptomic, proteomic, and metabolite/chemical. The computational component uses the predicted *A. alternata* protein sequences to demonstrate sequence feature engineering, local parallel processing, and exploratory unsupervised machine learning.

The work completed so far has two layers:

1. A provenance-aware preliminary dataset: selected public records were downloaded or exported, organized, counted, and described in metadata tables.
2. An exploratory AI/ML analysis of the fungal predicted proteome: 13,466 proteins were converted to numerical sequence features and grouped with K-means; PCA was used to visualize the feature space.

This is not yet an integrated biological interaction model. No result here establishes that a protein is an antifungal target, that a neem compound binds a protein, or that a compound has antifungal activity. Those would require additional annotation, evidence, and validation.

## 2. Repository layout and data handling

The project keeps source files and project-generated summaries separate.

- data/raw/alternaria/genome/: downloaded NCBI Datasets genome package and the user-supplied genome archive.
- data/raw/alternaria/proteome/: convenient copy of the fungal protein FASTA used by the feature script.
- data/raw/alternaria/transcriptome/: NCBI SRA Run Selector metadata export.
- data/raw/neem/genome/: neem assembly package and supplied compressed WGS/GenBank records.
- data/raw/neem/proteome/: UniProtKB neem FASTA export.
- data/raw/neem/transcriptome/: SRA Run Selector exports for the selected neem BioProjects.
- data/raw/compounds/: focused PubChem SDF and an invalid/incomplete COCONUT transfer retained with a .partial suffix.
- data/processed/: curated metadata tables and the master file inventory.
- data/processed/annotation/alternaria/: AI/ML feature table, benchmark reports, cluster assignments, cluster summary, and PCA image.
- scripts/objective2/: executable protein feature extraction and clustering scripts.
- notebooks/: created project folder; no notebook files are currently part of the workflow.
- tools/python311/: contains an embeddable Python archive downloaded during environment setup attempts; it is not biological input data and is not required by the completed analysis.
- .venv/: project-local Windows Python environment folder. Its configuration points to Python 3.14.3; package versions are not recorded in a requirements file.

Raw files preserve downloaded or source-exported data. Processed CSV files organize and summarize those data without replacing the raw files. The manifest gives a file-level index and marks the COCONUT partial archive as invalid.

## 3. Data collection workflow

### 3.1 *Alternaria alternata* genome and predicted proteins

The selected fungal assembly is NCBI RefSeq GCF_001642055.1, assembly name Altal1, strain SRC1lrK2f. It is a scaffold-level assembly. Its NCBI Datasets package is under data/raw/alternaria/genome/ncbi_dataset/data/GCF_001642055.1/.

The package provides:

- GCF_001642055.1_Altal1_genomic.fna: genome nucleotide FASTA, 79 sequence records.
- genomic.gff: feature/gene annotation on the same assembly.
- cds_from_genomic.fna: 13,466 coding-sequence FASTA records.
- protein.faa: 13,466 predicted protein FASTA records.
- assembly_data_report.jsonl and dataset_catalog.json: package metadata and file catalog.

The protein FASTA was copied without sequence modification to data/raw/alternaria/proteome/alternaria_alternata_GCF_001642055.1_protein.faa. That copy is the input to scripts/objective2/02_protein_features.py. The user-supplied alternaria_alternata_genome.zip is also retained. It was inspected as a valid archive containing both the matching GenBank GCA_001642055.1 and RefSeq GCF_001642055.1 package representations, with annotation files including GTF. These are two database representations of the same assembly lineage, not two independent fungal genomes.

### 3.2 Neem genome and annotation availability

The selected neem assembly is NCBI GenBank GCA_022749755.1, corresponding to WGS accession JAGQDM000000000.1 and BioProject PRJNA645650. It is reported as chromosome level, with 70 sequence records totaling 281,703,417 bases and 14 chromosome assignments.

The NCBI package is under data/raw/neem/genome/ncbi_dataset/data/GCA_022749755.1/. The user-provided files JAGQDM01.1.fsa_nt.gz and JAGQDM01.1.gbff.gz are retained beside that package. The compressed FASTA was checked against the assembly: it represents 70 sequence records with the same total base count. These files are alternate records/representations of the same genome project, not additional independent neem assemblies.

The selected neem package does not provide a GFF, CDS FASTA, or protein FASTA. The WGS record reports zero proteins. Therefore the neem genome is present, but it is not currently a genome-annotated neem gene/protein set. This distinction matters: the separate UniProt records below must not be described as proteins predicted from GCA_022749755.1.

### 3.3 Transcriptome metadata and selection

The project collected SRA run metadata rather than raw sequencing reads. The untouched Run Selector exports are:

- data/raw/alternaria/transcriptome/PRJNA813997_runinfo.csv
- data/raw/neem/transcriptome/PRJNA645650_runinfo.csv
- data/raw/neem/transcriptome/PRJNA176672_runinfo.csv
- data/raw/neem/transcriptome/PRJNA344442_runinfo.csv

The curated data/processed/transcriptomic_metadata.csv contains 34 run-level entries: 18 *A. alternata* and 16 neem. Each record preserves run, experiment, sample and BioSample identifiers where present, library strategy/source, layout, sequencing platform/instrument, source link, and notes. Rows are at the run level; a run is not automatically an independent biological sample.

The selection retains transcriptomic-source rows and preserves unusual strategies with explicit cautions:

- *A. alternata*, PRJNA813997: 18 RNA-Seq runs. Source sample/condition labels are retained; their biological interpretation should be checked against the study methods before treating labels as experimental groups.
- Neem, PRJNA645650: two RNA-Seq records, one Illumina paired-end run and one PacBio single-molecule run. They are distinct experiments and should not be assumed to be biological replicates of one another.
- Neem, PRJNA176672: 12 RNA-Seq runs plus one WGA-strategy row whose library source is marked TRANSCRIPTOMIC. The 12 RNA-Seq runs map to the same experiment SRX838026 and BioSample SAMN02715695; they are not 12 independent biological replicates. The WGA row is preserved as source metadata, but is not conventional RNA-Seq and should not enter an expression analysis without study-specific justification.
- Neem, PRJNA344442: one CLONE-strategy record with transcriptomic library source. It is a cDNA clone library, not conventional bulk RNA-Seq; its sample name contains a source spelling error.

The historical paper label SRA053330 is recorded as a legacy reference, but is not silently treated as equivalent to a modern BioProject accession. No SRA reads have been downloaded and no expression matrix or differential-expression results have been generated. The metadata table is a catalog for choosing runs in a later stage, not an expression dataset.

### 3.4 Protein data

The curated data/processed/proteomic_metadata.csv describes two protein resources:

- *A. alternata*: 13,466 predicted RefSeq proteins from GCF_001642055.1, matching the CDS record count.
- *A. indica*: 39 UniProtKB entries returned by a species/taxonomy-filtered query (taxonomy ID 124943).

The 39 neem records are a useful sequence subset with identifiers and UniProt annotations, but they are not a complete genome-derived neem proteome. Their evidence and completeness vary, and they are not tied to the selected chromosome-level assembly. The files and provenance are retained in data/raw/alternaria/proteome/ and data/raw/neem/proteome/.

Here, “proteomic data” means protein sequence records and their available annotations. No experimental mass-spectrometry abundance, peptide identification, or protein-expression dataset has been collected.

### 3.5 Neem compounds and broader natural-product coverage

The PubChem structure file data/raw/compounds/neem_limonoids_pubchem.sdf contains five records. The processed data/processed/compound_metadata.csv associates each structure with PubChem CID, title, molecular formula, molecular weight, canonical/connectivity SMILES, isomeric SMILES, InChIKey, source link, and local structure-file path.

| Compound | PubChem CID | Formula | Molecular weight |
| --- | ---: | --- | ---: |
| Azadirachtin | 5281303 | C35H44O16 | 720.7 |
| Nimbolide | 12313376 | C27H30O7 | 466.5 |
| Nimbin | 108058 | C30H36O9 | 540.6 |
| Salannin | 6437066 | C34H44O9 | 596.7 |
| Gedunin | 12004512 | C28H34O7 | 482.6 |

These five limonoids are a focused seed set, not an exhaustive neem metabolome or chemical library. COCONUT October 2026 2D-lite SDF was attempted as a broader natural-product resource. The resulting data/raw/compounds/coconut_sdf_2d_lite-10-2026.zip.partial failed ZIP integrity validation after concurrent writes. It remains marked invalid in the manifest and must not be extracted or used as a dataset. COCONUT is therefore not currently included as usable data.

The compound table contains chemical identity and structure metadata, not measured neem metabolite concentrations, tissue-specific abundance, or experimental antifungal potency.

### 3.6 Provenance tables

The processed metadata files serve different purposes:

- genomic_metadata.csv records genome/annotation file type, assembly accession, path, sequence count where measured, and assembly-specific notes.
- transcriptomic_metadata.csv describes each selected SRA run and makes replicate/strategy caveats explicit.
- proteomic_metadata.csv records source/query, accession, local FASTA path, sequence count and completeness caveat.
- compound_metadata.csv records chemical identifiers, structure strings, properties and source links.
- dataset_manifest.csv is the master raw/processed file inventory. It is a file catalog; detailed biological fields remain in the specialist metadata tables.

The specialist tables were assembled from the saved NCBI Run Selector exports, downloaded NCBI package reports/files, the UniProt FASTA query, the PubChem structure/property records, and direct file checks. The raw exports and accessions preserve the source evidence. There is currently no checked-in acquisition/curation script that rebuilds every metadata CSV from the remote services; the curated CSVs themselves are the saved project outputs.

Primary source portals and records include [NCBI RefSeq GCF_001642055.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001642055.1/), [NCBI GenBank GCA_022749755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022749755.1/), [SRA BioProjects PRJNA813997](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA813997), [PRJNA645650](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA645650), [PRJNA176672](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA176672), and [PRJNA344442](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA344442), the [UniProtKB taxonomy query](https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id%3A124943&format=fasta), the five CID-linked [PubChem records](https://pubchem.ncbi.nlm.nih.gov/), and the [COCONUT download portal](https://coconut.naturalproducts.net/download). Accession-specific PubChem links are stored row-by-row in compound_metadata.csv.

## 4. AI/ML workflow on the fungal proteome

### 4.1 Input and processing scope

The AI/ML input is the predicted *A. alternata* protein FASTA at data/raw/alternaria/proteome/alternaria_alternata_GCF_001642055.1_protein.faa. It is not the fungal genome FASTA, CDS FASTA, or the neem UniProt file. The feature table contains 13,466 protein rows, matching the 13,466 input protein records. The current feature CSV has 28 columns: protein ID plus 27 numerical features. There are no empty feature cells in the produced table. In this input, sequence_length and standard_aa_length match for all records, indicating no residues were discarded by standard-amino-acid filtering.

### 4.2 Feature engineering

scripts/objective2/02_protein_features.py uses Biopython FASTA parsing and calculates features independently for each sequence:

- sequence_length: original FASTA sequence length.
- standard_aa_length: sequence length after retaining the 20 standard amino acids ACDEFGHIKLMNPQRSTVWY.
- frac_A through frac_Y: fraction of each standard amino acid among retained residues.
- hydrophobic_fraction: fraction of residues from AVILMFWY.
- positive_fraction: fraction of residues from KRH.
- negative_fraction: fraction of residues from DE.
- aromatic_fraction: fraction of residues from FWY.
- sequence_entropy: Shannon entropy of the standard-amino-acid composition, in bits.

For an empty cleaned sequence, the script sets fraction/entropy values to zero. Protein IDs are retained to allow every feature row to be traced back to its FASTA record. These are composition and length descriptors; they do not encode domains, cellular locations, catalytic sites, expression, virulence, or experimentally measured activity.

### 4.3 Serial/parallel processing benchmark

The feature script calculates the features serially and with Python multiprocessing, then checks that the serial and parallel protein-ID sequences agree before writing the parallel result. The recorded configuration used four workers and three repeats on a machine reporting 20 CPU cores. The JSON report records the individual timings, averages, speedup, and efficiency.

| Measurement | Recorded value |
| --- | ---: |
| Proteins | 13,466 |
| Workers | 4 |
| Repeats per mode | 3 |
| Mean serial time | 0.5440 s |
| Mean parallel time | 0.8737 s |
| Parallel/serial speedup | 0.6226× |
| Parallel efficiency | 0.1556 |

The speedup below 1 means the measured parallel run was slower than serial for this workload on this machine. A likely explanation is that process startup and task-transfer overhead exceeded the savings from distributing these inexpensive calculations. This is a local CPU multiprocessing benchmark, not a distributed HPC-cluster benchmark. The output names include “hpc” because the project compares serial and parallel processing; the measurements do not demonstrate execution on an HPC system.

### 4.4 Unsupervised clustering and PCA visualization

scripts/objective2/03_protein_clustering.py reads the feature CSV from the project root. It excludes protein_id, replaces infinite values with missing values, and removes rows containing missing numerical features. It then standardizes each numerical feature with StandardScaler and fits K-means with four clusters, random_state 42, and n_init 10.

K-means is fit in the full standardized numerical feature space. PCA is separately fit to that same standardized feature matrix to produce two coordinates for the plot. Thus PCA is a two-dimensional visualization; it is not the clustering algorithm and is not used to reduce the matrix before K-means.

The resulting cluster sizes are:

| Cluster label | Protein count | Mean sequence length | Mean hydrophobic fraction | Mean aromatic fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 2,054 | 430.0 | 0.337 | 0.053 |
| 1 | 2,867 | 358.5 | 0.456 | 0.106 |
| 2 | 6,407 | 533.6 | 0.394 | 0.087 |
| 3 | 2,138 | 382.9 | 0.332 | 0.063 |

These summaries describe feature-space groups only. For example, cluster 2 is the largest group and has the longest mean sequence length; cluster 1 has the highest mean hydrophobic and aromatic fractions among the four. Such differences are not functional labels. Cluster numbers are arbitrary and may change if features, scaling, random seed, or cluster count change.

The feature design has deliberate simplifications and some redundancy: sequence_length and standard_aa_length are identical in this input, while aggregate hydrophobic/charge/aromatic fractions are derived from amino-acid fractions already included. Since all numerical columns are standardized and passed to K-means, overlapping features can effectively give related composition signals additional weight. Interpret the grouping as an exploratory demonstration, not a biological classification.

### 4.5 AI/ML outputs

All generated outputs are under data/processed/annotation/alternaria/:

- protein_features.csv: one row for each input protein and its sequence-derived features.
- hpc_benchmark.json: machine-readable input, protein count, worker/repeat settings, individual and mean runtimes, speedup, efficiency, and caveat.
- hpc_benchmark_summary.txt: human-readable benchmark summary.
- protein_clusters.csv: protein-level feature rows with K-means cluster label and PCA1/PCA2 coordinates.
- cluster_summary.csv: mean numerical features and protein count for each cluster.
- protein_clusters_pca.png: scatter plot of the two PCA coordinates, colored by arbitrary cluster label.

The feature and cluster tables preserve protein IDs. The benchmark records run-specific timing values. The PCA variance percentages are printed by the clustering script at run time but are not separately stored in the listed output files.

## 5. Reproduction

Run the commands from the neem_antifungal_project root with an environment containing Biopython, pandas, NumPy, scikit-learn, and Matplotlib. The project currently has a Windows .venv folder configured for Python 3.14.3; the exact dependency versions are not pinned in a requirements file.

PowerShell example:

    cd "C:\3 - Uni\Semester 5\1 - Projects\AIML_HPC_BiologicalDataFinal\neem_antifungal_project"
    .\.venv\Scripts\Activate.ps1
    python scripts\objective2\02_protein_features.py --workers 4 --repeats 3
    python scripts\objective2\03_protein_clustering.py

The feature script defaults to the selected fungal protein FASTA and project output directory. It also accepts an alternate FASTA path with --input and an alternate destination with --output-dir. The clustering script currently uses project-root-relative input/output paths, so run it with the project root as the working directory. Re-running either script replaces its corresponding generated output files.

## 6. Validation performed and evidence trail

The collected files and generated outputs were checked for consistency at the file/record level:

- fungal genome FASTA: 79 sequence records; fungal CDS/protein FASTA: 13,466 records each.
- neem assembly FASTA and supplied WGS FASTA: 70 records and 281,703,417 bases.
- focused neem SDF: five record terminators, matching the five curated compound rows.
- transcript metadata: 34 curated rows, broken down as 18 fungal and 16 neem records.
- protein features and cluster assignments: 13,466 rows, matching the fungal protein FASTA record count.
- cluster membership: cluster sizes sum to 13,466.
- COCONUT transfer: integrity check failed; it is explicitly excluded from usable data.

Primary provenance resides in the accessions, source links, raw SRA exports, NCBI Datasets package reports, and processed metadata tables. The master manifest is the index for locating files; it should be updated whenever a raw or generated artifact is added or replaced.

## 7. Limits and next-stage work

Completed work does not include the following:

1. A complete neem proteome derived from the selected GCA_022749755.1 assembly. The 39-entry UniProt subset is separate and incomplete.
2. Downloaded RNA-seq read files, read-quality checks, mapping/quantification, or differential-expression analysis.
3. A usable COCONUT natural-product library; the current local partial archive is invalid.
4. Functional enrichment or domain-based annotation of fungal clusters.
5. Any supervised model, labeled training data, prediction of antifungal activity, target prioritization, docking, or experimental validation.
6. Cross-species gene/protein mapping or a joined compound–protein interaction table.

A defensible next stage would select appropriate biological samples and independent replicates, obtain or construct an annotation compatible with the chosen neem assembly, replace the invalid COCONUT transfer with a verified archive if required, add functional annotations/domains to the fungal proteins, and only then define a validated biological question for cross-omics or compound–target modeling. Preserve the assembly and sample identifiers when joining records, and do not treat multiple runs from one BioSample as biological replication.

## 8. Short project summary

Public genomic resources were collected for both organisms; fungal gene/CDS/protein annotation is available for the selected *A. alternata* assembly, while the selected neem genome lacks annotation files. SRA run metadata were curated without downloading reads. Five neem limonoids were linked to PubChem structures and physicochemical identifiers; the broader COCONUT attempt is invalid. The AI/ML component extracted 27 numeric descriptors per fungal protein, benchmarked serial and four-worker local processing, standardized the feature matrix, generated four exploratory K-means groups, and plotted a two-dimensional PCA projection. The resulting artifacts are traceable and ready to support later analysis, but they do not yet support claims about protein function, antifungal targets, or compound activity.
