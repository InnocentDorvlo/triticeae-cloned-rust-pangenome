# Triticeae Wheat Rust Resistance Super-Pangenome Pipeline

Analysis, statistical, and visualization scripts separated by execution stage corresponding to **Results Sections 3.1–3.5** and **Methods Sections 5.1–5.8**.

## Repository Layout

### Root Configuration & Curated Reference Inputs
- `env_3.1.yml` — Conda/Mamba software environment specification.
- `config_3.1.json` — Pipeline parameters and directory settings.
- `catalogue_3.1.tsv` — Reference catalog of 43 cloned wheat rust resistance proteins (36 loci) + 2 comparator tandem kinases (`Rpg1`, `WTK4`).
- `panel_3.1.tsv` — 27 chromosome-scale Triticeae genome assemblies (56 haploid subgenomes) and source URLs.
- `centromeres_CSv21.tsv` — Chinese Spring RefSeq v2.1 CENH3 centromere coordinates.
- `queries_3.1.faa`, `queries_3.1.cds.fna` — Curated reference protein and coding sequences.

---

### `00_setup/` — Environment Verification & Dataset Utilities
- `check_tools_3.1.sh` — Verifies Conda dependencies (`miniprot`, `diamond`, `hmmer`, `bedtools`, `orthofinder`, `Rscript`) and data endpoints.
- `install_figure_packages.sh` / `install_figure_packages.R` — Installs required CRAN and Bioconductor plotting packages (`ggplot2`, `patchwork`, `ggtree`, `ComplexHeatmap`, etc.).
- `upgrade_incomplete_datasets.py` — Optional pre-fetch/upgrade utility for spliced genomic CDS across assemblies, AlphaFold DB/ESMFold 3D PDBs, Ensembl repeat windows, and BUSCO completeness tables.

### `01_section3.1_pangenome_inventory/` — Methods 5.1–5.3 & Results 3.1
*Annotation-independent homologue discovery, BUSCO species tree, pangenome partitioning, and saturation modeling (Fig. 1).*
- **Primary pipeline & runner:** `rust31.py`, `run_3.1.sh`
- **Methods 5.1 (Query domains, BUSCO QC & 42-tip subgenome species phylogeny):** `interproscan_queries.py`, `collect_busco.py`, `markers_from_proteomes.py`, `build_species_tree.py`
- **Methods 5.2 (Paralogue pruning & syntenic anchor verification):** `prune_paralogues.py`
- **Methods 5.3 & 5.8 (Heaps' law saturation, R audit & Fig. 1):** `heaps_3.1.R`, `stats_and_fig1_3.1.R`, `make_unified_master_fig1.R`

### `02_section3.2_chromosomal_synteny/` — Methods 5.4 & Results 3.2
*Centromere-to-telomere normalization (`d_norm`), hexaploid A:B:D triad balance, ±500 kb Synteny Retention Index (SRI), and flanking TE profiling (Fig. 2).*
- **Primary pipeline & runner:** `rust32.py` (`positions`, `triads`, `sri`, `te`, `stats`), `run_3.2.sh`
- **Methods 5.8 (R audit & Fig. 2):** `stats_and_fig2_3.2.R`, `make_fig2_review.R`

### `03_section3.3_tandem_kinases/` — Methods 5.5 & Results 3.3
*Class III tandem protein kinase and kinase-fusion dissection, ML kinase phylogenies, K1–K2 `pS`/`dS` divergence, progenitor distribution, and intron phase mapping (Fig. 3).*
- **Primary pipeline & runner:** `rust33.py` (`domains`, `pairs`, `phylogeny`, `origins`, `exons`, `stats`, `figure`), `run_3.3.sh`
- **Methods 5.8 (R audit & Fig. 3):** `stats_and_fig3_3.3.R`, `make_publication_fig3.R`

### `04_section3.4_selection_and_3d/` — Methods 5.6 & Results 3.4
*MACSE codon alignments, IQ-TREE gene trees, PAML `codeml` M0/site models, HyPhy FUBAR/MEME selection scans, Class III module asymmetry, and AlphaFold + DSSP/Shrake–Rupley 3D projection (Fig. 4).*
- **Primary pipeline, selection module & runners:** `rust34.py`, `rust34_real.py`, `run_3.4.sh`, `run_meme_3.4.sh`
- **3D structural rendering (PyMOL), R audit & Fig. 4:** `render_fig4_structures.py`, `stats_and_fig4_3.4.R`, `make_publication_fig4.R`

### `05_section3.5_haplotype_diversity/` — Methods 5.7 & Results 3.5
*Germplasm tier partitioning, orthogroup amino-acid diversity (`pi_aa`, `theta_W,aa`, Tajima's `D`, `FST`), MST haplotype networks, uncaptured candidate allele mining, and geographic origin analysis (Fig. 5).*
- **Primary pipeline & runner:** `rust35.py` (`extract`, `popgen`, `networks`, `diversity`, `mining`, `geography`, `bottleneck`, `stats`, `figure`), `run_3.5.sh`
- **Methods 5.8 (R audit & Fig. 5):** `stats_and_fig5_3.5.R`, `make_publication_fig5.R`

### `06_tables_and_downstream/` — Methods 5.8 Table Export & Downstream Runner
- `make_tables.py` — Exports standardized Table 1 and Supplementary Tables (`tables/`).
- `build_complete_supplementary_workbook.py` — Builds `Supplementary_Data_Tables.xlsx`.
- `reproduce_all_figures_and_tables.sh` — Runs downstream R statistical audits, Figures 1–5, table export, and supplementary Excel workbook generation from `results/`.

---

## Execution Order

Run all commands from the repository root directory:

```bash
# 0. Environment setup and tool verification
conda env create -f env_3.1.yml
conda activate rust31
bash 00_setup/check_tools_3.1.sh
bash 00_setup/install_figure_packages.sh
# Optional pre-fetch of remote Ensembl/NCBI CDS, AlphaFold PDBs, and repeat windows:
# python3 00_setup/upgrade_incomplete_datasets.py all --threads 4

# 1. Section 3.1: Pangenome inventory & homologue rescue
bash 01_section3.1_pangenome_inventory/run_3.1.sh
python3 01_section3.1_pangenome_inventory/interproscan_queries.py --email you@institute.edu
python3 01_section3.1_pangenome_inventory/markers_from_proteomes.py --subset 500
python3 01_section3.1_pangenome_inventory/build_species_tree.py
python3 01_section3.1_pangenome_inventory/collect_busco.py
python3 01_section3.1_pangenome_inventory/prune_paralogues.py
Rscript 01_section3.1_pangenome_inventory/stats_and_fig1_3.1.R
Rscript 01_section3.1_pangenome_inventory/make_unified_master_fig1.R

# 2. Section 3.2: Chromosomal position, triads & microsynteny
bash 02_section3.2_chromosomal_synteny/run_3.2.sh

# 3. Section 3.3: Tandem kinases & kinase-fusions
bash 03_section3.3_tandem_kinases/run_3.3.sh
Rscript 03_section3.3_tandem_kinases/stats_and_fig3_3.3.R
Rscript 03_section3.3_tandem_kinases/make_publication_fig3.R

# 4. Section 3.4: Selection analysis & 3D structural projection
bash 04_section3.4_selection_and_3d/run_3.4.sh
python3 04_section3.4_selection_and_3d/render_fig4_structures.py
Rscript 04_section3.4_selection_and_3d/make_publication_fig4.R

# 5. Section 3.5: Haplotype networks, diversity & candidate allele mining
bash 05_section3.5_haplotype_diversity/run_3.5.sh
Rscript 05_section3.5_haplotype_diversity/stats_and_fig5_3.5.R
Rscript 05_section3.5_haplotype_diversity/make_publication_fig5.R

# 6. Standardized tables & workbook
bash 06_tables_and_downstream/reproduce_all_figures_and_tables.sh
```
