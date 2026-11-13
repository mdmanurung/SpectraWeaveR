# Cytometry analysis on GitHub — studies, reproducible workflows, gating, QC, and normalization

**Catalog assembled:** 12 November 2026  
**Purpose:** A broad, **unranked** inventory for later screening, not a recommendation or an assertion that every repository runs without modification.  
**Scope:** Flow, spectral flow, and mass cytometry. R implementations are emphasized; hybrid R/Nextflow/Python projects are retained when useful. Repositories from the preceding discussion are included, plus additional study pipelines, QC studies, benchmarks, and reference software.  
**Evidence standard:** Entries link to actual public repositories and/or published papers. Descriptions distinguish *verified code content* from broader paper claims. **Real data available** means a publicly described data source, not a guarantee that it is currently downloadable. No repositories were executed or independently benchmarked for this catalog.

## Navigation

- [A. Study-specific analyses and manuscript repositories](#a-study-specific-analyses-and-manuscript-repositories)
- [B. Standalone and adaptable multi-sample pipelines](#b-standalone-and-adaptable-multi-sample-pipelines)
- [C. QC, normalization, and gating evaluations / benchmarks](#c-qc-normalization-and-gating-evaluations--benchmarks)
- [D. Core R cytometry tools and integration frameworks](#d-core-r-cytometry-tools-and-integration-frameworks)
- [E. Tutorials and reproducible demonstrations with real datasets](#e-tutorials-and-reproducible-demonstrations-with-real-datasets)
- [F. Published applications with no confirmed public analysis repository](#f-published-applications-with-no-confirmed-public-analysis-repository)
- [G. Study-review worksheet and screening criteria](#g-study-review-worksheet-and-screening-criteria)

**Shorthand:** `AG` = automated algorithmic gates; `SG` = semi-automated or expert-guided gates; `QC` = event/sample quality control; `BC` = batch correction/normalization; `CL` = high-dimensional clustering; `ICS` = intracellular cytokine staining. “Package” and “training” entries are explicitly separated from study applications.

---

## A. Study-specific analyses and manuscript repositories

These are the strongest matches to the original request: R code used on an experimental dataset, rather than merely a new R package.

### A01 — SATVILab: TB antigen-specific T-cell activation assay (2025)

- **GitHub:** https://github.com/SATVILab/DataTidyMusvosviPredictTBTASA
- **Relevant code:** [Script to perform TB-TASA.R](https://github.com/SATVILab/DataTidyMusvosviPredictTBTASA/blob/main/Script%20to%20perform%20TB-TASA.R); manuscript analysis script in repository root.
- **Experiment:** Tuberculosis antigen-specific activation and treatment-monitoring assay; stimulated versus unstimulated T-cell readouts. **Methods:** `openCyto::gatingTemplate()`, `gt_gating()`, hierarchical scatter/singlet/viability/CD3 gates, HLA-DR, IFNγ and TNF detection, donor-grouped gating, Fisher tests and exported summary tables. **Tags:** AG, SG, ICS, real-study.
- **Why inspect:** Explicitly constructs gating templates in R and applies them to study samples. Strong example for marker-specific thresholds and matched stimulation controls.
- **Caveat:** Requires local FCS/workspace inputs and custom gating methods; not a plug-and-play pipeline. **Paper:** *Optimization and Application of a T cell Antigen-Specific Activation Assay as Diagnostic and Treatment Monitoring Tool for Tuberculosis* (Musvosvi et al., 2025; title as supplied by repository).

### A02 — SoegaardLab: spectral-flow ICS and T-cell polyfunctionality (2026)

- **GitHub:** https://github.com/SoegaardLab/2026_Fisher_Garcia_Frattari_Naasz_et_al
- **Code:** [spectral_flow_ics](https://github.com/SoegaardLab/2026_Fisher_Garcia_Frattari_Naasz_et_al/tree/main/spectral_flow_ics); `2_pregate.R`, `3_qualitycontrol.Rmd`, `4_normalize.Rmd`, `5_gateMarkers_cytokines.R`, `6a`/`6b` clustering, `7_Polyfunctionality.Rmd`.
- **Experiment:** Spectral ICS with IFN, TNF, IL-2, CD107a and immune subset phenotyping. **Methods:** `CytoExploreR` drawn/reusable gating templates, `PeacoQC`, `cyCombine`, `CATALYST`/FlowSOM, boolean cytokine calls, polyfunctionality, unstimulated-background subtraction. **Tags:** SG, QC, BC, CL, ICS, spectral, real-study.
- **Why inspect:** One of the clearest complete spectral-flow workflows, from pregating to functional endpoints.
- **Caveat:** Marker gate CSV files are referenced but were not visible in the top-level `spectral_flow_ics` directory listing; manual construction/review remains integral. **Paper:** https://doi.org/10.1038/s41590-026-02448-z

### A03 — Fox et al.: Cyto-Feature Engineering, BCG vaccination / tuberculosis challenge (2020)

- **GitHub:** https://github.com/aef1004/cyto-feature_engineering
- **Code:** [BCG_v_PBS.Rmd](https://github.com/aef1004/cyto-feature_engineering/blob/master/writing/BCG_v_PBS.Rmd), [Human_clinical.Rmd](https://github.com/aef1004/cyto-feature_engineering/blob/master/writing/Human_clinical.Rmd).
- **Experiment:** Cytek Aurora mouse lung phenotyping following BCG or PBS vaccination and subsequent *M. tuberculosis* infection; second dataset is human whole blood. **Methods:** `openCyto` templates and `gt_gating`, initial live-lymphocyte gates, FMO-based positivity, marker-combination features and downstream association/visualization. **Tags:** AG, SG, spectral, real-study.
- **Why inspect:** Actual biological comparisons after hierarchical gating, including FMO-derived feature definitions.
- **Data:** Authors identify FlowRepository deposit in README; raw FCS not bundled on GitHub. **Paper:** Fox et al., *Scientific Reports* 10, 7651 (2020), https://doi.org/10.1038/s41598-020-64516-0

### A04 — Kaur Alasoo: flow cytometry genetics / macrophage cell-surface expression

- **GitHub:** https://github.com/kauralasoo/flow_cytomtery_genetics
- **Code:** [estimateSamplePurity.R](https://github.com/kauralasoo/flow_cytomtery_genetics/blob/master/analysis/preprocessing/estimateSamplePurity.R); `analysis/preprocessing/gating_parameters_2.csv`.
- **Experiment:** Human macrophage CD14, CD16 and CD206 surface expression, linked to genotype analyses. **Methods:** `gatingTemplate()` applied to transformed samples with a `GatingSet`; exclusion of previously identified failed samples; population purity estimation. **Tags:** AG, QC, real-study.
- **Why inspect:** Small, focused example of study-specific CSV templates and integrating cell purity into subsequent genomics analyses.
- **Caveat:** Data are referenced from processed local files; FCS availability not verified.

### A05 — Seth Masters Lab: ASC-speck flow assay and NLRP3 activity

- **GitHub:** https://github.com/Seth-Masters-Lab/Speck-Assay-Automated-Workflow
- **Code:** [Functions.R](https://github.com/Seth-Masters-Lab/Speck-Assay-Automated-Workflow/blob/main/Functions.R).
- **Experiment:** Detection of ASC specks and quantitative comparison of NLRP3 mutants. **Methods:** `openCyto::gate_flowclust_2d()` and `gate_mindensity()` for cell/gate selection, `flowAI`, FCS import and fitting EC50/ASC50-type summaries. **Tags:** AG, QC, real-assay.
- **Why inspect:** Function-based gates with a control-reference option and an assay-specific quantitative endpoint.
- **Caveat:** Not immune lineage phenotyping, and helper functions partly rely on legacy/internal interfaces.

### A06 — Davis McCarthy: fibroblast clonality / differentiation time course

- **GitHub:** https://github.com/davismcc/fibroblast-clonality
- **Code:** [facs_analysis_template.Rmd](https://github.com/davismcc/fibroblast-clonality/blob/master/src/Rmd/facs_analysis_template.Rmd).
- **Experiment:** FACS analysis of fibroblast differentiation at days 0–3. **Methods:** `openCyto::gatingTemplate()` using CSV files, `gating()` applied per day, `flowWorkspace`, `flowClean`, `ggcyto`; different gating templates where appropriate. **Tags:** AG, QC, longitudinal, real-study.
- **Why inspect:** Explicitly reuses gate hierarchies over timepoints and produces a documented report.
- **Caveat:** The analysis references `data_raw/facs/` inputs that were not present at the inspected GitHub path.

### A07 — Yachie Lab: CloneSelect / reporter-based FACS experiments

- **GitHub:** https://github.com/yachielab/CloneSelect_v1
- **Code example:** [C-to-T 6x6 assay](https://github.com/yachielab/CloneSelect_v1/blob/master/FACS/C-to-T_6x6_assay/process_fcs_files.R) and multiple `FACS/` notebooks.
- **Experiment:** HEK293 and mESC reporter/readout experiments. **Methods:** `CytoExploreR` for interactive FSC/SSC, singlet and GFP gates; saved CSV `gatingTemplate()` replayed by `gt_gating()` across samples; sample metadata and GFP-positive proportions. **Tags:** SG, real-study.
- **Why inspect:** Excellent applied example of drafting gates interactively, exporting a reusable template, and applying it systematically to all tubes/wells.
- **Caveat:** Reporter biology, not broad immune immunophenotyping.

### A08 — Chuong et al.: experimental evolution and GAP1 copy-number dynamics

- **GitHub:** https://github.com/chuongj/ee_gap1_archmuts ; newer study repository https://github.com/GreshamLab/local_arch_variants
- **Code example:** `lineplots_to_g240_031523.R`, `RawFlowPlots_Pub.R`, `Figure2A_pub.R`.
- **Experiment:** Yeast experimental evolution followed through multiple generations/timepoints. **Methods:** `CytoExploreR::cyto_gatingTemplate_apply()` across timepoints; gated frequencies, fluorescence measurements, time-series graphics. **Tags:** SG, longitudinal, real-study.
- **Data:** OSF and SRA links given in the final repository README. **Paper:** https://elifesciences.org/articles/98934
- **Why inspect:** Repeatable gating across a longitudinal, multi-file experiment.

### A09 — Bin He Lab: Pho4 transcription-factor evolution

- **GitHub:** https://github.com/binhe-lab/E013-Pho4-evolution
- **Code example:** `03-Pho4-chimera/03-flow-full-data-lfs/input/PHO5-data/20230216-gate-normalize-data.Rmd`.
- **Experiment:** Yeast fluorescent reporter expression after transcription-factor manipulation. **Methods:** `openCyto::gate_flowclust_2d()` for multivariate reporter/singlet gates, manual corrections, `ggcyto` reporting. **Tags:** AG, SG, real-study.
- **Why inspect:** Explicit troubleshooting and correction of poor automated gates, with repeated experimental measurements.

### A10 — Marine bacterial extracellular vesicles across the South Pacific

- **GitHub:** https://github.com/edfadeev/Marine_BEVs_across_South_Pacific
- **Code:** [R/process_FCM.R](https://github.com/edfadeev/Marine_BEVs_across_South_Pacific/blob/main/R/process_FCM.R).
- **Experiment:** Marine microbial cells and extracellular vesicle fractions. **Methods:** `flowWorkspace`, `flowAI`, `ggcyto`, `openCyto::gate_flowclust_2d()` and FCM-derived abundance processing. **Tags:** AG, QC, real-study.
- **Why inspect:** Real study uses data-driven gates for cell/particle measurement, distinct from mammalian immune phenotyping.
- **Paper:** Publication linked from repository README (ISME Journal); see https://github.com/edfadeev/Marine_BEVs_across_South_Pacific

### A11 — CYCLOPS: phenotyping within a cyclic imaging study

- **GitHub:** https://github.com/alkhalidis/CYCLOPS
- **Code:** [Cell LDA phenotyping / Gating.R](https://github.com/alkhalidis/CYCLOPS/blob/main/Cell%20LDA%20phenotyping/Gating.R).
- **Experiment:** Cell phenotyping support for a cyclic multiplex immunofluorescence workflow. **Methods:** `cyto_setup()`, `CytoExploreR::cyto_gate_draw()` and saved `gatingTemplate.csv`, `cyto_gatingTemplate_apply()`, downstream LDA scripts. **Tags:** SG, real-project.
- **Caveat:** The repository's primary purpose is imaging; the flow gating script is an ancillary analysis.

### A12 — den Braanker et al.: practical spectral-flow analysis (2021)

- **GitHub:** https://github.com/HdBraanker/Spectral_Flow_Workflow
- **Code:** [script.Rmd](https://github.com/HdBraanker/Spectral_Flow_Workflow/blob/main/script.Rmd).
- **Experiment:** Two 31-/36-color spectral flow case studies. **Methods:** `flowAI` and `PeacoQC` QC comparisons, transformation, CytoNorm batch normalization, FlowSOM/`CATALYST`, visual inspection and data quality diagnostics. **Tags:** SG, QC, BC, CL, spectral, public-FCS.
- **Data:** FlowRepository **FR-FCM-Z4KT**, **FR-FCM-Z3WR**. **Paper:** den Braanker, Bongenaar & Lubberts (2021), https://doi.org/10.3389/fimmu.2021.768113
- **Why inspect:** Full R Markdown script tied to real spectral FCS data, including alternative QC settings.

### A13 — Bonte et al.: 10 years of multicenter EuroFlow AML diagnostics (2026)

- **GitHub:** https://github.com/saeyslab/BonteEtAl_AML_Dx_EuroFlow
- **Code:** `Rscripts/1_DataPreprocessing.R`, `Rscripts/2_FlowSOM.R`, `Rscripts/3_BatchEffects.R`, later clinical modeling scripts.
- **Experiment:** Multicenter AML diagnostic immunophenotyping over 10 years; repository describes paper in *Cytometry Part B* (2026). **Methods:** FCS harmonization, compensation, `RemoveMargins`, `PeacoQC`, doublet removal, compensation-outlier statistics, `flowDensity::deGate`, FlowSOM reference models and CytoBatchFlagR-style batch diagnostics. **Tags:** QC, BC, CL, SG, real-study.
- **Why inspect:** Detailed metrics of how many events are excluded at each QC step and explicit cross-center compensation/batch checks.
- **Caveat:** Project-specific source file names and metadata paths; data access may be restricted.

### A14 — Saeys Lab: computational EuroFlow AML pipeline (earlier companion/research code)

- **GitHub:** https://github.com/saeyslab/Computational_pipeline_EuroFlow_AML
- **Code:** [0_preprocessing_all_files.R](https://github.com/saeyslab/Computational_pipeline_EuroFlow_AML/blob/main/0_preprocessing_all_files.R), `1_preprocessing_blasts.R`, `2_Complete_models_blasts.R`.
- **Experiment:** Real AML flow panels and clinical metadata. **Methods:** FCS and channel harmonization, `PeacoQC`, margin exclusion, `RemoveDoublets`, `flowDensity`, FlowSOM, CytoNorm, response/cell-state analysis. **Tags:** QC, BC, CL, real-study.
- **Why inspect:** A second full preprocessing implementation with applied per-file QC and normalization. Distinct repository from A13; may have overlapping clinical data and methods.

### A15 — Sodium chloride, CD8 T-cell function and tumor killing (*Nature Immunology*, 2024)

- **GitHub:** https://github.com/mBiocoder/Nature_immunology_code_repository
- **Code:** [CD8_Tm_FSFC](https://github.com/mBiocoder/Nature_immunology_code_repository/tree/main/CD8_Tm_FSFC), notably `Preprocessing_functions.R`, `Preprocess.R`, `Main.R`.
- **Experiment:** Full-spectrum flow cytometry of CD8 T-cell responses under altered NaCl conditions. **Methods:** Flexible `flowAI`/`PeacoQC` event QC, transformation, downstream high-dimensional phenotyping, UMAP, functions for reproducibility. **Tags:** QC, CL, spectral, real-study.
- **Data:** Scripts and dependent data linked via https://doi.org/10.5281/zenodo.12201093.
- **Why inspect:** Direct comparison of two automated QC algorithms applied to biological experiments.

### A16 — Lynall et al.: B cells and stress (2019 study code)

- **GitHub:** https://github.com/maryellenlynall/2019_bcell_stress
- **Code:** `bcellstress06.Rmd`, `bcellstress07.Rmd`, `bcellstress08a.Rmd`, `bcellstress08b.Rmd`.
- **Experiment:** Murine meningeal and splenic immune changes in stress paradigms. **Methods:** Read FlowJo workspaces, `flowAI::flow_auto_qc` event cleaning, FlowSOM, ConsensusClusterPlus, downstream population proportions and marker comparisons. **Tags:** SG, QC, CL, real-study.
- **Why inspect:** A transparent analysis that connects manual gates, automated QC, clustering and tests for biological changes.

### A17 — RGLab: openCyto paper analyses, HIV vaccine trials

- **GitHub:** https://github.com/RGLab/paper-opencyto
- **Code:** `pipeline-HVTN080.r`, `pipeline-HVTN065.r`, `gt-HVTN065.csv`, `generate-GatingTemplate-HVTN065.r`.
- **Experiment:** HVTN HIV vaccine trials; ICS gating and comparisons with manual analysis. **Methods:** `openCyto` templates, `flowWorkspace`, stimulation-specific analyses (negative/ENV/GAG), recovered manual gating sets and population statistics. **Tags:** AG, ICS, benchmark, real-study.
- **Why inspect:** Original application-level code with real vaccine-response gating and hierarchy configuration.
- **Caveat:** Older APIs and lab-local paths, so adaptation is required. **Paper:** Finak et al. (2014), https://doi.org/10.1371/journal.pcbi.1003806

### A18 — Ivison et al.: multicenter DuraClone immune phenotyping

- **GitHub:** https://github.com/mehrnoushmalek/DuraClone-gating
- **Code:** `Run_automated.R`, `OrganizeFiles.R` and panel-specific scripts.
- **Experiment:** Multicenter immune phenotyping, with Basic and B-cell DuraClone panels. **Methods:** `flowDensity` sequential supervised gates calibrated against manual hierarchies, comparison across sites and instrument settings. **Tags:** AG, SG, real-study, public-FCS.
- **Data:** FlowRepository **FR-FCM-ZYQT**. **Paper:** Ivison et al. (2018), https://doi.org/10.1172/jci.insight.121867
- **Why inspect:** Real multisite validation; tests automated versus manual reproducibility and details how instrument-specific adjustments were made.

### A19 — Rybakowska et al.: large-scale CyTOF immune monitoring (2021)

- **GitHub:** https://github.com/prybakowska/CyTOF_analysis_Pipeline1
- **Code:** `pipeline.R`, `functions.R`, `installation.R`.
- **Experiment:** Large-scale mass cytometry immune monitoring. **Methods:** `CATALYST` for normalization/debarcoding, `flowAI` and `flowCut` artifact cleaning, flowClean/flowDensity gating, CytoNorm, FlowSOM and UMAP. **Tags:** QC, AG, BC, CL, mass, real-study.
- **Data:** FlowRepository **FR-FCM-Z3YR** (reported in paper). **Paper:** Rybakowska et al. (2021), https://doi.org/10.1016/j.csbj.2021.05.032
- **Why inspect:** One of the most comprehensive published multibatch QC-and-normalization application pipelines.

### A20 — Immunosenescence cytometry analysis

- **GitHub:** https://github.com/masashi-CU/ISC-flow-cytometry
- **Code:** `immunosenescence.Rmd`.
- **Experiment:** Immunosenescence-associated cytometry analysis. **Methods:** `flowAI::flow_auto_qc`, QC metadata, FlowSOM metaclustering and visualizations. **Tags:** QC, CL, project-analysis.
- **Caveat:** Large exploratory R Markdown; associated publication and raw data not independently verified.

### A21 — NK manuscript: HIV-associated NK phenotyping

- **GitHub:** https://github.com/codeneeded/NK_Manuscript
- **Code:** `FLOWSOM/Scripts/QC_TARA.R`, `QC_FLORAH.R`.
- **Experiment:** NK-cell phenotyping in HIV-related experiments. **Methods:** FlowSOM, CytoNorm and cyCombine exploratory batch-effect evaluations; QC and cluster visualizations. **Tags:** QC, BC, CL, study-code.
- **Caveat:** Some `flowAI` QC calls are **commented out** in the inspected `QC_TARA.R`; do not count them as executed workflow steps without checking the specific script.

### A22 — Liechti et al.: twin-study cytometry analysis

- **GitHub:** https://github.com/saeyslab/Liechti_et_al_TwinStudy
- **Code:** `210708_loadDemographics.R`, `210708_originalFlowSOMmodels.R`, `210708_figures.R`.
- **Experiment:** Cytometry phenotyping in a twin study. **Methods:** Demographic metadata integration, FlowSOM model building, results and figure generation; code search also shows CytoNorm training. **Tags:** CL, BC, real-study.
- **Why inspect:** Applied FlowSOM and batch-normalization choices in a biological population-comparison setting.

### A23 — PIDgeon: automated cytometry for primary immunodeficiency

- **GitHub:** https://github.com/saeyslab/PIDgeon
- **Code:** `R_scripts/1. Healthy_controls/`, `R_scripts/3. Validation/`.
- **Experiment:** Clinical immune screening with training and three independent multicenter validation cohorts. **Methods:** Reference FlowSOM tree, population abundance features, disease-prediction models, interpretable patient reports; CytoNorm is used in analysis scripts. **Tags:** CL, BC, automated annotation, real-study.
- **Why inspect:** Reference-based cell assignment and reproducibility across independent sites rather than threshold-only gating.

### A24 — Pepper et al.: mucosal immune priming and influenza-related protection (2026)

- **GitHub:** https://github.com/jackspepper/Pepper-2026
- **Code:** `Analysis Pipeline.Rmd`, `ClusterID_App.R`, `data/`.
- **Experiment:** Nasal immune effects of *Haemophilus haemolyticus* inoculation and heterologous protection. **Methods:** Study-level R cytometry analysis, CytoNorm/FlowSOM-related workflows and interactive cluster support (inspect individual notebooks before assuming automation). **Tags:** CL, real-study.
- **Paper:** Pepper et al. (2026), https://doi.org/10.3389/fimmu.2026.1855464

### A25 — A scientific cytometry compendium applying automated gating to public OMIP data

- **GitHub:** https://github.com/wguesdon/everything-flow-cytometry
- **Code:** `scripts/01_automated_gating_pbmc.R`, `02_omip39_compensation_and_gating.R`, `06_omip43_asc_analysis.R`; `gating/` CSV templates.
- **Experiments:** OMIP-39 and OMIP-43 public flow repositories, among others. **Methods:** `openCyto` templates, instrument versus calculated compensation, FlowJo workspace comparison, FlowSOM and marker-positive population analysis. **Tags:** AG, QC, CL, reproducible public-case analyses.
- **Why inspect:** Reconstructs published manual gating hierarchies and compares them with automated results on the same FCS data. Not itself a single original biological-study repository.

### A26 — Additional historical application: Hally / personal R FCS project scripts

- **GitHub:** https://github.com/hally166/Cytometry-R-scripts
- **Code:** `BFP4Fiona.r`, `flowWorkspace.R`, `FCSCleaningComparisons/cleaningTest.R`.
- **Experiment:** Smaller laboratory FCS analyses and practical preprocessing comparisons. **Methods:** `flowAI`, `flowCut`, `openCyto`, `ggcyto`, gating sets and exported visualization. **Tags:** QC, SG, practical-project scripts.
- **Caveat:** Mixed demonstration and exploratory code; no single manuscript can be assigned to the entire repository.

---

## B. Standalone and adaptable multi-sample pipelines

This section intentionally includes pipeline/infrastructure repositories (some are package projects, unlike Section A).

### B01 — FredHutch/immune-flow-gating **[explicitly requested]**

- **GitHub:** https://github.com/FredHutch/immune-flow-gating
- **Code:** [main workflow](https://github.com/FredHutch/immune-flow-gating/blob/main/main.nf), [`bin/run.R`](https://github.com/FredHutch/immune-flow-gating/blob/main/bin/run.R), [`bin/process.R`](https://github.com/FredHutch/immune-flow-gating/blob/main/bin/process.R), [`bin/gate.R`](https://github.com/FredHutch/immune-flow-gating/blob/main/bin/gate.R).
- **Purpose:** Nextflow workflow automatically gates immune-cell populations from FCS files, adapting `HIPCCyto` processing steps into a configurable project workflow. **Methods:** FCS/channel standardization, `cytoqc`, spillover compensation, asinh transformation, `openCyto::gs_add_gating_method()`, `flowClust` lymphocyte gates and saved `GatingSet`s. **Tags:** AG, QC, Nextflow+R.
- **Key technical detail:** The `gate.R` file defines per-sample FSC density debris detection, batch-aware viability gating and FlowClust-derived lymphocyte ellipses.
- **Caveat:** Pre-gates rather than a complete lineage-to-rare-population immunophenotyping hierarchy. Some optional gates are commented out in `gate_gs()`.

### B02 — FredHutch/nf-flow-cytometry-qc **[explicitly requested]**

- **GitHub:** https://github.com/FredHutch/nf-flow-cytometry-qc
- **Code:** [main.nf](https://github.com/FredHutch/nf-flow-cytometry-qc/blob/main/main.nf), [`bin/clean.Rscript`](https://github.com/FredHutch/nf-flow-cytometry-qc/blob/main/bin/clean.Rscript), `bin/cluster.Rscript`, `bin/summarize.py`, `modules/`.
- **Purpose:** Sample-level and cohort-level QC workflow. **Methods:** Nextflow stages for sample-sheet parsing → FCS cleaning → unsupervised clustering → metrics extraction → summaries → HTML report. The cleaning script runs `flowClean::clean()` and `flowAI::flow_auto_qc()` with exception handling. **Tags:** QC, CL, Nextflow+R+Python.
- **Methods basis:** README cites Rahim et al., *Methods* (2018), https://doi.org/10.1016/j.ymeth.2017.12.015
- **Caveat:** A QC pipeline, not a deep immune-cell hierarchical gating solution. The cleaning step catches errors and may skip a failed algorithm; review those logs before treating outputs as fully QC-passed.

### B03 — DrCytometer/ImmunePhenotyping

- **GitHub:** https://github.com/DrCytometer/ImmunePhenotyping
- **Code:** `source/preprocess_flow_files.R`, `source/perform_batch_correction.R`.
- **Purpose:** Automated preprocessing and exploratory immune phenotyping. **Methods:** Optional spectral unmixing optimization, FlowJo-defined target population, PeacoQC, CytoNorm with internal reference controls, FlowSOM, UMAP, statistical tests, machine learning. **Tags:** SG, QC, BC, CL, spectral.
- **Caveat:** Requires a FlowJo workspace defining initial gated target populations; not end-to-end *fully automatic* gating.

### B04 — DavideCapozzi/pbmc-fcs-batch-correction

- **GitHub:** https://github.com/DavideCapozzi/pbmc-fcs-batch-correction
- **Code:** `src/functions/qc.R`, `main.R`, `config/`.
- **Purpose:** PBMC data normalization/batch-correction workflow with programmatic QC. **Methods:** PeacoQC → singlet/debris/viability gates → QC metric summaries → FlowSOM-related downstream analysis. Includes adaptive retry logic when a QC pass excludes unusually many events. **Tags:** AG, QC, BC, CL.
- **Caveat:** Active workflow engineering; no specific published study/data availability was verified.

### B05 — RGLab/HIPCCyto

- **GitHub:** https://github.com/RGLab/HIPCCyto
- **Purpose:** Standardization of HIPC ImmPort FCS files. **Methods:** harmonize channels, merge clinical/technical metadata, compensate, asinh transform, automatic pregating to lymphocytes, save gating sets. **Tags:** AG, QC, multi-study.
- **Example:** README uses ImmPort `SDY820` via `fetch_files()`/`process_study()`.
- **Caveat:** Package/pipeline infrastructure, and primarily pregating. The FredHutch/immune-flow-gating repository reproduces and modifies parts of this approach.

### B06 — ppreshant/flow_cytometry

- **GitHub:** https://github.com/ppreshant/flow_cytometry
- **Code:** `analyze_fcs.R`, `11-manual_gating_workflow.R` and Python FlowCal wrappers.
- **Purpose:** Laboratory FCS plate/well analysis; data cleanup and calibration followed by R gating and frequency visualization. **Methods:** Python `FlowCal` density and singlet cleaning, MEFL, R `openCyto` threshold gates, `flowWorkspace` counts and `ggcyto` plots. **Tags:** SG, QC, hybrid.
- **Caveat:** **Not full-R**. Author explicitly states R gating is not fully automated yet.

### B07 — 00berst33/flowFun (includes real-world analysis scripts)

- **GitHub:** https://github.com/00berst33/flowFun
- **Code:** [scripts/openCyto.R](https://github.com/00berst33/flowFun/blob/main/scripts/openCyto.R), `scripts/template.R`.
- **Purpose:** Experimental FCS preprocessing and reusable gating. **Methods:** import cytoset, compensate, logicle transform, programmatic template for non-margins/debris/singlets/live cells, `gt_gating()`, manually adjust/reapply gates, multi-sample plots. **Tags:** SG, AG, QC.
- **Caveat:** Code explicitly documents encountered errors and lab-local paths; consider it an applied working notebook, not a validated reproducible release.

### B08 — peterleary/MARMOT

- **GitHub:** https://github.com/peterleary/MARMOT
- **Code:** `inst/pipeline/MARMOT_Pipeline.qmd`.
- **Purpose:** Multifaceted R pipeline for spectral flow data. **Methods:** configurable `PeacoQC` or `flowAI` cleaning, dimensionality reduction, FlowSOM and alternative clustering, data integration and report generation. **Tags:** QC, CL, spectral, pipeline.
- **Caveat:** Framework/application rather than a manuscript-specific code deposit; confirm installation and example inputs.

### B09 — camellialambert/CARROT

- **GitHub:** https://github.com/camellialambert/CARROT
- **Code:** `pipeline_cytometrie.R`, `R/module_qc.R`.
- **Purpose:** Cytometry data analysis software with conventional/spectral modes. **Methods:** flowAI and PeacoQC per-sample QC, QC report archiving, FlowSOM clustering and plotting, interactive/manual gates. **Tags:** SG, QC, CL.
- **Caveat:** Tool/application development repository, not evidence of performance in a published study.

### B10 — paiva-s-lab/FlowCT

- **GitHub:** https://github.com/paiva-s-lab/FlowCT
- **Purpose:** Semi-automated pipeline for large immunophenotyping studies; standardized compensated inputs, simultaneous many-file analysis, clustering and biological statistics. **Methods:** `flowAI` sample cleaning, clustered phenotypes and summary exports. **Tags:** SG, QC, CL, clinical-scale.
- **Paper:** Botta et al. (2021), *Blood Advances*, https://doi.org/10.1182/bloodadvances.2021005198
- **Caveat:** This is a package/workflow rather than scripts of one specific clinical dataset.

### B11 — prybakowska/CytoQP

- **GitHub:** https://github.com/prybakowska/CytoQP
- **Purpose:** Cytometry Quality Pipeline for multibatch, large-scale blood studies. **Methods:** bead-based normalization, flow rate and signal cleaning (`flowAI`, `flowCut` or PeacoQC), debarcoding, viable/singlet/intact-cell selection, reference-sample CytoNorm, FlowSOM features, visualization. **Tags:** AG, QC, BC, CL, mass.
- **Caveat:** Packaged method, rather than a study analysis. Compare with the Rybakowska 2021 application repository (A19).

### B12 — UCLouvain-CBIO/CytoPipeline

- **GitHub:** https://github.com/UCLouvain-CBIO/CytoPipeline
- **Purpose:** R framework to explicitly configure, execute and visualize cytometry preprocessing/QC steps. **Methods:** S4 `CytoPipeline`/`CytoProcessingStep`, wrappers around `flowAI` and PeacoQC and other preprocessing routines. **Tags:** QC, infrastructure.
- **Caveat:** Reusable framework, not a particular biological study or deeply automated gating pipeline.

### B13 — lorenzobonaguro/cyCONDOR

- **GitHub:** https://github.com/lorenzobonaguro/cyCONDOR
- **Purpose:** High-dimensional cytometry analysis ecosystem. **Methods:** normalization including CytoNorm, clustering, differential abundance / analysis visualizations. **Tags:** BC, CL, package/workflow.
- **Caveat:** Evaluate as a broad interface/framework, not evidence of study-specific QC completeness.

---

## C. QC, normalization, and gating evaluations / benchmarks

### C01 — saeyslab/PeacoQC_evaluation

- **GitHub:** https://github.com/saeyslab/PeacoQC_evaluation
- **Code:** `R_Scripts/FlowCAP4/5_Quality_Control_Algorithms.R` and time/algorithm comparisons.
- **Purpose:** Analysis to evaluate PeacoQC relative to existing acquisition QC algorithms using real/reference datasets. **Methods:** PeacoQC, flowAI, flowCut, comparative exclusions/metrics. **Tags:** QC, benchmark.
- **Why inspect:** Allows looking beyond a package's vignette to algorithm-level validation and disagreement.

### C02 — 1ordinateur/FlowMOP_paper

- **GitHub:** https://github.com/1ordinateur/FlowMOP_paper
- **Code:** `benchmarks/run_default_biological_qc.R`, `benchmarks/parallel_peacoqc.sh`, `parallel_flowcut.sh`, `algo_validation/`.
- **Purpose:** QC benchmarking and synthetic-real mixing experiments for a FlowMOP paper. **Methods:** R implementations of PeacoQC/flowAI/flowCut as comparators, algorithm validation and figure generation. **Tags:** QC, benchmark, R+Python.
- **Caveat:** Benchmark repository, not an independent biological analysis pipeline; confirm datasets and publication status.

### C03 — hung-ching-chang/GatingMethod_evalutation

- **GitHub:** https://github.com/hung-ching-chang/GatingMethod_evalutation
- **Code:** `autogating.R`, `unsupervised_methods.R`, `python_code/`.
- **Purpose:** Quantitative comparison of manual, unsupervised and automated high-dimensional mass cytometry gating. **Methods:** R code for automated gating tools, comparison with 23 unsupervised methods; supplemental Python examples. **Tags:** AG, CL, benchmark, mass.
- **Caveat:** Mixed-language method evaluation rather than full FCS QC.

### C04 — openproblems-bio/task_cyto_batch_integration

- **GitHub:** https://github.com/openproblems-bio/task_cyto_batch_integration
- **Code:** `src/methods/cytonorm/script.R`, metric scripts, benchmark configuration.
- **Purpose:** Reproducible multi-method cytometry batch-integration benchmark. **Methods:** CytoNorm and other method adapters; metrics for batch removal versus biological abundance/label preservation; matched technical replicate designs. **Tags:** BC, benchmark.
- **Why inspect:** The most direct testing framework for whether normalization improves batch mixing without erasing cell-state/population structure.

### C05 — saeyslab/CytoNorm_Figures

- **GitHub:** https://github.com/saeyslab/CytoNorm_Figures
- **Code:** `2_normalisation_of_validation_samples.Rmd`, `4_evaluation_patient_samples.Rmd`.
- **Purpose:** Reproduce the analyses and figures from the CytoNorm publication. **Methods:** control-based quantile/spline normalization, comparisons with validation and patient data. **Tags:** BC, real-data method evaluation.
- **Why inspect:** Understands the assumptions and validation plots of reference-sample CytoNorm.

### C06 — RGLab/flowcap3

- **GitHub:** https://github.com/RGLab/flowcap3
- **Code:** `pipeline-DC.r`, `pipeline-bcell.r`, `pipeline-treg.r`, `pipeline-tcell.r`.
- **Purpose:** Historical FlowCAP-style automated gating benchmarks for different lineages. **Methods:** `openCyto` templates and flow population detection on standardized challenge data. **Tags:** AG, benchmark.
- **Caveat:** Old R infrastructure; valuable for historical algorithms and ground-truth comparison, not immediate re-execution.

### C07 — RGLab/faust_flowcap4_analysis

- **GitHub:** https://github.com/RGLab/faust_flowcap4_analysis
- **Purpose:** FlowCAP4 analysis code associated with the FAUST approach to cytometry population discovery. **Methods:** FlowCAP preprocessing and cross-method gating/annotation benchmarks; `flowDensity` referenced. **Tags:** AG, CL, benchmark.
- **Caveat:** Method-specific and historically dated; verify execution and input sources.

### C08 — AUMC-HEMA/CytoScan-manuscript

- **GitHub:** https://github.com/AUMC-HEMA/CytoScan-manuscript
- **Code example:** `02-Liechti_28-experiments.Rmd`.
- **Purpose:** Manuscript-associated analytical experiments applying FlowSOM to cytometry. **Tags:** CL, manuscript/benchmark.
- **Caveat:** Specific assay/domain and complete metadata context were not independently reviewed; inspect README before reuse.

### C09 — KyleeCJ/UNITO

- **GitHub:** https://github.com/KyleeCJ/UNITO
- **Purpose:** Bivariate-segmentation-based automated cytometric gating; publicly released data and expert consensus annotations. **Tags:** AG, benchmark, mixed-language.
- **Paper:** *Automated cytometric gating with human-level performance using bivariate segmentation*, *Nature Communications* (2025), https://www.nature.com/articles/s41467-025-56622-2
- **Data:** FlowRepository **FR-FCM-Z8XU** (mass), **FR-FCM-Z8P9** (flow), and Pennsieve linked from paper.
- **Caveat:** Primarily a new method rather than a full-R applied study; useful as a gating baseline.

---

## D. Core R cytometry tools and integration frameworks

These are **not** study-specific applications. They are retained because they can help reconstruct complete analysis pipelines or supply algorithm alternatives for later comparison.

| ID | Repository | What it offers | Key role / caveat |
|---|---|---|---|
| D01 | [RGLab/openCyto](https://github.com/RGLab/openCyto) | Hierarchical `gatingTemplate`, `gt_gating`, density/flowClust gates and user methods | Core automated gate engine; not a study |
| D02 | [RGLab/flowWorkspace](https://github.com/RGLab/flowWorkspace) | Persistent `GatingSet`s, sample metadata, transformations and population statistics | Reproducibility and persistence |
| D03 | [RGLab/CytoML](https://github.com/RGLab/CytoML) | FlowJo/Cytobank/GatingML workspace import/export | Manual-to-automated comparison/interoperability |
| D04 | [mehrnoushmalek/flowDensity](https://github.com/mehrnoushmalek/flowDensity) | Data-driven sequential gating using marker-density structures | Automated threshold-based lineage hierarchies |
| D05 | [DillonHammill/CytoExploreR](https://github.com/DillonHammill/CytoExploreR) | Interactive R gates, review/plots, templates and population summaries | Semi-automated gating; graphical intervention |
| D06 | [DillonHammill/CytoGate](https://github.com/DillonHammill/CytoGate) | Interactive gate drawing / openCyto plugin | Assisted gate editing |
| D07 | [saeyslab/PeacoQC](https://github.com/saeyslab/PeacoQC) | Peak-based time-course QC and event removal | Acquisition anomaly filtering; inspect plots and removed fractions |
| D08 | [giannimonaco/flowAI](https://github.com/giannimonaco/flowAI) | Flow rate, signal stability, dynamic range and acquisition QC | Alternative/complement to PeacoQC |
| D09 | [jmeskas/flowCut](https://github.com/jmeskas/flowCut) | Event-level anomaly exclusion | QC comparator and option in CyTOF workflows |
| D10 | [saeyslab/FlowSOM](https://github.com/saeyslab/FlowSOM) | Self-organizing maps and meta-clustering | High-dimensional cell phenotyping, not a direct gate hierarchy |
| D11 | [saeyslab/CytoNorm](https://github.com/saeyslab/CytoNorm) | Reference-control, FlowSOM-assisted quantile/spline normalization | Batch correction; requires suitable controls and assumptions |
| D12 | [biosurf/cyCombine](https://github.com/biosurf/cyCombine) | Batch alignment and high-dimensional cytometry integration | Alternative batch correction; biological conservation must be tested |
| D13 | [HelenaLC/CATALYST](https://github.com/HelenaLC/CATALYST) | CyTOF preprocessing, debarcoding, normalization, FlowSOM integration | Mass cytometry workflow backbone |
| D14 | [saeyslab/FlowSOM_protocol](https://github.com/saeyslab/FlowSOM_protocol) | Complete real-data FlowSOM analysis with QC/model checks and plots | *Nature Protocols* article-linked demonstrator |
| D15 | [prybakowska/CytoQP](https://github.com/prybakowska/CytoQP) | Multi-step mass cytometry QC and normalization | Full-stack reusable workflow (see B11) |
| D16 | [UCLouvain-CBIO/CytoPipeline](https://github.com/UCLouvain-CBIO/CytoPipeline) | Composable S4 cytometry QC processing steps | Pipeline architecture (see B12) |
| D17 | [paiva-s-lab/FlowCT](https://github.com/paiva-s-lab/FlowCT) | Semi-automated multi-file immunophenotyping and statistics | Clinical immunophenotyping (see B10) |
| D18 | [aef1004/cytotypr](https://github.com/aef1004/cytotypr) | FMO-based marker combinations and immune phenotype features | Linked to Fox et al. real-data study, but is itself a package |
| D19 | [lorenzobonaguro/cyCONDOR](https://github.com/lorenzobonaguro/cyCONDOR) | High-dimensional cytometry analysis ecosystem | BC + CL framework (see B13) |
| D20 | [DrCytometer/AutoSpectral](https://github.com/DrCytometer/AutoSpectral) | Automated/supervised support for spectral unmixing and QC | Spectral acquisition/unmixing-specific |
| D21 | [Olsman/CompensAID](https://github.com/Olsman/CompensAID) | Compensation calculation/diagnostics | Study compensation QC design |
| D22 | [saeyslab/iidx](https://github.com/saeyslab/iidx) | Differential expression and normalization/FlowSOM-feature workflows | Analysis beyond gating |
| D23 | [RGLab/HIPCCyto](https://github.com/RGLab/HIPCCyto) | ImmPort-focused study harmonization | Application framework (see B05) |
| D24 | [openproblems-bio/task_cyto_batch_integration](https://github.com/openproblems-bio/task_cyto_batch_integration) | Normalization benchmark methods and metrics | BC reference benchmark (see C04) |
| D25 | [saeyslab/FlowSOM_workshop](https://github.com/saeyslab/FlowSOM_workshop) | Practical FlowSOM examples and QC instruction | Training; not a manuscript data deposit |

**Core tool citations:** Finak et al. (2014) [openCyto](https://doi.org/10.1371/journal.pcbi.1003806); Quintelier et al. (2021) [FlowSOM protocol](https://doi.org/10.1038/s41596-021-00550-0); Rybakowska et al. (2021) [full mass cytometry QC](https://doi.org/10.1016/j.csbj.2021.05.032). For additional methods, obtain the appropriate citation from each repository's README/DESCRIPTION rather than inferring one.

---

## E. Tutorials and reproducible demonstrations with real datasets

These are lower-priority for **study-specific application** searches but provide hands-on patterns, images, benchmarks, and scripts.

| ID | Repository | What to inspect | What is and is not established |
|---|---|---|---|
| E01 | [UMGCCCFCSR/CytometryInR](https://github.com/UMGCCCFCSR/CytometryInR) | `course/08_WaysToGate`, `course/11_DataForStats` | OpenCyto template application with immune samples/cytokine statistics; instructional course |
| E02 | [hammerlab/t-cell-data](https://github.com/hammerlab/t-cell-data) | `Flow.Rmd` | FCS, `openCyto::gate_flowClust_2d`, gating template and customized plugin; Bookdown teaching/reproducibility resource |
| E03 | [RGLab/cytogallery](https://github.com/RGLab/cytogallery) | `content/en/gallery/opencyto/neg_control.Rmd` | Example control-sample gate and custom preprocessing; educational method gallery |
| E04 | [RGLab/openCyto](https://github.com/RGLab/openCyto) | `vignettes/HowToAutoGating.Rmd`, `HowToWriteCSVTemplate.Rmd` | Canonical template syntax, algorithm demonstrations; package vignette |
| E05 | [gfinak/OpenCytoTutorial](https://github.com/gfinak/OpenCytoTutorial) | `OpenCytoTutorial.Rmd` | Historical full gate hierarchy demonstration; instructional |
| E06 | [aef1004/cytotypr](https://github.com/aef1004/cytotypr) | `protocol/full_example.Rmd`, examples | Reuses published mouse data to describe FMO-based phenotype calling; tutorial/package |
| E07 | [saeyslab/ComputationalCytometry](https://github.com/saeyslab/ComputationalCytometry) | `ComputationalCytometry.Rmd` | Applied methods and normalisation/FlowSOM examples; training |
| E08 | [saeyslab/FlowSOM_protocol](https://github.com/saeyslab/FlowSOM_protocol) | `R/FlowSOM_protocol.R`, figure scripts | Real public FCS and publication-linked reproducibility pipeline, includes model QC |
| E09 | [OzetteTech/cytoverse-bioc-2023](https://github.com/OzetteTech/cytoverse-bioc-2023) | `vignettes/5_Gating_gatingTemplate.Rmd` | Reusable openCyto gate-hierarchy exercise; training workshop |
| E10 | [hally166/R_flowcytometry_course](https://github.com/hally166/R_flowcytometry_course) | `R handout.R`, `Your data.R` | flowCore/openCyto R course snippets; not a peer-reviewed application |

---

## F. Published applications with no confirmed public analysis repository

These **are studies**, but I could not verify a complete accessible study-code repository for each. Keeping them prevents losing useful validation evidence while screening GitHub projects.

| Study | Application / method | Paper |
|---|---|---|
| Hunter-Schlichting et al. (2020), *Health and Retirement Study* | Hybrid openCyto automated template + FlowAnnotator manual review for 24 immune subsets in 931 participants. Eighteen of 24 achieved correlation ≥0.80 against manual methods; rare populations were hardest. | https://doi.org/10.1038/s41598-020-65016-x |
| Lee et al. (2019), *High-Throughput Analysis of Clinical Flow Cytometry Data by Automated Gating* | Automated FMO-control-guided gating in large clinical T-cell panels (1,698 effector/memory and 1,908 Treg samples), benchmarked versus manual analysis. | https://doi.org/10.1177/1177932219838851 |
| Conrad et al. (2019), *Implementation and Validation of an Automated Flow Cytometry Analysis Pipeline for Human Immune Profiling* | Real immuno-oncology profiling; adoption of `flowDensity` for high-throughput immune phenotyping with manual validation. | https://doi.org/10.1002/cyto.a.23664 |
| Quintelier et al. (2026), NSCLC immune checkpoint blockade | Automated flow-cytometry analysis linked neutrophil:lymphocyte ratio changes with outcomes; uses QC, FlowSOM, normalization and density thresholds. Paper-linked code URL was not accessible during the preceding search, so do not treat as a verified public repository. | https://doi.org/10.1158/2326-6066.CIR-25-0662 |

---

## G. Study-review worksheet and screening criteria

### Screening table template (copy into a notes file or a spreadsheet later)

| Repository | Priority (later) | Data accessible? | End-to-end R? | Actual automatic gates? | Manual QC/override? | Reference controls? | Batch structure? | Gate validity checks? | Reproducible outputs? | Limitations/notes |
|---|---|---|---|---|---|---|---|---|---|---|
| _To be filled after review_ | — | — | — | — | — | — | — | — | — | — |

### What to record when opening each repository

1. **Inputs and provenance:** raw versus processed FCS; transformations, compensation/unmixing, channel map, stable sample and donor IDs, acquisition metadata.
2. **Gating design:** actual CSV `gatingTemplate` or coded gates; density, flowClust, 2D gating, reference/FMO gates; `collapseDataForGating` and `groupBy`; exactly which upstream nodes are manually drawn.
3. **QC auditability:** event/time artifacts, margin/saturation events, debris, doublets and viability; plot/report per file, number and proportion excluded, fallback thresholds and explicit failure flags.
4. **Batch and control logic:** reference sample available per batch, repeated control calibration, frozen vs refitted normalization, cross-batch validation, how donor and timepoint pairing is respected.
5. **Failure review:** empty or anomalously large gates, biologically implausible counts, rare-subset instability, unmixing residuals and autofluorescence in spectral FCS; whether a reviewer can override and record the reason.
6. **Biological validity:** comparison to expert gates and meaningful readouts, cell frequencies, MFIs, cytokine positivity, polyfunctionality, mixed-effects inference; whether population identities remain stable across conditions.
7. **Reproducibility:** package lockfile or container, unambiguous workflow entry point, tests, access to FCS/workspaces/templates, export of gate boundaries/indices and QC metrics.

### Limits of this catalog

- A GitHub repository mentioning `openCyto`, `CytoNorm`, `FlowSOM` or `PeacoQC` does **not** prove those methods were applied successfully in the published final analysis; code paths may be exploratory or commented out.
- GitHub often omits raw patient FCS data, FlowJo workspaces, compensation controls, or gating-template CSV files. Mark those as missing rather than assuming full replication.
- Some listed programs are package-development, demos or benchmarking resources. They are deliberately segregated from real study repositories.
- A workflow using Nextflow + R scripts is **not** a full-R workflow; an R Markdown notebook that manually draws gates is **not** fully automated.
- Versions of `openCyto`, `flowWorkspace`, `CytoML`, `CytoExploreR` and underlying flowCore/cytolib APIs have changed. Code from 2013–2021 will frequently require modernization.
- URLs were verified as public links by repository/site inspection or code search where possible; no code was executed, and data accessibility, completeness, security and quantitative performance were not independently tested.

