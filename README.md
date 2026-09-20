**Distinct Genetic Models of Psychosis Exhibit Limited Neurodevelopmental Convergence in Telencephalic Organoids**

Deepak A. Kaji, Gabriel E. Hoffman, Panos Roussos

Overview

This repository contains the computational workflows used to construct and analyze a single-cell RNA-sequencing atlas of human telencephalic organoids modeling idiopathic schizophrenia and rare genetic forms of psychosis. The study integrates scRNA-seq data from 9 independent studies, comprising approximately 1.38 million cells from 100 libraries and 40 biological donors, spanning idiopathic schizophrenia and four rare genetic models of psychosis:

22q11.2 deletion syndrome
NRXN1 (2p16.3) deletion
15q13.3 deletion syndrome
3q29 deletion syndrome

The analyses include read alignment, quality control, ambient RNA removal, genotype-based demultiplexing, single-cell integration and annotation, donor-aware differential expression analysis, fixed-effects meta-analysis, cell-type compositional analysis, multivariate adaptive shrinkage (mash), and power/sensitivity analyses.This repository contains analysis code only. Raw sequencing data and processed data objects are not included.

Repository organization

The repository is organized according to the major stages of the analysis pipeline.

download_files/

Scripts used to obtain and prepare sequencing data from the contributing studies.

align_fastqs/

Scripts used to align raw sequencing reads using STARsolo against the GRCh38 reference genome and transcript annotation.

cellbender_filtration/

CellBender workflows used for ambient RNA correction and empty-droplet identification across datasets.

qc_datasets/

Dataset-specific quality-control workflows used after alignment and CellBender processing.

demuxing/

Genotype-based demultiplexing workflows, including SNP calling and donor assignment for datasets requiring genotype-based demultiplexing. This is relevant for the Shin et al dataset.

integrate_annotate/

Workflows used to construct the joint single-cell atlas, including:

scVI integration
hierarchical clustering
cell-type and subtype annotation
iterative subclustering
CellTypist query to reference mapping
differential expression analyses supporting cluster annotation
disease_analyses/

Primary statistical analyses of disease-associated transcriptional effects, including:

donor-aware differential expression using dreamlet
fixed-effects meta-analysis
cell-type compositional analysis using crumblr
variance partitioning
direct differential-expression overlap analyses
multivariate adaptive shrinkage (mashR)
mash sensitivity analyses

The primary analyses supporting the manuscript's Results and Figures 3–4 are contained within this directory.

Figures/

Jupyter notebooks used to generate the manuscript figures and supplementary figures.

Figure-generation notebooks are named according to the corresponding manuscript figure.

Abstracts_Prelim/

Preliminary analyses and exploratory work. These analyses are not required to reproduce the primary results presented in the manuscript.

network_analyses/

Exploratory transcriptional network analyses that are not part of the primary analyses presented in the manuscript.

trajectory_analyses/

Exploratory developmental trajectory analyses that are not part of the primary analyses presented in the manuscript.

Primary analysis workflow

The major computational workflow is:

Raw sequencing data
        ↓
STARsolo alignment
        ↓
CellBender ambient RNA correction
        ↓
Quality control
        ↓
Genotype-based demultiplexing
        ↓
scVI integration
        ↓
Hierarchical clustering and annotation
        ↓
Donor-aware differential expression
        ↓
Fixed-effects meta-analysis
        ↓
Disease-specific transcriptional signatures
        ↓
Direct overlap and pathway analyses
        ↓
mash shared-effect analysis
        ↓
Power and sensitivity analyses
        ↓
Manuscript figures and tables

Statistical framework

Differential expression analyses were performed using pseudobulk-based linear mixed models implemented in dreamlet, with biological donor and study-associated variables incorporated into the statistical framework.

Fixed-effects meta-analysis was used to derive consensus disease-associated transcriptional effects across independent manuscripts studying the same psychosis model.

Shared transcriptional effects across psychosis models were evaluated using multivariate adaptive shrinkage (mashR) and composite posterior testing.

Cell-type compositional differences were assessed using crumblr.

Power and sensitivity analyses were performed using simulations incorporating the empirically observed standard errors of gene-level effect estimates.

Reproducibility

The repository contains the scripts and notebooks used for the analyses reported in the manuscript.

Because the contributing studies were generated using different sequencing platforms, differentiation protocols, and data-access procedures, the complete raw datasets are not redistributed in this repository. Instead, scripts for dataset acquisition and preprocessing are provided where permitted.

Several analyses require large intermediate data objects and/or access-controlled source datasets. These files are therefore not included in the repository.

Where necessary, users should modify dataset-specific input and output paths to reflect their local computational environment.

Software

Major computational tools used in this study include:

STARsolo
CellBender
Scanpy
Pegasus
scVI/scANVI
CellTypist
dreamlet
crumblr
mashR
GSEApy
Enrichr
R
Python
Jupyter

Exact package versions and computational environment information are provided in the accompanying environment specification.

Figures and tables

The Figures/ directory contains notebooks corresponding to the manuscript figures and supplementary figures.

Figure 1	Figures/Figure_1_Subclass_Composition.ipynb
Figure 2	Figures/Figure_2_Subtype_Composition.ipynb
Figure 3	Figures/Figure_3-MetaAnalysis.ipynb
Figure 4	Figures/Figure_4-Shared_Divergent.ipynb
Supplementary Figure 1	Figures/Supplementary_Figure_1_Mixing.ipynb
Supplementary Figure 2	Figures/Supplementary_Figure_2_Atlas_Consensus.ipynb
Supplementary Figure 3	Figures/Supplementary_Figure_3_logFC.ipynb
Supplementary Figure 4	Figures/Supplementary_Figure_4_Deleted_Regions.ipynb
Supplementary Figures 5–6	Figures/Supplementary_Figure_5_6_Defining_Disease_Signatures.ipynb
Supplementary Figure 7	Figures/Supplementary_Figure_7_Cell_Type_Composition.ipynb
Supplementary Figure 8	Figures/Supplementary_Figure_8_mashR_Sensitivity_Analysis.ipynb
Supplementary Figure 9	Figures/Supplementary_Figure_9_mash_Power_Analysis.ipynb

Data availability

Raw sequencing data from the contributing studies are available through the repositories and accession numbers described in the manuscript and Supplementary Table 2, with the exception of the Khan et al. dataset, which was not publicly available because of privacy restrictions.

The integrated atlas and associated processed data will be made available through CELLxGENE.

Contact: For questions regarding the computational analyses, please contact Deepak Kaji (deepak.kaji@mountsinai.org).
