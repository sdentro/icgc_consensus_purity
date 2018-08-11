# ICGC PCAWG-11 consensus purity pipeline

This repository contains the pipeline used to create the consensus purity estimates of the ICGC PCAWG data set. The procedure has been described in detail in https://www.biorxiv.org/content/early/2018/05/07/312041

## General procedure

The pipeline takes purity estimates from six individual callers (ABSOLUTE, Aceseq, Battenberg, cloneHD, JaBbA and Sclust) and four subclonal architecture callers (Ccube, CliP, CTPSingle and PhyloWGS) and combines them into a single consensus estimate by first removing outliers and subsequently taking the value that corresponds to the highest density across the estimates.

The runtime of this procedure is around 30 minutes.

## Dependencies

Software packages used to develop the code and run the pipeline on the PCAWG dataset. Installation of these packages should normally take a few minutes via Bioconductor.

```
R (version 3.1.0)
```

R libraries (all installed via Bioconductor)
```
Bioconductor (version 3.0)
BiocInstaller (version 1.16.5)
readr
ggplot2
gridExtra
reshape2
MASS
```

## Data bundle

This pipeline takes input data that has been bundled and is available TODO (access restricted). Beyond purity estimates from the 10 methods, the bundle contains these reference files:

| File | Description |
| --- | --- |
| consensus.20170119.purity.ploidy.annotated.txt | Older consensus table, used to obtain ploidy values |
| aberration_summary.txt | Summary information about CNAs used for figures |
| WGD-info_20170216.txt | Estimates of whole genome duplication per tumour, with timing |
| WGD-info_20170217.txt | Update of whole genome duplication estimates for selected tumours |

## How to run the pipeline

```
Rscript obtain_purity_estimates.R /path/to/data_bundle/ /path/to/output_dir
```

## Produced output
The pipeline produces a single file with the purity estimates.

| Column | Description |
| --- | --- |
| samplename | ICGC tumour aliquot id |
| purity | Estimated purity |
| ploidy | Calculated ploidy from the consensus copy number profile |
| purity_conf_mad | Confidence on the purity estimate |
| wgd_status | Contains `wgd` or `no_wgd` depending on whether the tumour is thought to have had a whole genome duplication |