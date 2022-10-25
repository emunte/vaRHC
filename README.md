# vaRHC
R package to automate, as far as possible, the variant classification process for hereditary cancer genes. 


## Introduction
Variant classification is a manual complex long process that combines information of distinct nature. An accurate classification is necessary to ensure proper genetic counselling and personalised risk estimation. 
In 2015, the American College of Molecular Genetics and Genomics (ACMG) together with the Association of Molecular Pathologists published generic guidelines to standardize and provide an objective framework to evaluate variant pathogenicity in Mendelian disease.
Later, specific guidelines have been published for some genes by collaborative groups. 

`vaRHC` has been developed to automate as much as possible the process of variant classification in hereditary cancer.
It follows gene-specific guidelines for ATM, CDH1, CHEK2, MLH1, MSH2, MSH6, PMS2, PTEN and TP53, and the updated general ACMG/AMP rules for other cancer susceptibility genes. The final classification is obtained according to Tavtigian’s natural scoring Bayesian-based metastructure (Tavtigian et al., 2020) but also considering some of the proposed CanVIG-UK incompatibilities.

The current version of the package is based on the GRCh37 assembly of the human genome and works for single substitutions, deletions and insertions up to 25 bp, intronic variants and 5’ or 3’-UTR variants 25 bp beyond the coding sequence. 

## How to use it
Documentation (vignette) is available at the vignette folder.

## Citation
vaRHC was developed by Elisabet Munté Roca. The paper is now being submitted. 

