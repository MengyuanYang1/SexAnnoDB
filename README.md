# SexAnnoDB
![Image text](https://github.com/MengyuanYang1/SexAnnoDB/blob/main/figures/Figures1.png)
Here we built SexAnnoDB(A comprehensive sex annotation database in cancer) to provide a resource and reference for intensive functional annotations of sex effects in cancer to facilitate the development of personalized precision medicine. In this database, we analyzed different levels of molecules, including molecules of genome, epigenome, transcriptome, and proteome of TCGA across 27 cancer types. For each molecule, we performed multiple systematic and bioinformatic analyses, including male-specific differential molecules, female-specific differential molecules, and sex-biased molecules. What's more, we analyzed transcript factor-gene regulation network, RNA binding protein-exon skipping regulation network, ceRNA of female and male groups, diverse variant information (Sex-biased Expression Quantitative Trait Loci (eQTL), Sex-biased Splicing Quantitative Trait Loci (sQTL)), and methylation information (Sex-biased Expression Quantitative Trait Methylation(eQTM), Sex-biased Splicing Quantitative Trait Methylation (sQTM)). SexAnnoDB is the first database that systematically and intensively annotates the sex effect across differential biomolecules in cancer, which will be a unique resource for cancer and drug research communities to develop personalized precision medicine plans for patients.

# Web server

https://ccsm.uth.edu/SexAnnoDB

# Data 

TCGA muti-omics data	https://portal.gdc.cancer.gov/

RNA editing events	https://ccsm.uth.edu/CAeditome/

Alternative splicing	https://gdc.cancer.gov/about-data/publications/PanCanAtlas-Splicing-2018

# Software and algorithms

LncBase v3.0	https://diana.e-ce.uth.gr/lncbasev3

maftools 3.19	https://bioconductor.org/packages/release/bioc/html/maftools.html

GISTIC2.0	https://cloud.genepattern.org/gp/pages/index.jsf?

lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00125:6.15.30

GeoTcgaData	https://github.com/YuLab-SMU/GeoTcgaData

miRTarBase 9.0	https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php

miRDB 6.0	https://mirdb.org/

TargetScan huam 8.0	https://www.targetscan.org/vert_80/

CIS-BP 2.0	https://cisbp.ccbr.utoronto.ca/

StringDB v11.5	https://string-db.org/

CisBP-RNA 0.6	http://cisbp-rna.ccbr.utoronto.ca/

FIMO 5.5.3	https://meme-suite.org/meme/doc/fimo.html

netZooR 1.5.0	https://netzoo.github.io/netZooR/

MatrixeQTL 2.3	https://cran.r-project.org/web/packages/MatrixEQTL/index.html

DrugBank 5.0	https://go.drugbank.com/

DisGeNET v7.0	https://www.disgenet.org/

Enrichr	https://maayanlab.cloud/Enrichr/



# Usage:

        python SexAnnoDB.py

# Cancers
cancer_types  =  ['BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL']


#Sex-biased DNA mutation

      Rscript Mutation_Analysis.r

#Sex-biased DNA Methylation and Sex-biased DNA methylation affects gene expression

        Methyl_analysis.py 

#Sex-biased mRNA

        DiffGene_summary.py

#Sex-biased Exon skipping events
  
        AS_anlaysis.py

#Sex-biased RNA A-to-I editing

        Editing_anlaysis.py

#Sex-biased protein

        Proteindata_Ana.py

#Sex-biased sQTL

        Rscript sQTL.r

#Sex-biased eQTL

         Rscript eQTL.r

#Sex-biased eQTM

         Rscript eQTM.r

#Sex-biased sQTM

         Rscript sQTM.r
         
#Sex-biased RBP-ES network
        
        Rscript RBP-ES-network.r
        
#Sex-biased TF-Coding gene network

        Rscript TF-Gene-network.r


