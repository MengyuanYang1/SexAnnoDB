# SexAnnoDB
![Image text](https://github.com/MengyuanYang1/SexAnnoDB/blob/main/figures/Figures1.png)
Here we built SexAnnoDB(A comprehensive sex annotation database in cancer) to provide a resource and reference for intensive functional annotations of sex effects in cancer to facilitate the development of personalized precision medicine. In this database, we analyzed different levels of molecules, including molecules of genome, epigenome, transcriptome, and proteome of TCGA across 27 cancer types. For each molecule, we performed multiple systematic and bioinformatic analyses, including male-specific differential molecules, female-specific differential molecules, and sex-biased molecules. What's more, we analyzed transcript factor-gene regulation network, RNA binding protein-exon skipping regulation network, ceRNA of female and male groups, diverse variant information (Sex-biased Expression Quantitative Trait Loci (eQTL), Sex-biased Splicing Quantitative Trait Loci (sQTL)), and methylation information (Sex-biased Expression Quantitative Trait Methylation(eQTM), Sex-biased Splicing Quantitative Trait Methylation (sQTM)). SexAnnoDB is the first database that systematically and intensively annotates the sex effect across differential biomolecules in cancer, which will be a unique resource for cancer and drug research communities to develop personalized precision medicine plans for patients.


# Usage:

        python SexAnnoDB.py

# Cancers
cancer_types  =  ['BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL']


#Sex-biased DNA mutation

      Rscript Mutation_Analysis.r

#Sex-biased DNA Methylation

#Sex-biased DNA methylation affects gene expression

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


