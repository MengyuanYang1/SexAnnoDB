import os,sys

cancer_types  =  ['BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL']


def ProteinAnanlysis():
    import Proteindata_Ana as PA
    #PA.ProceessingPro()
    PA.DiffProAnalysis()
def ASAnanlysis():
    import ASdata_anlaysis as AA
    AA.AS_diffalaysis()

def EditingAnalysis():
    import Editing_anlaysis as EA
    EA.Editing_diff()

def MehtyAnalysis(cancer_types):
    import Methyl_analysis as MA
    #MA.MergeMethylation(cancer_types)
    #MA.MethyAnalysis(cancer_types)
    #MA.sQTM_position()
    #MA.GeteQTMinput_mutiple(cancer_types)
    #MA.sQTM_cmd(cancer_types)
    #MA.Person_CpGandGene_mut(cancer_types)
    ###Promoter Analysis
    #MA.PromoterMethylationsite()
    #MA.GetsQTMinput_multip(cancer_types)
    #MA.sQTM_cmd2(cancer_types)
    #MA.Person_CpGandEX_mut(cancer_types)

def MutationAnalysis(cancer_types):
    import Mutation_analysis as MutA
    #MutA.GetMutation(cancer_types)
    #MutA.MafFileterPatients()
    #MutA.GetSexSpecificMutation()
    #MutA.SexbiasedMutationEnrichr()
    MutA.Mutation_Gene()

def QTLAnalysis(cancer_types):
    import QTL_analysis as SA
    SA.MakeSNPmatrix(cancer_types)
    #SA.SNP_matrix_idchange_multip()
    #SA.GeteQTLinput_multip(cancer_types)
    #SA.SNP_hg19tohg38_multip()
    #SA.Getposition()
    #SA.eQTL_cmd()
    #SA.EX_position()
    #SA.GetsQTLinput_multip(cancer_types)
    #SA.sQTL_cmd()

def Summary():
    import Gene_Summary as GS
    #GS.Gene_Summary()
    import DiffGene_summary as DS
    #DS.CodingGene_mut(cancer_types,'CodingGene')
    #DS.Coding_boxplot(cancer_types)
    #DS.PCheeatmap(cancer_types)
    #DS.CodingGene_mut(cancer_types,'lncRNA')
    #DS.CodingGene_mut(cancer_types,'miRNA')
    #DS.Summary_CanType(cancer_types)
    import Proteindata_Ana as PA
    #PA.Protein(cancer_types)
    #PA.Proteinheatmap(cancer_types)
    import ASdata_anlaysis as AA
    #AA.ASSummary(cancer_types,'ASevents')
    #AA.ASheeatmap(cancer_types)
    import Editing_anlaysis as EA
    #EA.Editing_preprocessing()
    #EA.Editing_diff()
    #print(len(cancer_types),'AtoIediting')
    #EA.Editing_Summary(cancer_types,'AtoIediting')
    #EA.Edheeatmap(cancer_types)
    import Network_Summary as NS
    #NS.RBP_EX()
    #NS.TF_Gene_Summary(cancer_types)
    NS.RBP_ES_Summary(cancer_types)
    #NS.ceRNA()
    #NS.eQTL_Summary(cancer_types)
    #NS.eQTL_boxplot(cancer_types)
    #NS.sQTL_Summary(cancer_types)
    #NS.sQTL_boxplot(cancer_types)
    import Methyl_analysis as MA
    #MA.PromoterME(cancer_types,'Methy')
    #MA.MEheeatmap(cancer_types)
    #MA.sQTMSummary_multip(cancer_types)
    #MA.sQTM_boxplot(cancer_types)
    #MA.eQTM_Summary(cancer_types)
    #MA.eQTM_plot(cancer_types)
    import Dis_Drug_Summary as DD
    #DD.Dis_info()
    import Cancer_Page as CP
    #CP.DiffSummary(cancer_types)
    #CP.Mergefile(cancer_types)
    import SexGene as SG
    #SG.GeneIDtoSym()

def main():
    #ProteinAnanlysis()
    #ASAnanlysis()
    #EditingAnalysis()
    #MehtyAnalysis(cancer_types)
    #MutationAnalysis(cancer_types)
    #QTLAnalysis(cancer_types)
    #Summary()

main()



        
