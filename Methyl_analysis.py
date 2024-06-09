import os,sys
import numpy as np
import math 
import pandas as pd
import Proteindata_Ana as PA
import ASdata_anlaysis as AA
import warnings
from scipy import stats
warnings.filterwarnings("ignore")
import multiprocessing



def GetMethySamInfo():
    f = open('/data2/myang9/SexAnnoDB/data/gdc_sample_sheet.2022-04-11_450k.tsv')
    d_SamInfo = {}
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        patient_id = line[6][0:15]
        cancer_type = line[4].split('-')[1]
        m_f = line[0]+'/'+line[1]
        if cancer_type in d_SamInfo:
            d_SamInfo[cancer_type].append((m_f,patient_id))
        else:
            d_SamInfo[cancer_type] = [(m_f,patient_id)]
    return d_SamInfo


def MergeMethylation(cancer_types):
    d_SamInfo = GetMethySamInfo()
    path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    for cancer_type in cancer_types:
        m_path = path + 'TCGA-'+cancer_type + '/harmonized/DNA_Methylation/Methylation_Beta_Value/'
        i = 0
        l_df = []
        print(cancer_type)
        for f in d_SamInfo[cancer_type]:
            m_f = m_path + f[0]
            patient_id = f[1]
            df = pd.read_table(m_f,sep='\t',header = None,index_col=0)
            df.columns=[patient_id]
            l_df.append(df)
            #df.to_csv(path+'TCGA-'+cancer_type+'-Methylation_matrix.txt',sep='\t', header=1, index=True)
        df_out = pd.concat(l_df,axis=1)
        df_out.to_csv(path+'TCGA-'+cancer_type+'-Methylation_matrix.txt',sep='\t', header=1, index=True)
def PromoterMethylationsite():
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    d_pro = {}
    r = open('/data2/myang9/SexAnnoDB/data/QTL_position/Gene_position.txt','w')
    head = ['geneid','chr','s1','s2']
    r.write('\t'.join(head)+'\n')
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 5:
            continue
        if line[2] != 'gene':
            continue
        if 'protein_codin' not in line[8]:
            continue
        chrom = line[0]
        geneid = line[8].split('"')[1].split('.')[0]
        strand = line[6]
        start = line[3]
        end = line[4]
        newl = [geneid,chrom,start,end]
        r.write('\t'.join(newl)+'\n')
        if strand == '+':
            s = line[3]
            promoter = (int(s)-2000,int(s)+2000)
        else:
            s = line[4]
            promoter = (int(s)-2000,int(s)+2000)
        if chrom in d_pro:
            d_pro[chrom][geneid]=promoter
        else:
            d_pro[chrom]={}
            d_pro[chrom][geneid]=promoter
    f = open('/data2/myang9/SexAnnoDB/data/DNA_Methylation_anno.txt') 
    f.readline()
    d_ME = {}
    i = 0
    for line in f:
        line = line.strip().split('\t')
        cgid = line[0]
        i = i + 1
        if i>1000:
            break
        chrom = line[2]
        if chrom[0:3] != 'chr':
            continue
        cg_pos = int(line[3])-1
        for geneid in d_pro[chrom]:
            if cg_pos >= d_pro[chrom][geneid][0] and cg_pos <=  d_pro[chrom][geneid][1]:
                if geneid not in d_ME:
                    d_ME[geneid] = []
                d_ME[geneid].append(cgid)
    r = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/Methy/CpgSiteinPromoter.txt','w')
    for geneid in d_ME:
        l = [geneid,','.join(d_ME[geneid])]
        r.write('\t'.join(l)+'\n')
def MethyAnalysis(cancer_types):
    d_Patientinfo = PA.PatientSubgroup()
    path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    print(path)
    l_task = []
    for cancer_type in cancer_types:
        task = (cancer_type,path,d_Patientinfo)
        l_task.append(task)
        if cancer_type != 'SARC':
            continue
        MethyAnalysis_step1(cancer_type,path,d_Patientinfo)
    #cores =  30
    #pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(MethyAnalysis_step1,l_task)
        
def MethyAnalysis_step1(cancer_type,path,d_Patientinfo):

        me_f = path + 'TCGA-{0}-Methylation_matrix.txt'.format(cancer_type)
        df = pd.read_csv(me_f,sep='\t',index_col=0)
        patient_ids =[x[0:15] for x in df.columns.values.tolist()]
        if cancer_type+'_turmol_female' not in d_Patientinfo:
            l_female_tumor = []
        else:
            l_female_tumor = [x for x in d_Patientinfo[cancer_type+'_turmol_female'] if x in patient_ids]
        if cancer_type+'_normal_female' not in d_Patientinfo:
            l_female_normal = []
        else:
            l_female_normal = [x for x in d_Patientinfo[cancer_type+'_normal_female'] if x in patient_ids ]
        if cancer_type+'_turmol_male' not in d_Patientinfo:
            l_male_tumor = []
        else:
            l_male_tumor =[x for x in  d_Patientinfo[cancer_type+'_turmol_male'] if x in patient_ids]
        if cancer_type+'_normal_male' not in d_Patientinfo:
            l_male_normal = []
        else:
            l_male_normal =[x for x in  d_Patientinfo[cancer_type+'_normal_male'] if x in patient_ids]
        print(cancer_type,len(l_female_tumor),len(l_female_normal),len(l_male_tumor),len(l_male_normal))
        df_female_tumor = df[l_female_tumor]
        df_female_normal = df[l_female_normal]
        df_male_tumor = df[l_male_tumor]
        df_male_normal = df[l_male_normal]
        outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/Methy/'
        df1 = df_female_tumor
        df2 = df_female_normal
        df_female_tumor = df_female_tumor.T
        df_female_tumor= df_female_tumor[['cg16998810','cg16626088','cg15727409','cg19736094','cg00179598','cg26039373','cg23750556','cg14759284','cg17407511','cg27485646','cg03509565']]
        df_female_tumor['group'] = 'F_T'
        df_male_tumor = df_male_tumor.T
        df_male_tumor = df_male_tumor[['cg16998810','cg16626088','cg15727409','cg19736094','cg00179598','cg26039373','cg23750556','cg14759284','cg17407511','cg27485646','cg03509565']]
        df_male_tumor['group'] = 'M_T'
        #print(df_female_tumor,df_male_tumor)
        print(l_male_tumor,l_female_tumor)
        d_out = pd.concat([df_female_tumor,df_male_tumor])
        print(d_out)
        d_out = d_out.T
        d_out.to_csv('/data2/myang9/SexAnnoDB/output/paper/Figure2_GPC3_CpGsites.txt',sep='\t')
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_male_tumor
        df2 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_female_tumor
        df2 = df_male_tumor
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_Tumor_femaleVSmale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_female_normal
        df2 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_Normal_femaleVSmale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
       

def sQTM_position():
    f = open('/data2/myang9/SexAnnoDB/data/DNA_Methylation_anno.txt')
    f.readline()
    r = open('/data2/myang9/SexAnnoDB/data/sQTM_SNP_Matrix_position.txt','w')
    head = ['snp','chr','pos']
    r.write('\t'.join(head)+'\n')
    for line in f:
        line = line.strip().split()
        cgid = line[0]
        pos =[line[2],str(int(line[3])-1)]
        newl = [cgid]+pos
        r.write('\t'.join(newl)+'\n')
    r.close()
    d_pro = {}
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    for line in f:
        line = line.strip().split('\t')
        if len(line)<5:
            continue
        if line[2] != 'gene':
            continue
        chrom = line[0]
        geneid = line[8].split('"')[1]
        strand = line[6]
        if strand == '+':
            s = line[3]
            promoter = (int(s)-2000,int(s)+2000)
        else:
            s = line[4]
            promoter = (int(s)-2000,int(s)+2000)
        promoter = [chrom]+list(promoter)
        promoter = [str(x) for x in promoter]
        d_pro[geneid]=promoter

    r = open('/data2/myang9/SexAnnoDB/data/sQTM_PSI_Matrix_position.txt','w')
    head = ['geneid','chr','s1','s2']
    r.write('\t'.join(head)+'\n')
    for geneid in d_pro:
        newl = [geneid]+d_pro[geneid]
        r.write('\t'.join(newl)+'\n')
    r.close()

def Pearson_cor(df1,df2):
    key = df1.columns.values.tolist()[0]
    l_cor = []
    for index2, col2 in df2.iteritems():
        l1 = df1[key].sort_index().values.tolist()
        l2 = col2.sort_index().values.tolist()
        newd = pd.DataFrame([l1,l2])
        newd = newd.T
        newd.columns = ['rbp','as']
        newd = newd.dropna()
        l1 = newd['as']
        l2 = newd['rbp']
        cor,p = stats.stats.pearsonr(l1,l2)
        l_cor.append(cor)
    return l_cor

def GetPromoterCpGsite():
    f = open('/data2/myang9/SexAnnoDB/output/Network/PromoterCpGsiteCor/CgGinPromoter_Summary.txt')
    d_cpg = {}
    for line in f:
        line = line.strip().split()
        geneid = line[0]
        cpgs = line[1].split(',')
        d_cpg[geneid] = cpgs 
    return d_cpg


def drop_col(df, cutoff=0.8):
    n = len(df)
    cnt = df.count()
    cnt = cnt / n
    return df.loc[:, cnt[cnt >= cutoff].index]


def coefficient_of_variation(data):
    mean = np.mean(data)
    std = np.std(data, ddof=0)
    cv = std/mean
    return cv

def Person_CpGandGene_mut(cancer_types):
    sexs = ['female','male']
    l_task = []
    for cancer_type in cancer_types:
        for sex in sexs:
            parameter = (cancer_type,sex)
            l_task.append(parameter)
            #print(cancer_type,sex)
            #Person_CpGandGene(cancer_type,sex)
    cores = 40
    pool = multiprocessing.Pool(processes=cores)
    pool.starmap(Person_CpGandGene,l_task)


def Person_CpGandGene(cancer_type,sex):
    #d_cpg = GetPromoterCpGsite()
    #geneids =  d_cpg.keys()
    outf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/correlation/{0}_{1}_tumor_Cor.txt'.format(cancer_type,sex)
    #f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    #pcgeneids = [line.strip().split()[1] for line in f]
    f_index = '/data2/myang9/SexAnnoDB/output/Network/eQTM/tmp/{0}_{1}_tumor_eQTM_cis.txt'.format(cancer_type,sex)
    if os.path.exists(f_index) == False:
        return 'Error'
    #if os.path.exists(outf):
    #    return 'Done'
    print(cancer_type,sex,'Run')
    print(outf)
    f = open(f_index)
    f.readline()
    d_index = {}
    for line in f:
        line = line.strip().split()
        cgid = line[0]
        FDR = float(line[-1])
        if FDR>0.05:
            continue
        geneid = line[1]
        if geneid in d_index:
            d_index[geneid].append(cgid)
        else:
            d_index[geneid] = [cgid]
    geneids = d_index.keys()
    #geneids =  [x for x in geneids if x in pcgeneids]
    exp_f = '/data2/myang9/SexAnnoDB/data/eQTM_input/tmp/{0}_{1}_tumor_eQTM_FPKM_scale.txt'.format(cancer_type,sex)
    cg_f = '/data2/myang9/SexAnnoDB/data/eQTM_input/tmp/{0}_{1}_tumor_eQTM_CpG.txt'.format(cancer_type,sex)
    df_gene = pd.read_csv(exp_f,sep='\t',index_col=0)
    df_cpg = pd.read_csv(cg_f,sep='\t',index_col=0)
    df_cpg = df_cpg.T
    df_cpg = drop_col(df_cpg)
    df_gene = df_gene.T
    #df_gene = np.log10(df_gene)# math.log10(df_gene)
    geneids = [gene for gene in  geneids if gene in df_gene.columns.values.tolist()]
    df_gene = df_gene[geneids]
    df_out = []
    #print(df_cpg,df_gene)
    r = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/correlation/{0}_{1}_tumor_Cor.txt'.format(cancer_type,sex),'w')
    head = ['CpG','Geneid','cor','p']
    r.write('\t'.join(head)+'\n')
    for geneid in geneids:
        l_exp = df_gene[geneid].values.tolist()
        for cpg in d_index[geneid]:
            if cpg not in df_cpg.columns.values.tolist():
                continue
            l_cpg = df_cpg[cpg]
            newd = pd.DataFrame([l_exp,l_cpg])
            newd = newd.T
            newd.columns = ['l1','l2']
            newd = newd.dropna()
            l1 = newd['l1']
            l2 = newd['l2']
            cor,p = stats.stats.pearsonr(l1,l2)
            l = [cpg,geneid,cor,p]
            l = [str(x) for x in l]
            r.write('\t'.join(l)+'\n')
    r.close()


def GencodeV22_SymToID():
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    d_SymToID = {}
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 5 :
            continue
        if line[2] != 'gene':
            continue
        #print(line[8].split('"'))
        geneid = line[8].split('"')[1]
        genesymbol = line[8].split('"')[7]
        d_SymToID[genesymbol] = geneid.split('.')[0]
    return d_SymToID

def GencodeV22_strucutre():
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    d_Stru = {}
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 5:
            continue
        if line[2] == 'transcript':
            continue
        geneid = line[8].split('"')[1].split('.')[0]
        Type = line[2]
        pos = line[0]+':'+'-'.join(line[3:5])
        if geneid in d_Stru:
            if Type in d_Stru[geneid]:
                d_Stru[geneid][Type].append(pos)
            else:
                d_Stru[geneid][Type]=[pos]
        else:
            d_Stru[geneid]={}
            d_Stru[geneid][Type] = [pos]
        #print(d_Stru[geneid])
        if line[2] != 'gene':
            continue
        strand = line[6]
        chrom = line[0]
        start = line[3]
        end = line[4]
        if strand == '+':
            s = line[3]
            promoter = (int(s)-2000,int(s)+2000)
            enhancer1 = (int(s)-2000,int(s)-5000)
            enhancer2 = (int(s)+2000,int(s)+5000)
        else:
            s = line[4]
            promoter = (int(s)-2000,int(s)+2000)
            enhancer1 = (int(s)-2000,int(s)-5000)
            enhancer2 = (int(s)+2000,int(s)+5000)
        d_Stru[geneid]['promoter'] = [line[0]+':'+str(promoter[0])+'-'+str(promoter[1])]
        d_Stru[geneid]['enhancer'] = [line[0]+':'+str(enhancer1[0])+'-'+str(enhancer1[1]),line[0]+':'+str(enhancer2[0])+'-'+str(enhancer2[1])]
    return d_Stru



def eQTM_Summary(cancer_types):
    cancer_types = ['THYM','LAML','PCPG','PAAD','SARC','THCA','ESCA','COAD','KIRP','KIRC','LUAD','STAD','LIHC','LUSC','LGG','SKCM','BLCA','HNSC']
    #d_cginfo = Methylation_anno()
    #d_Stru = GencodeV22_strucutre()
    path = '/data2/myang9/SexAnnoDB/output/Network/eQTM/tmp/'
    task_l = [] 
    for cancer_type in cancer_types:
        l = (cancer_type,'eQTM')
        #eQTM_summary1(cancer_type,'eQTM')
        task_l.append(l)
    cores = 20
    #pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(eQTM_summary1,task_l)
    #pool.starmap(eQTMSummary_1_select,task_l)
    #eQTM_Summary2(cancer_types)
    #eQTM_Summary3(cancer_types)
    eQTM_Summary4(cancer_types)

def eQTM_Summary4(cancer_types):
    inf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2.txt'
    f = open(inf)
    head = f.readline().strip().split()
    i = 1
    head = head[1:4]+head[5:]
    outf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2_Malespecific.txt'
    r1 = open(outf,'w')
    newhead = ['Index', 'Gene_ID', 'CpG_Site', 'CpG_Postion', 'Position_to_Gene', 'Effect_score', 'FDR', 'Cor.r', 'Cor.Pvalue', 'Cancer_type']
    r1.write('\t'.join(newhead)+'\n')
    outf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2_Femalespecific.txt'
    r2 = open(outf,'w')
    r2.write('\t'.join(newhead)+'\n')
    outf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2_opposite.txt'
    r3 = open(outf,'w')
    r3.write('\t'.join(head)+'\n')
    print(len(head),len(newhead))
    for line in f:
        line = line.strip().split()
        if line[6] == '-':
            wl =[str(i)] + line[1:4]+[line[5]]+line[8:10]+line[12:14]+[line[-1]]# ,'Female-specific']
            r2.write('\t'.join(wl)+'\n')
            #print(line,'\n',wl,'\n',len(wl),'\n\n')
        elif line[8] == '-':
            wl = [str(i)] + line[1:4]+line[5:8]+line[10:12]+[line[-1]]#,'Male-specific']
            #print(line,'\n',len(wl),wl,'\n\n')
            r1.write('\t'.join(wl)+'\n')
        else:
            #print(line)
            wl = line[1:4]+line[5:]
            #print(len(wl))
            r3.write('\t'.join(line)+'\n')
        i = i + 1
    

def eQTM_Summary3(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2.txt')
    head = f.readline()
    #print(head)
    d_opp = {}
    d_male = {}
    d_female = {}
    for line in f:
        line = line.strip().split()
        geneid = line[1]
        #genename = line[2]
        #print(line)
        cancer_type = line[-1]
        key = '|'.join([cancer_type,geneid])
        if '-'  not in line[6:10]:
            if key  in d_opp:
                d_opp[key].append(line[2])
            else:
                d_opp[key]= [line[2]]
        elif line[8] == '-':
            if key in d_male:
                d_male[key].append(line[2])
            else:
                d_male[key]=[line[2]]
        elif line[6] == '-':
            if key  in d_female:
                d_female[key].append(line[2])
            else:
                d_female[key]=[line[2]]
        else:
            print(line[10:14])
    r = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_CancerPage.txt','w')
    i = 0
    head = ['Index','Cancer_type','GeneID','CpGID','Type']
    r.write('\t'.join(head)+'\n')
    for key in d_female:
        i = i+1
        print(key,len(d_female[key]))
        if len(d_female[key])>10:
            d_female[key] = d_female[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_female[key]),'model2']
        r.write('\t'.join(wl)+'\n')
    for key in d_male:
        i = i+1
        print(key,len(d_male[key]))
        if len(d_male[key])>10:
            d_male[key] = d_male[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_male[key]),'model1']
        r.write('\t'.join(wl)+'\n')
    for key in d_opp:
        i = i+1
        print(key,len(d_opp[key]))
        if len(d_opp[key])>10:
            d_opp[key] = d_opp[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_opp[key]),'model3&4']
        r.write('\t'.join(wl)+'\n')
    r.close()

def eQTM_Summary2(cancer_types):
    d_cginfo = Methylation_anno()
    d_Stru = GencodeV22_strucutre()
    r = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary.txt','w')
    r2 = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2.txt','w')
    head = ['Index','Gene_ID','CpG_Site','CpG_Postion','CpG_Island','Position_to_Gene','Male.effect','Male.FDR','Female.effect','Female.FDR','Male.Cor','Male.Pvalue','Female.Cor','Female.Pvalue','Sex_biased','Cancer_type']
    print(head)
    r.write('\t'.join(head)+'\n')
    r2.write('\t'.join(head)+'\n')
    i = 0
    for cancer_type in cancer_types:
        sumf = '/data2/myang9/SexAnnoDB/output/Network/eQTM/{0}_eQTM_summary_select.txt'.format(cancer_type)
        print(sumf)
        f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/CodingGene/Summary/SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type))
        geneids = [line.strip().split()[0] for line in f]
        f.close()
        if os.path.exists(sumf):
            f = open(sumf)
            f.readline()
            for line in f:
                line = line.strip().split()
                cgid = line[1]
                cginfo = d_cginfo[cgid]
                pos = [x.split('_')[1] for x in cginfo[3].split('|') if x.split('_')[0]== line[0] ]
                pos = ','.join(pos)
                if len(pos)==0:
                    #print(pos)
                    pos = '-'
                if line[0] in geneids:
                    anno = 'Yes'
                else:
                    anno = '-'
                i = i+1
                line = [str(i)]+line[0:2]+cginfo[1:3]+[pos]+line[2:]+[anno,cancer_type]   
                r.write('\t'.join(line)+'\n')
                if pos== '-':
                    continue
                if line[7] != '-' and line[9]=='-':
                    if float(line[7])<0.0001 and line[11] != '-':
                        if abs(float(line[10]))>0.3 and float(line[11])<0.05:
                            r2.write('\t'.join(line)+'\n')
                elif line[7] =='-' and line[9] != '-':
                    if float(line[9])<0.0001 and line[13] != '-':
                        if abs(float(line[12]))>0.3 and float(line[13])<0.05:
                            r2.write('\t'.join(line)+'\n')
                else:
                    if '-' in line[10:14]:
                        print(line)
                        continue
                    if abs(float(line[10]))>0.3 and float(line[11])<0.05:
                        if abs(float(line[12]))>0.3 and float(line[13])<0.05:
                            r2.write('\t'.join(line)+'\n')
    r.close()

def eQTMSummary_1_select(cancer_type,Anatype):
    sumf = '/data2/myang9/SexAnnoDB/output/Network/{1}/{0}_sQTM_summary.txt'.format(cancer_type,Anatype)
    if os.path.exists(sumf)==False:
        return 'file not exist'
    fn = '/data2/myang9/SexAnnoDB/output/Network/{1}/tmp/{0}_female_tumor_{1}_cis_nosign.txt'.format(cancer_type,Anatype)
    mn = '/data2/myang9/SexAnnoDB/output/Network/{1}/tmp/{0}_male_tumor_{1}_cis_nosign.txt'.format(cancer_type,Anatype)
    f = open(fn)
    fn_l = [line.strip().split()[0]+'\t'+line.strip().split()[1] for line in f]
    f.close()
    f = open(mn)
    mn_l =  [line.strip().split()[0]+'\t'+line.strip().split()[1] for line in f]
    f.close()
    f = open(sumf)
    outf = '/data2/myang9/SexAnnoDB/output/Network/{1}/{0}_{1}_summary_select.txt'.format(cancer_type,Anatype)
    r = open(outf,'w')
    print(outf)
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        male_e = line[2]
        female_e = line[4]
        key = line[1] + '\t' + line[0]
        if male_e != '-' and female_e == '-':
            if key not in fn_l:
                continue
            r.write('\t'.join(line)+'\n')
        elif male_e == '-' and female_e !='-':
            if key not in mn_l:
                continue
            r.write('\t'.join(line)+'\n')
        else:
            r.write('\t'.join(line)+'\n')

def eQTM_summary1(cancer_type,Analysis):
    print(cancer_type)
    path = '/data2/myang9/SexAnnoDB/output/Network/eQTM/correlation/'
    cis_f = path +cancer_type+'_female_tumor_Cor.txt'
    cis_m = path +cancer_type+'_male_tumor_Cor.txt'
    if os.path.exists(cis_f) and os.path.exists(cis_m):
        f = open(cis_m)
        f.readline()
        d_out = {}
        for line in f:
            line = line.strip().split()
            key = line[0]+'_'+line[1]
            if line[2] == 'nan':
                continue
            d_out[key] = [line[2]+'p'+line[3],'-p-']
        f.close()
        f = open(cis_f)
        f.readline()
        for line in f:
            line = line.strip().split()
            key = line[0]+'_'+line[1]
            if line[2] == 'nan':
                continue
            if key in d_out:
                d_out[key][1]=line[2]+'p'+line[3]
            else:
                d_out[key] = ['-p-',line[2]+'p'+line[3]]
        f.close()
    path = '/data2/myang9/SexAnnoDB/output/Network/eQTM/tmp/'
    cis_f = path +cancer_type+'_female_tumor_eQTM_cis_select.txt'
    cis_m = path +cancer_type+'_male_tumor_eQTM_cis_select.txt'
    outdir =  '/data2/myang9/SexAnnoDB/output/Network/eQTM/'
    if os.path.exists(cis_f) and os.path.exists(cis_m):
        f = open(cis_m)
        f.readline()
        d_effect = {}
        for line in f:
            line = line.strip().split()
            cgid = line[0]
            geneid = line[1]
            key = cgid+'_'+geneid
            beta = line[2]
            d_effect[key] = [beta+'p'+line[-1],'-p-']
        f.close()
        f = open(cis_f)
        f.readline()
        for line in f:
            line = line.strip().split()
            cgid = line[0]
            geneid = line[1]
            beta = line[2]
            key = cgid + '_' + geneid
            if key not in d_effect:
                d_effect[key] = ['-p-',beta+'p'+line[-1]]
            else:
                d_effect[key][1] = beta+'p'+line[-1]
        f.close()
        r = open('{1}/{0}_sQTM_summary.txt'.format(cancer_type,outdir),'w')
        head = ['Gene ID','CpG Site','Male.effect','Male.FDR','Female.effect','Female.FDR','Male.Cor','Male.Pvalue','Female.Cor','Female.Pvalue']
        r.write('\t'.join(head)+'\n')
        for key in d_effect:
            if key not in d_out:
                continue
            wl =[key.split('_')[1],key.split('_')[0]]+d_effect[key][0].split('p')+d_effect[key][1].split('p') +d_out[key][0].split('p')+d_out[key][1].split('p')
            if '-' not in wl[2:6]:
                print(wl)
                if float(wl[2])>0 and float(wl[4])<0:
                    r.write('\t'.join(wl)+'\n')
                if float(wl[2])<0 and float(wl[4])>0:
                    r.write('\t'.join(wl)+'\n')
            else:
                r.write('\t'.join(wl)+'\n')



def eQTM_plot(cancer_types):
    #cancer_types = ['THYM','LAML','PCPG','PAAD','SARC','THCA','ESCA','COAD','KIRP','KIRC','LUAD','STAD','LIHC','LUSC','LGG','SKCM','BLCA','HNSC']
    eQTM_boxplot_makefile(cancer_types)
    eQTM_boxplot_Draw()

def eQTM_boxplot_makefile(cancer_types):
    import Proteindata_Ana as PA
    import QTL_analysis as SA
    Ana_type = 'eQTM'
    d_Patientinfo = PA.PatientSubgroup()
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2.txt')
    d_eQTM = {}
    f.readline()
    geneids = []
    cgids = []
    for line in f:
        line = line.strip().split()
        cancer_type = line[-1]
        info = line[1:]
        #print(line)
        geneid = line[1]
        cgid = line[2]
        geneids.append(geneid)
        cgids.append(cgid)
        if cancer_type in d_eQTM:
            d_eQTM[cancer_type].append(info)
        else:
            d_eQTM[cancer_type] = [info]
    geneids = list(set(geneids))
    cgids = list(set(cgids))
    print(len(geneids),len(cgids))
    path = '/data2/myang9/SexAnnoDB/data/eQTM_input/tmp/'
    i = 0
    for cancer_type in d_eQTM:
        d_group= SA.GetQTLinput_group(cancer_type,Ana_type,d_Patientinfo)
        snpf = '/data2/myang9/SexAnnoDB/data/GDCdata/TCGA-{0}-Methylation_matrix.txt'.format(cancer_type)
        expf = '/data2/myang9/SexAnnoDB/data/GDCdata/TCGA-{0}_expression_FPKM.csv'.format(cancer_type)
        df_snp = []
        df_sum = pd.read_csv(snpf,sep='\t',index_col=0,chunksize=200000)
        for df_tmp in df_sum:
            df_tmp = df_tmp.T
            df_tmp.columns = [rsid for rsid in df_tmp.columns.values.tolist()]
            df_tmp_rsids = [rsid for rsid in df_tmp.columns.values.tolist() if rsid in cgids]
            df_tmp = df_tmp[df_tmp_rsids]
            df_tmp = df_tmp.T
            df_snp.append(df_tmp)
        df_snp = pd.concat(df_snp)
        df_snp.columns = [x[0:12] for x in df_snp.columns.values.tolist()]
        df_snp = df_snp.loc[:,~df_snp.columns.duplicated()]
        df_snp = df_snp.T
        df_exp = pd.read_csv(expf,sep = ',',index_col=0)
        df_exp.columns = [x[0:12] for x in df_exp.columns.values.tolist()]
        df_exp = df_exp.loc[:,~df_exp.columns.duplicated()]
        df_exp = df_exp.T
        ESids_s = [x for x in df_exp.columns.values.tolist() if x in geneids]
        df_exp = df_exp[ESids_s]
        m_t = d_group['male_tumor']
        f_t = d_group['female_tumor']
        df_snp_f = df_snp.T[f_t].T
        df_snp_m = df_snp.T[m_t].T
        df_exp_f = df_exp.T[f_t].T
        df_exp_m = df_exp.T[m_t].T
        #print(df_exp_f,df_exp_m,df_snp_f,df_snp_m)
        for line in d_eQTM[cancer_type]:
            cgid = line[1]
            geneid = line[0]
            i = i+1
            snpf_f = df_snp_f[cgid]
            snpf_m = df_snp_m[cgid]
            psif_m = df_exp_m[geneid]
            psif_f = df_exp_f[geneid]
            f_tmp_f = pd.concat([snpf_f,psif_f],axis=1)
            f_tmp_m = pd.concat([snpf_m,psif_m],axis=1)
            f_tmp_f['Group'] ='Female'
            f_tmp_m['Group']='Male'
            f_tmp = pd.concat([f_tmp_m,f_tmp_f],ignore_index=True)
            if f_tmp.shape[1]>3:
                continue
            f_tmp.columns = ['beta_value','y','Group']
            f_tmp.loc[f_tmp['beta_value']<0.2,'x'] = 'Low'
            f_tmp.loc[(f_tmp['beta_value']>=0.2)&(f_tmp['beta_value']<0.6),'x'] = 'Middle'
            f_tmp.loc[f_tmp['beta_value']>=0.6,'x'] = 'High'
            outf = '/data2/myang9/SexAnnoDB/output/Figures/eQTM_boxplot/file/{0}-{1}_{2}-{3}-eQTM.txt'.format(geneid,cgid,geneid,cancer_type)
            f_tmp.index.name = 'Pids'
            f_tmp.dropna(axis=0,inplace=True)
            f_tmp.to_csv(outf,sep = '\t')


def eQTM_boxplot_Draw():
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTM/SexAnnoDB_eQTM_summary2.txt')
    f.readline()
    d_score = {}
    for line in f:
        line = line.strip().split('\t')
        cancer_type = line[-1]
        rsid = line[2]
        geneid = line[1]
        key = rsid + '_' + geneid+'_'+cancer_type
        if line[7] == '-':
            sign_m = 'ns'
        else:
            fdr_m = float(line[7])
            if fdr_m <0.0001:
                sign_m = '"****"'
            elif fdr_m <0.001:
                sign_m = '"***"'
            elif fdr_m <0.01:
                sign_m = '"**"'
            elif fdr_m <0.05:
                sign_m = '"*"'
            else:
                sign_m = 'ns'
        if line[9] == '-':
            sign_f = 'ns'
        else:
            fdr_f = float(line[9])
            if fdr_f <0.0001:
                sign_f = '"****"'
            elif fdr_f <0.001:
                sign_f = '"***"'
            elif fdr_f <0.01:
                sign_f = '"**"'
            elif fdr_f <0.05:
                sign_f = '"*"'
            else:
                sign_f = '"-"'
        l = [sign_f,sign_m]
        d_score[key] = l
    r = open('eQTM_Draw.sh','w')
    print('Please run nohup bash eQTM_Draw.sh &')
    i = 0
    path = '/data2/myang9/SexAnnoDB/output/Figures/eQTM_boxplot/file/'
    for inf in os.listdir(path):
        cancer_type = inf.split('-')[-2]
        key = inf.split('-')[1]+'_'+cancer_type
        if key not in d_score:
            print(inf)
        #print(key)
        sign_f = d_score[key][0]
        sign_m = d_score[key][1]
        inf = path + inf
        outf = inf.replace('file','figure').replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/boxplot_eQTM.r {0} {1} {2} {3} {4}'.format(inf,outf,sign_f,sign_m,cancer_type)
        i = i + 1
        r.write('nohup '+cmd+'  &\n')
        if i%10 ==0:
            r.write('wait;\n')
    os.system('nohup bash /data2/myang9/SexAnnoDB/output/eQTM_Draw.sh &')

def Methylation_anno():
    d_Stru = GencodeV22_strucutre() 
    d_SymToID = GencodeV22_SymToID()
    f = open('/data2/myang9/SexAnnoDB/data/DNA_Methylation_anno.txt')
    head = f.readline()
    d_cginfo = {}
    for line in f:
        line = line.strip().split()
        cgid = line[0]
        pos = ':'.join(line[2:4])
        genes = list(set(line[5].split(';')))
        geneids = [d_SymToID[gene] for gene in genes if gene in d_SymToID]
        geneids = list(set(geneids))
        Island = line[-2]
        position = line[-1]
        l = [cgid,pos,Island]
        cgpos = int(line[3])
        pos_l = []
        for geneid in geneids:
            gene_pos = []
            for gene_stru in d_Stru[geneid]:
                for pos in d_Stru[geneid][gene_stru]:
                    gene_stru_s = int(pos.split(':')[1].split('-')[0])
                    gene_stru_e = int(pos.split(':')[1].split('-')[1])
                    if cgpos > gene_stru_s and cgpos < gene_stru_e:
                        gene_pos.append( gene_stru)
            gene_pos = list(set(gene_pos))
            gene_pos = geneid + '_' + ','.join(gene_pos)
            pos_l.append(gene_pos)
            pos_l  = list(set(pos_l))
            if len(pos_l) == 0:
                pos_l = ['-']
                continue
        l = l + ['|'.join(pos_l)]
        d_cginfo[cgid] = l
    f.close()
    return d_cginfo        

def PromoterME(cancer_types,genekey):
    from statsmodels.stats.multitest import multipletests
    sum_path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    task_l = []
    for cancer_type in cancer_types:
        print(cancer_type)
        #PormoterME_TS(cancer_type,genekey)
        #PromoterME_FMS(cancer_type,genekey)
        l = (cancer_type,genekey)
        task_l.append(l)
        pass
    cores = 30
    #pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(PormoterME_TS,task_l)
    #pool.starmap(PromoterME_FMS,task_l)

    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    for t in keys:
        pass
        print(t)
        #PromoterME_sum1(sum_path,t,genekey)
        PromoterME_sum2(t,genekey)
    PromoterME_sum3(keys,genekey)

def PormoterME_TS(cancer_type,genekey):
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    f1 = path + '{0}_Tumor_femaleVSmale_diffResult.txt'.format(cancer_type)
    f2 = path + '{0}_Normal_femaleVSmale_diffResult.txt'.format(cancer_type)
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            ##Tumor specific sex biased gene
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            df1['dBeta'] = df1['ave1'] - df1['ave2']
            df2['dBeta'] = df2['ave1'] - df2['ave2']
            df1 = df1.loc[abs(df1['dBeta'])>0.1,]
            df2 = df2.loc[abs(df2['dBeta'])>0.1,]
            #print(df1.shape,df2.shape)
            tumor_ids = [x for x in df1._stat_axis.values.tolist() if x not in df2._stat_axis.values.tolist()]
            normal_ids = [x for x in df2._stat_axis.values.tolist() if x not in df1._stat_axis.values.tolist()]
            print(tumor_ids[0:10],normal_ids[0:10])
            df1 = df1.T
            df1 = df1[tumor_ids]
            df1 = df1.T
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
        else:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df1['dBeta'] = df1['ave1'] - df1['ave2']
            df1 = df1.loc[abs(df1['dBeta'])>0.1,]
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')

def PromoterME_FMS(cancer_type,genekey):
    print(cancer_type)
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    fdr_cut = 0.05
    f1 = path + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type)
    f2 = path + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            df1['dBeta'] = df1['ave1'] - df1['ave2']
            df2['dBeta'] = df2['ave1'] - df2['ave2']
            df1 = df1.loc[abs(df1['dBeta'])>0.1,]
            df2 = df2.loc[abs(df2['dBeta'])>0.1,]
            female_ids = [x for x in df1._stat_axis.values.tolist() if x not in df2._stat_axis.values.tolist()]
            male_ids = [x for x in df2._stat_axis.values.tolist() if x not in df1._stat_axis.values.tolist()]
            df1 = df1.T
            df1 = df1[female_ids]
            df1 = df1.T
            df2 = df2.T
            df2 = df2[male_ids]
            df2 = df2.T
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
            df2.to_csv(outf2,sep='\t')
        else:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df1['dBeta'] = df1['ave1'] - df1['ave2']
            df1 = df1.loc[abs(df1['dBeta'])>0.1,]
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
    else:
        if os.path.exists(f2) == True:
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df2 = df2.loc[df2['p_adj']<0.05,]
            df2['dBeta'] = df2['ave1'] - df2['ave2']
            df2 = df2.loc[abs(df2['dBeta'])>0.1,]
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df2.to_csv(outf2,sep='\t')

def PromoterME_sum1(sum_path,t,genekey):
    #d_cginfo = Methylation_anno()
    df_l = []
    for sumf in os.listdir(sum_path):
        if t not in sumf:
            continue
        print(sumf)
        cancer_type = sumf.split('_')[-1].split('.')[0]
        sumf = sum_path + sumf
        df = pd.read_csv(sumf,sep='\t') #,index_col=0)
        df['cancer_type'] = cancer_type
        df_l.append(df)
    df_out = pd.concat(df_l)
    df_out.index.name='Index'
    newindex = [str(x) for x in range(1,df_out.shape[0]+1)]
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    df_out.to_csv(outf,sep='\t')


def PromoterME_sum2(t,genekey):
    d_cginfo = Methylation_anno()
    sumf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    f = open(sumf)
    head = f.readline()
    r = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,genekey),'w')
    head = ['Index','Gene_ID','CpG_site','Position','CpG_Island','Position_to_Gene','ave1','ave2','wilcoxon.w','wilcoxon.p','p_adj','dBeta','cancer_type']
    r.write('\t'.join(head)+'\n')
    i = 0
    for line in f:
        line = line.strip().split()
        cgid = line[1]
        ave1 = float(line[2])
        ave2 = float(line[3])
        dvalue = ave1-ave2
        if abs(dvalue)<0.1:
            continue
        if cgid not in d_cginfo:
            continue
        infos = d_cginfo[cgid]
        if len(infos) < 4:
            continue
        for info in infos[-1].split('|'):
            if len(infos[-1])<1:
                continue
            i = i + 1
            gene_p = info.split('_')[1].replace('gene','gene body')
            l = [str(i)] + [info.split('_')[0]]+[line[1]]+ infos[1:3]+[gene_p] +line[2:]
            r.write('\t'.join(l)+'\n')
    f.close()
    r.close()

def PromoterME_sum3(keys,genekey):
    d_out = {}
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/SexAnnoDB_{0}_Summary3.txt'.format(genekey)
    for i in range(0,len(keys)):
        t = keys[i]
        sum2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,genekey)
        print(i,sum2)
        f = open(sum2)
        head = f.readline().strip().split('\t')
        for line in f:
            line = line.strip().split()
            key = '|'.join(line[1:6])
            cancer_type = line[-1]
            if key in d_out:
                if d_out[key][i] == '-':
                    d_out[key][i] = [cancer_type]
                else:
                    d_out[key][i].append(cancer_type)
            else:
                d_out[key] = ['-']*3
                d_out[key][i] = [cancer_type]
    i = 0
    r = open(outf,'w')
    newhead = head[0:6] + ['sex_biased', 'Male_specific', 'Female_specific']
    r.write('\t'.join(newhead)+'\n')
    for key in d_out:
        i = i + 1
        d_out[key] = [','.join(can) for can in d_out[key]]
        if d_out[key].count('-')==3:
            print(newl)
        newl =[str(i)]+ key.split('|')+ d_out[key]
        r.write('\t'.join(newl)+'\n')

def MEheeatmap(cancer_types):
    #MEheeatmap_Summaryfile(cancer_types)
    MEheeatmap_makefile(cancer_types)
    #MEheeatmap_Drawfile(cancer_types)


def MEheeatmap_Summaryfile(cancer_types):
    d_Patientinfo = PA.PatientSubgroup()
    f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/Methy/SexAnnoDB_Methy_Summary3.txt')
    f.readline()
    cgids = [line.strip().split('\t')[2] for line in f]
    f.close()
    df_out = []
    path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    colname = []
    for cancer_type in cancer_types:
        exp_f = path + 'TCGA-{0}-Methylation_matrix.txt'.format(cancer_type)
        df = pd.read_csv(exp_f,sep='\t',index_col=0)
        patient_ids =[x[0:15] for x in df.columns.values.tolist()]
        df.columns = patient_ids
        print(df)
        if cancer_type+'_turmol_female' not in d_Patientinfo:
            l_female_tumor = []
        else:
            l_female_tumor = [x for x in d_Patientinfo[cancer_type+'_turmol_female'] if x in patient_ids]
        if cancer_type+'_normal_female' not in d_Patientinfo:
            l_female_normal = []
        else:
            l_female_normal = [x for x in d_Patientinfo[cancer_type+'_normal_female'] if x in patient_ids ]
        if cancer_type+'_turmol_male' not in d_Patientinfo:
            l_male_tumor = []
        else:
            l_male_tumor =[x for x in  d_Patientinfo[cancer_type+'_turmol_male'] if x in patient_ids]
        if cancer_type+'_normal_male' not in d_Patientinfo:
            l_male_normal = []
        else:
            l_male_normal =[x for x in  d_Patientinfo[cancer_type+'_normal_male'] if x in patient_ids]
        print(cancer_type,len(l_female_tumor),len(l_female_normal),len(l_male_tumor),len(l_male_normal))
        df_female_tumor = df[l_female_tumor]
        df_female_normal = df[l_female_normal]
        df_male_tumor = df[l_male_tumor]
        df_male_normal = df[l_male_normal]
        l = [cancer_type+'_'+'female_tumor',cancer_type+'_'+'female_normal',cancer_type+'_'+'male_tumor',cancer_type+'_'+'male_normal']
        colname = colname + l
        df_l = [df_female_tumor,df_female_normal,df_male_tumor,df_male_normal]
        for df_tmp in df_l:
            df_tmp = df_tmp.T
            df_mean = df_tmp.mean(skipna=True)
            df_out.append(df_mean)
    result = pd.concat(df_out,axis=1)
    result.columns = colname
    outf = '/data2/myang9/SexAnnoDB/output/Figures/Methylation_PCheatmap_summary.txt'
    result.to_csv(outf)

def MEheeatmap_makefile(cancer_types):
    d_cginfo = Methylation_anno()
    d_cgids = {}
    for cgid in d_cginfo:
        for x in  d_cginfo[cgid][-1].split('|'):
            geneid = x.split('_')[0]

            if geneid in d_cgids:
                d_cgids[geneid].append(cgid)
            else:
                d_cgids[geneid] = [cgid]
    inf = '/data2/myang9/SexAnnoDB/output/Figures/Methyheatmap/Methylation_PCheatmap_summary.txt'
    df = pd.read_csv(inf,sep=',',index_col=0)
    
    df = df.T
    for geneid in d_cgids:
        #print(d_cgids[geneid])
        outf = '/data2/myang9/SexAnnoDB/output/Figures/Methyheatmap/file/{0}_heatmap.txt'.format(geneid)
        df_tmp = df[d_cgids[geneid]]
        df_tmp = df_tmp.T
        #print(df_tmp)
        df_tmp.to_csv(outf,sep=',')

def MEheeatmap_Drawfile(cancer_types):
    path = '/data2/myang9/SexAnnoDB/output/Figures/Methyheatmap/file/'
    r = open('/data2/myang9/SexAnnoDB/output/Figures/Methyheatmap/Methyheatmap_run.sh','w')
    i = 0
    for inf in os.listdir(path):
        inf = path + inf
        svgf = inf.replace('file','figure')
        svgf = svgf.replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/heatmap_Methy.r {0} {1}'.format(inf,svgf)
        r.write('nohup ' + cmd + ' &\n')
        i = i + 1
        if i%15 == 0:
            r.write('wait;\n')


def Person_CpGandEX_mut(cancer_types):
    sex_types = ['female','male']
    l_task = []
    for cancer_type in cancer_types:
        for sex in sex_types:
            l = (cancer_type,sex)
            l_task.append(l)
            #Person_CpGandEX(cancer_type,sex)
    cores = 40
    #print(l_task)
    pool = multiprocessing.Pool(processes=cores)
    pool.starmap(Person_CpGandEX,l_task)

def Person_CpGandEX(cancer_type,sex):
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/correlation/{0}_{1}_tumor_Cor.txt'.format(cancer_type,sex)
    f_index = '/data2/myang9/SexAnnoDB/output/Network/sQTM/tmp/{0}_{1}_tumor_sQTM_cis.txt'.format(cancer_type,sex)
    if os.path.exists(f_index) == False:
        return 'Error'
    #if os.path.exists(outf):
    #    return 'Done'
    print(cancer_type,sex)
    print(outf)
    f = open(f_index)
    f.readline()
    d_index ={}
    for line in f:
        line = line.strip().split()
        cgid = line[0]
        exid = line[1]
        if float(line[-1])>=0.05:
            continue
        if exid not in d_index:
            d_index[exid] = [cgid]
        else:
            d_index[exid].append(cgid)
    ex_f = '/data2/myang9/SexAnnoDB/data/sQTM_input/tmp/{0}_{1}_tumor_sQTM_PSI.txt'.format(cancer_type,sex)
    cg_f = '/data2/myang9/SexAnnoDB/data/sQTM_input/tmp/{0}_{1}_tumor_sQTM_CpG.txt'.format(cancer_type,sex)
    df_ex = pd.read_csv(ex_f,sep='\t',index_col=0)
    df_cpg = pd.read_csv(cg_f,sep='\t',index_col=0)
    df_cpg = df_cpg.T
    df_cpg = drop_col(df_cpg)
    df_ex = df_ex.T
    exids = d_index.keys()
    exids =[exid for exid in exids if exid in  df_ex.columns.values.tolist()]
    df_out = []
    r = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/correlation/{0}_{1}_tumor_Cor.txt'.format(cancer_type,sex),'w')
    head = ['CpG','exid','cor','p']
    r.write('\t'.join(head)+'\n')
    for exid in exids:
        l_ex = df_ex[exid].values.tolist()
        for cpg in d_index[exid]:
            if cpg not in df_cpg.columns.values.tolist():
                continue
            l_cpg = df_cpg[cpg].values.tolist()
            newd = pd.DataFrame([l_ex,l_cpg])
            newd = newd.T
            newd.columns = ['l1','l2']
            newd = newd.dropna()
            l1 = newd['l1']
            l2 = newd['l2']
            cor,p = stats.stats.pearsonr(l1,l2)
            l = [cpg,exid,cor,p]
            l = [str(x) for x in l]
            r.write('\t'.join(l)+'\n')
    r.close()

def gene_strand():
    d_strand = {}
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    for line in f:
        line = line.strip().split('\t')
        if len(line) <5:
            continue
        if line[2]!='gene':
            continue
        strand = line[6]
        if 'protein_coding' not in line[8]:
            continue
        geneid = line[8].split('"')[1].split('.')[0]
        d_strand[geneid] = strand
    return d_strand

def sQTM_Summary(cancer_types):
    #d_strand = gene_strand()
    #d_cginfo = Methylation_anno()
    #d_exidchange = AA.EXid_gene()
    task_l = []
    task_l2 = []
    for cancer_type in cancer_types:
        pass
        task2 = (cancer_type,'sQTM')
        task_l2.append(task2)
        #l = (cancer_type,d_cginfo,d_exidchange,d_strand)
        #task_l.append(l)
        #sQTMSummary_1_select(cancer_type,'sQTM')
        #sQTMSummary_1(cancer_type,d_cginfo,d_exidchange,d_strand)
    import multiprocessing
    cores = 20 
    pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(sQTMSummary_1,task_l)
    #pool.starmap(sQTMSummary_1_select,task_l2)
    #sQTMSummary_2(cancer_types)
    #sQTMSummary_3(cancer_types)
    #sQTMSummary_4(cancer_types)
    sQTMSummary_5(cancer_types)

def sQTMSummary_5(cancer_types):
    inf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2.txt'
    f = open(inf)
    head = f.readline().strip().split()
    head = head[0:4]+head[5:8]+head[9:]
    newhead = ['Index', 'Gene_ID', 'Gene_symbol', 'EX_ID', 'SkippedExon', 'CpG_Site', 'CpG_Postion', 'Position_to_EX', 'Effect_score', 'FDR','Cor.r', 'Cor.Pvalue', 'ORF_anno', 'Cancer_Type']
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2_Malespecific.txt'
    r2 = open(outf,'w')
    r2.write('\t'.join(newhead)+'\n')
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2_Femalespecific.txt'
    r1 = open(outf,'w')
    r1.write('\t'.join(newhead)+'\n')
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2_opposite.txt'
    r3 = open(outf,'w')
    r3.write('\t'.join(newhead)+'\n')
    print(len(head),len(newhead))
    for line in f:
        line = line.strip().split('\t')
        line[5] = line[4].split(':')[0] + ':'+line[5]
        if line[10] == '-':
            wl = line[0:4]+line[5:8]+line[9:10]+line[12:14]+line[16:]
            r1.write('\t'.join(wl)+'\n')
        elif line[12] == '-':
            #print(line)
            wl = line[0:4]+line[5:8]+[line[9]]+line[10:12]+line[14:16] + line[18:]
            r2.write('\t'.join(wl)+'\n')
            #print(len(wl),wl)
        else:
            print(line)
            wl = line[0:4] + line[5:8]+line[9:]
            r3.write('\t'.join(wl)+'\n')
            print(len(wl))
    r1.close()
    r2.close()
    r3.close()
        



def sQTMSummary_4(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2.txt')
    head = f.readline()
    d_male = {}
    d_female = {}
    d_opp = {}
    for line in f:
        line = line.strip().split('\t')
        geneid = line[1]
        genename = line[2]
        cancer_type = line[-1]
        exid = line[3]
        key = key = '|'.join([cancer_type,geneid,genename,exid])
        if '-'  not in line[10:14]:
            if key  in d_opp:
                d_opp[key].append(line[6])
            else:
                d_opp[key]= [line[6]]
        elif line[12] == '-':
            if key in d_male:
                d_male[key].append(line[6])
            else:
                d_male[key]=[line[6]]
        elif line[10] == '-':
            if key  in d_female:
                d_female[key].append(line[6])
            else:
                d_female[key]=[line[6]]
        else:
            print(line[10:14])
    r = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_CancerPage.txt','w')
    i = 0
    head = ['Index','Cancer_type','GeneID','GeneName','ESID','CpGID','Type']
    r.write('\t'.join(head)+'\n')
    for key in d_female:
        i = i+1
        print(key,len(d_female[key]))
        if len(d_female[key])>10:
            d_female[key] = d_female[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_female[key]),'model2']
        r.write('\t'.join(wl)+'\n')
    for key in d_male:
        i = i+1
        print(key,len(d_male[key]))
        if len(d_male[key])>10:
            d_male[key] = d_male[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_male[key]),'model1']
        r.write('\t'.join(wl)+'\n')
    for key in d_opp:
        i = i+1
        print(key,len(d_opp[key]))
        if len(d_opp[key])>10:
            d_opp[key] = d_opp[key][0:10]
        wl = [str(i)]+key.split('|')+[','.join(d_opp[key]),'model3&4']
        r.write('\t'.join(wl)+'\n')
    r.close()

def sQTMSummary_1(cancer_type,d_cginfo,d_exidchange,d_strand):
    path = '/data2/myang9/SexAnnoDB/output/Network/sQTM/correlation/'
    cis_f = path +cancer_type+'_female_tumor_Cor.txt'
    cis_m = path +cancer_type+'_male_tumor_Cor.txt'
    if os.path.exists(cis_f) and os.path.exists(cis_m):
        f = open(cis_m)
        f.readline()
        d_out = {}
        for line in f:
            line = line.strip().split()
            key = line[0]+':'+line[1]
            if 'nan' in line[2:]:
                continue
            d_out[key] = [line[2]+'p'+line[3],'-p-']
        f.close()
        f = open(cis_f)
        f.readline()
        for line in f:
            line = line.strip().split()
            key = line[0]+':'+line[1]
            if 'nan' in line[2:]:
                continue
            if key in d_out:
                d_out[key][1]=line[2]+'p'+line[3]
            else:
                d_out[key] = ['-p-',line[2]+'p'+line[3]]
        f.close()
    path = '/data2/myang9/SexAnnoDB/output/Network/sQTM/tmp/'
    cis_f = path +cancer_type+'_female_tumor_sQTM_cis_select.txt'
    cis_m = path +cancer_type+'_male_tumor_sQTM_cis_select.txt'
    outdir =  '/data2/myang9/SexAnnoDB/output/Network/sQTM/'

    if os.path.exists(cis_f) and os.path.exists(cis_m):
        f = open(cis_m)
        f.readline()
        d_effect = {}
        d_cg_to_ex = {}
        for line in f:
            line = line.strip().split()
            cgid = line[0]
            exid = line[1]
            key = cgid+':'+exid
            if key not in d_out:
                continue
            beta = line[2]
            d_effect[key] = [beta+'q'+line[-1],'-q-']
        f.close()
        f = open(cis_f)
        f.readline()
        for line in f:
            line = line.strip().split()
            cgid = line[0]
            exid = line[1]
            beta = line[2]
            key = cgid + ':' + exid
            if key not in d_out:
                continue
            if key not in d_effect:
                d_effect[key] = ['-q-',beta+'q'+line[-1]]
            else:
                d_effect[key][1] = beta+'q'+line[-1]
        f.close()
        r = open('{1}/{0}_sQTM_summary.txt'.format(cancer_type,outdir),'w')
        head = ['Gene_ID','Gene_symbol','ES_ID','ES_info','SkippedExon','CpG Site','CpG Postion','CpG Island','Position to ES events','Male.effect','Male.FDR','Female.effect','Female.FDR','Male.Cor','Male.Pvalue','Female.Cor','Female.Pvalue']
        r.write('\t'.join(head)+'\n')
        for key in d_effect:
            cgid = key.split(':')[0]
            exid = key.split(':')[1]
            cgpos = int(d_cginfo[cgid][1].split(':')[1])
            pos_l = []
            if exid not in d_exidchange:
                continue
            info = d_exidchange[exid]
            if info[0] not in d_strand:
                continue
            strand = d_strand[info[0]]
            chrom = info[3].split(':')[0]
            l = info[3].split(':')[1:]
            l.sort()
            l = [int(x) for x in l]
            s = int(l[2])
            e = int(l[3])
            if strand == '+':
                if cgpos<=l[0]:
                    pos_l.append('Distant upstream')
                elif  cgpos>=l[5]:
                    pos_l.append('Distant downstream')
                elif cgpos>=s-2 and cgpos<=s+2:
                    pos_l.append('Acceptor splice site')
                elif cgpos>=e-2 and cgpos<=e+2:
                    pos_l.append('Donor splice site')
                elif cgpos > l[0] and cgpos <s-2:
                    pos_l.append('Proximal upstream')
                elif cgpos > e+2 and cgpos <l[5] :
                    pos_l.append('Proximal downstream')
                elif cgpos >s+2 and cgpos <e-2:
                    pos_l.append('Skipped exon')
                else:
                    print('Error',cgpos,l,strand)
            else:
                if cgpos<=l[0]:
                    pos_l.append('Distant downstream')
                elif  cgpos>=l[5]:
                    pos_l.append('Distant upstream')
                elif cgpos>=s-2 and cgpos<=s+2:
                    pos_l.append('Donor splice site')
                elif cgpos>=e-2 and cgpos<=e+2:
                    pos_l.append('Acceptor splice site')
                elif cgpos > l[0] and cgpos <s-2:
                    pos_l.append('Proximal downstream')
                elif cgpos > e+2 and cgpos <l[5] :
                    pos_l.append('Proximal upstream')
                elif cgpos >s+2 and cgpos <e-2:
                    pos_l.append('Skipped exon')
                else:
                    print('Error',cgpos,l,strand)
            info = d_exidchange[exid]
            wl =info+d_cginfo[key.split(':')[0]][0:3]+pos_l+d_effect[key][0].split('q')+d_effect[key][1].split('q')+d_out[key][0].split('p')+d_out[key][1].split('p')
            if '-' not in wl[9:13]:
                if float(wl[9])>0 and float(wl[11])<0:
                    r.write('\t'.join(wl)+'\n')
                if float(wl[11])<0 and float(wl[9])>0:
                    r.write('\t'.join(wl)+'\n')
            else:
                r.write('\t'.join(wl)+'\n')
def sQTMSummary_1_select(cancer_type,Anatype):
    sumf = '/data2/myang9/SexAnnoDB/output/Network/{1}/{0}_{1}_summary.txt'.format(cancer_type,Anatype)
    if os.path.exists(sumf)==False:
        return 'file not exist'
    fn = '/data2/myang9/SexAnnoDB/output/Network/{1}/tmp/{0}_female_tumor_{1}_cis_nosign.txt'.format(cancer_type,Anatype)
    mn = '/data2/myang9/SexAnnoDB/output/Network/{1}/tmp/{0}_male_tumor_{1}_cis_nosign.txt'.format(cancer_type,Anatype)
    f = open(fn)
    fn_l = [line.strip().split()[0]+'\t'+line.strip().split()[1] for line in f]
    f.close()
    f = open(mn)
    mn_l =  [line.strip().split()[0]+'\t'+line.strip().split()[1] for line in f]
    f.close()
    f = open(sumf)
    outf = '/data2/myang9/SexAnnoDB/output/Network/{1}/{0}_{1}_summary_select.txt'.format(cancer_type,Anatype)
    r = open(outf,'w')
    print(outf)
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        male_e = line[10]
        female_e = line[12]
        key = line[5] + '\t' + line[2]
        if male_e != '-' and female_e == '-':
            if key not in fn_l:
                continue
            r.write('\t'.join(line)+'\n')
        elif male_e == '-' and female_e !='-':
            if key not in mn_l:
                continue
            r.write('\t'.join(line)+'\n')
        else:
            r.write('\t'.join(line)+'\n')

def sQTMSummary_2(cancer_types):
    import ASdata_anlaysis as AA
    d_ORF = AA.ORF_anno()
    r = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary.txt','w')
    head = ['Index','Gene_ID','Gene_symbol','EX_ID','EX_info','SkippedExon','CpG_Site','CpG_Postion','CpG_Island','Position_to_EX','Male.effect','Male.FDR','Female.effect','Female.FDR','Male.Cor','Male.Pvalue','Female.Cor','Female.Pvalue','ORF_anno','Cancer_Type']
    r.write('\t'.join(head)+'\n')
    i = 0
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/ASevents/Summary/SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type))
        exids = [exid.split()[0] for exid in f]
        f.close()
        if os.path.exists('/data2/myang9/SexAnnoDB/output/Network/sQTM/{0}_sQTM_summary_select.txt'.format(cancer_type))==False:
            continue
        f = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/{0}_sQTM_summary_select.txt'.format(cancer_type))
        head = f.readline().strip().split('\t')
        for line in f:
            i = i+1
            line = line.strip().split('\t')
            if 'cg' not in line[5]:
                continue
            orf_anno = d_ORF[line[2]]
            line = [str(i)] + line + [orf_anno,cancer_type]
            r.write('\t'.join(line)+'\n')
            '''
            if line[10] != '-' and line[14]!='-':
                if abs(float(line[14]))>0.3 and float(line[15])<0.05:
                    r.write('\t'.join(line)+'\n')
            if line[12] != '-' and line[16] != '-':
                if abs(float(line[16]))>0.3 and float(line[17])<0.05:
                    r.write('\t'.join(line)+'\n')
            #print(len(line),line)
            '''
def sQTMSummary_3(cancer_types):
    inf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary.txt'
    print(inf)
    f = open(inf)
    head = f.readline()
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2.txt'
    r = open(outf,'w')
    r.write(head)
    for line in f:
        line = line.strip().split('\t')
        if line[9] not in ['Donor splice site','Acceptor splice site','Skipped exon','Proximal upstream','Distant upstream']:
            continue
        if line[11] != '-' and line[13]=='-':
            if float(line[11])<0.0001 and line[14] != '-':
                if abs(float(line[14]))>0.3 and float(line[15])<0.05:
                    r.write('\t'.join(line)+'\n')
        elif line[11] == '-' and line[13] != '-':
            if float(line[13])<0.0001 and line[16]!= '-':
                if abs(float(line[16]))>0.3 and float(line[17])<0.05:
                    r.write('\t'.join(line)+'\n')
        else:
            if '-' in line[14:18]:
                print(line)
                continue
            #r.write('\t'.join(line)+'\n')
            if abs(float(line[14]))>0.3 and float(line[15])<0.05:
                if abs(float(line[16]))>0.3 and float(line[17])<0.05:
                    r.write('\t'.join(line)+'\n')
    



def sQTM_boxplot(cancer_types):
    sQTM_boxplot_makefile(cancer_types)
    sQTM_boxplot_Draw()


def sQTM_boxplot_makefile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2.txt')
    d_sQTM = {}
    f.readline()
    import Proteindata_Ana as PA
    import QTL_analysis as SA
    Ana_type = 'sQTM'
    d_Patientinfo = PA.PatientSubgroup()
    snpids = []
    ESids = []
    for line in f:
        line = line.strip().split()
        cancer_type = line[-1]
        ESids.append(line[3])
        snpids.append(line[6])
        info = line[1:]
        if cancer_type in d_sQTM:
            d_sQTM[cancer_type].append(info)
        else:
            d_sQTM[cancer_type] = [info]
    ESids = list(set(ESids))
    snpids = list(set(snpids))
    print(len(ESids),len(snpids))
    i = 0
    for cancer_type in d_sQTM:
        d_group= SA.GetQTLinput_group(cancer_type,Ana_type,d_Patientinfo)
        snpf = '/data2/myang9/SexAnnoDB/data/GDCdata/TCGA-{0}-Methylation_matrix.txt'.format(cancer_type)
        expf = '/data2/myang9/DriverRBP/data/ASevents/{0}_exon_skip_ASmatrix.txt'.format(cancer_type)
        df_snp = []
        df_sum = pd.read_csv(snpf,sep='\t',index_col=0,chunksize=200000)
        for df_tmp in df_sum:
            df_tmp = df_tmp.T
            df_tmp.columns = [rsid for rsid in df_tmp.columns.values.tolist()]
            df_tmp_rsids = [rsid for rsid in df_tmp.columns.values.tolist() if rsid in snpids]
            df_tmp = df_tmp[df_tmp_rsids]
            df_tmp = df_tmp.T
            df_snp.append(df_tmp)
        df_snp = pd.concat(df_snp)
        df_snp.columns = [x[0:12] for x in df_snp.columns.values.tolist()]
        df_snp = df_snp.loc[:,~df_snp.columns.duplicated()]
        df_snp = df_snp.T
        df_exp = pd.read_csv(expf,sep = '\t',index_col=0)
        df_exp.columns = [x[0:12] for x in df_exp.columns.values.tolist()]
        df_exp = df_exp.loc[:,~df_exp.columns.duplicated()]
        df_exp = df_exp.T
        ESids_s = [x for x in df_exp.columns.values.tolist() if x in ESids]
        df_exp = df_exp[ESids_s]
        m_t = d_group['male_tumor']
        f_t = d_group['female_tumor']
        df_snp_f = df_snp.T[f_t].T
        df_snp_m = df_snp.T[m_t].T
        df_exp_f = df_exp.T[f_t].T
        df_exp_m = df_exp.T[m_t].T
        #print(df_exp_f,df_exp_m,df_snp_f,df_snp_m)
        for line in d_sQTM[cancer_type]:
            rsid = line[5]
            exid = line[2]
            i = i+1
            snpf_f = df_snp_f[rsid]
            snpf_m = df_snp_m[rsid]
            psif_m = df_exp_m[exid]
            psif_f = df_exp_f[exid]
            f_tmp_f = pd.concat([snpf_f,psif_f],axis=1)
            f_tmp_m = pd.concat([snpf_m,psif_m],axis=1)
            f_tmp_f['Group'] ='Female'
            f_tmp_m['Group']='Male'
            f_tmp = pd.concat([f_tmp_m,f_tmp_f],ignore_index=True)
            if f_tmp.shape[1]>3:
                continue
            f_tmp.columns = ['beta_value','y','Group']
            f_tmp.loc[f_tmp['beta_value']<0.2,'x'] = 'Low'
            f_tmp.loc[(f_tmp['beta_value']>=0.2)&(f_tmp['beta_value']<0.6),'x'] = 'Middle'
            f_tmp.loc[f_tmp['beta_value']>=0.6,'x'] = 'High'
            #print(f_tmp)
            outf = '/data2/myang9/SexAnnoDB/output/Figures/sQTM_boxplot/file/{0}-{1}_{2}-{3}-sQTM.txt'.format(line[0],rsid,exid,cancer_type)
            #print(i,outf)
            f_tmp.index.name = 'Pids'
            f_tmp.dropna(axis=0,how='any',inplace=True)
            f_tmp.to_csv(outf,sep = '\t')

def sQTM_boxplot_Draw():
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTM/SexAnnoDB_sQTM_summary2.txt')
    f.readline()
    d_score = {}
    for line in f:
        line = line.strip().split('\t')
        rsid = line[6]
        exid = line[3]
        cancer_type = line[-1]
        key = rsid + '_' + exid +  '_' + cancer_type
        if line[11] == '-':
            sign_m = 'ns'
        else:
            fdr_m = float(line[11])
            if fdr_m <0.0001:
                sign_m = '"****"'
            elif fdr_m <0.001:
                sign_m = '"***"'
            elif fdr_m <0.01:
                sign_m = '"**"'
            else:
                sign_m = '"*"'
        if line[13] == '-':
            sign_f = 'ns'
        else:
            fdr_f = float(line[13])
            if fdr_f <0.0001:
                sign_f = '"****"'
            elif fdr_f <0.001:
                sign_f = '"***"'
            elif fdr_f <0.01:
                sign_f = '"**"'
            else:
                sign_f = '"*"'
        l = [sign_f,sign_m]
        d_score[key] = l
    r = open('sQTM_Draw.sh','w')
    i = 0
    path = '/data2/myang9/SexAnnoDB/output/Figures/sQTM_boxplot/file/'
    for inf in os.listdir(path):
            cancer_type = inf.split('-')[-2]
            key = inf.split('-')[1]+'_'+inf.split('-')[2]
            if key not in d_score:
                continue
            sign_f = d_score[key][0]
            sign_m = d_score[key][1]
            inf = path + inf
            outf = inf.replace('file','figure').replace('.txt','.svg')
            cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/boxplot_sQTM.r {0} {1} {2} {3} {4}'.format(inf,outf,sign_f,sign_m,cancer_type)
            i = i + 1
            r.write('nohup '+cmd+'  &\n')
            if i%10 ==0:
                r.write('wait;\n')
    print('Please run nohup bash /data2/myang9/SexAnnoDB/output/sQTM_Draw.sh &\n')
    os.system('nohup bash /data2/myang9/SexAnnoDB/output/sQTM_Draw.sh &')
