import os,sys
import DiffGene_summary as DS
import Proteindata_Ana as PA
import pandas as pd
from scipy import stats
import numpy as np
import math
from statsmodels.stats.multitest import multipletests


def ProceessingPro():
    datapath = '/data2/myang9/SexAnnoDB/data/Protein_Data/'
    f = open(datapath + 'gdc_sample_sheet.2022-03-22.tsv')
    d_files = {}
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        Cancer_type = line[4]
        filepath = datapath + 'OriginalData/' + line[0] +'/' + line[1]
        if Cancer_type in d_files:
            d_files[Cancer_type].append(filepath)
        else:
            d_files[Cancer_type] = [filepath]
    f.close()
    for Cancer_type in d_files:
        print(Cancer_type,len(d_files[Cancer_type]))
        r = open(datapath+'{0}_RPPA_expression.txt'.format(Cancer_type),'w')
        d_out = {}
        head = []
        for i in range(0,len(d_files[Cancer_type])):
            filepath = d_files[Cancer_type][i]
            patientid = filepath.split('/')[-1][0:15]
            head.append(patientid)
            f = open(filepath)
            f.readline()
            for line in f:
                line = line.strip().split()
                proid = line[4]
                proexp = line[5]
                if proid not in d_out:
                    d_out[proid] = ['NA']*len(d_files[Cancer_type])
                    d_out[proid][i] = proexp
                else:
                    d_out[proid][i] = proexp
            f.close()
            #print(patientid,i,proid,proexp)
        r.write('\t'.join(['Pro_ID']+head)+'\n')
        for proid in d_out:
            #print(proid,d_out[proid])
            l_w = [proid] + d_out[proid]
            r.write('\t'.join(l_w)+'\n')

def PatientSubgroup():
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/cancer_type_info_27.txt')
    d_Patientinfo = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        cancer_type = line[3]
        sex_type = line[1]
        tumortype = line[2]
        key = cancer_type + '_' + tumortype +'_' + sex_type
        if key in d_Patientinfo:
            d_Patientinfo[key].append(line[0][0:15])
        else:
            d_Patientinfo[key] = [line[0][0:15]]
    return d_Patientinfo
def differential_Wilcoxontest(df1,df2):
    #print(df1,df2)
    d_out = []
    rownames = []
    proids = df1._stat_axis.values.tolist()
    for proid in proids:
        l1 = df1.loc[proid].tolist()
        l2 = df2.loc[proid].tolist()
        l1 = [x for x in l1 if str(x) != 'nan']
        l2 = [x for x in l2 if str(x) != 'nan']
        if len(l1) <5 or len(l2) < 5:
            continue
        #print(l1,l2)
        w,p = stats.ranksums(l1,l2)
        ave1 = np.nanmean(l1)
        ave2 = np.nanmean(l2)
        #print(ave1,ave2)
        #log2fc= math.log(2,((ave1+1)/(ave2+1)))
        newl=[ave1,ave2,w,p]
        d_out.append(newl)
        rownames.append(proid)
    col_name = ['ave1','ave2','wilcoxon.w','wilcoxon.p']
    df = pd.DataFrame(d_out,columns=col_name,index=rownames)
    df = df.loc[df['wilcoxon.p']<0.05]
    df = df.dropna()
    p_adj = multipletests(df['wilcoxon.p'], method='fdr_bh')
    df['p_adj'] = p_adj[1]
    df = df.loc[df['p_adj']<0.05]
    #print(df)
    return df


def DiffProAnalysis():
    d_Patientinfo = PatientSubgroup()
    cancer_types = set([x.split('_')[0] for x in d_Patientinfo.keys()])
    for cancer_type in cancer_types:
        pro_f = '/data2/myang9/SexAnnoDB/data/TCGA/Protein_Data/TCGA-{0}_RPPA_expression.txt'.format(cancer_type)
        if os.path.exists(pro_f) == False:
            continue
        df = pd.read_csv(pro_f,sep='\t',index_col=0)
        #head = f.readline().strip().split()
        patient_ids = df.columns.values.tolist()
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

        df_female_tumor = df[l_female_tumor]
        df_female_normal = df[l_female_normal]
        df_male_tumor = df[l_male_tumor]
        df_male_normal = df[l_male_normal]
        print(cancer_type,len(l_female_tumor),len(l_female_normal),len(l_male_tumor),len(l_male_normal))
        #print(df_female_tumor,df_female_normal,df_male_tumor,df_male_normal)
        #print('\n\n\n')
        outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/Protein/'
          
        #-----group1-------
        df1 = df_female_tumor
        df2 = df_female_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_male_tumor
        df2 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_female_tumor
        df2 = df_male_tumor
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_Tumor_femaleVSmale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df1 = df_female_normal
        df2 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = differential_Wilcoxontest(df1,df2)
            outf = outdir + '{0}_Normal_femaleVSmale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
       
def Proteinheatmap(cancer_types):
    d_Patientinfo = PatientSubgroup()
    df_out = []
    colname = []
    path = '/data2/myang9/SexAnnoDB/data/TCGA/Protein_Data/'
    for cancer_type in cancer_types:
        if cancer_type=='LAML':
            continue
        exp_pro = path + 'TCGA-{0}_RPPA_expression.txt'.format(cancer_type)
        df = pd.read_csv(exp_pro,sep='\t',index_col=0)
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
        l = [cancer_type+'_'+'female_tumor',cancer_type+'_'+'female_normal',cancer_type+'_'+'male_tumor',cancer_type+'_'+'male_normal']
        colname = colname + l
        df_l = [df_female_tumor,df_female_normal,df_male_tumor,df_male_normal]
        for df_tmp in df_l:
            df_tmp = df_tmp.T
            print(df_tmp)
            df_mean = df_tmp.mean(skipna=False)
            print(df_mean)
            df_out.append(df_mean)
        result = pd.concat(df_out,axis=1)
        result.columns = colname
        outf = '/data2/myang9/SexAnnoDB/output/Figures/Proheatmap/Proheatmap_summary.txt'
        result.to_csv(outf)

def Protein(cancer_types):
    d_gene = DS.PCIDchange()
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/Protein/'
    d_sum = {}
    for cancer_type in cancer_types:
        diff_f = path + '{0}_Tumor_femaleVSmale_diffResult.txt'.format(cancer_type)
        if os.path.exists(diff_f)==False:
            continue
        f = open(diff_f)
        f.readline()
        for line in f:
            line = line.strip().split()
            genename = line[0]
            padj = float(line[-1])
            if padj <0.0001:
                key = cancer_type+'_'+'****'
            elif padj <0.001:
                key = cancer_type+'_'+'***'
            elif padj <0.01:
                key = cancer_type+'_'+'**'
            elif padj <0.05:
                key = cancer_type+'_'+'*'
            else:
                key = cancer_type+'_'+'ns'
            if genename in d_sum:
                d_sum[genename].append(key)
            else:
                d_sum[genename]= [key]
        f.close()
    i = 0
    r = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/Protein/SexAnnoDB_Tumor_femaleVSmale_Summary.txt','w')
    head = ['Index','Gene ID','Protein','Gene symbol','Cancer_type','Sign']
    r.write('\t'.join(head)+'\n')
    for genename in d_sum:
        cancer_types = [x.split('_')[0] for x in d_sum[genename]]
        signs = [x.split('_')[1] for x in d_sum[genename]]
        if genename not in d_gene:
            if genename.split('_')[0] in d_gene:
                geneid = d_gene[genename.split('_')[0]]
                gene_sym = genename.split('_')[0]
            else:
                geneid = '-'
                gene_sym = '-'
        else:
            geneid= d_gene[genename]
            gene_sym = genename
        i = i + 1
        l = [str(i),geneid,genename,gene_sym,','.join(cancer_types),','.join(signs)]
        r.write('\t'.join(l)+'\n')

