import os,sys
import Proteindata_Ana as PA
import pandas as pd


def Editing_preprocessing():
    path = '/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/'
    infs = [x for x in os.listdir(path) if 'RNAediting' in x]
    for inf in infs:
        inf = path + inf
        df = pd.read_csv(inf,sep='\t',index_col=0)
        df['NaNCount'] = df.isnull().sum(axis=1)
        n = df.shape[1]*0.3
        df_filtered = df[df['NaNCount'] < n]
        outf = inf+'.filter'
        df_filtered.to_csv(outf,sep='\t')

def Editing_diff():
    d_Patientinfo = PA.PatientSubgroup()
    cancer_types = set([x.split('_')[0] for x in d_Patientinfo.keys()])
    path = '/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/'
    for  cancer_type in cancer_types:
        editing_f = path + cancer_type + '.RNAediting'
        df = pd.read_csv(editing_f,sep='\t',index_col=0)
        patient_ids =[x[0:15] for x in df.columns.values.tolist()]
        df.columns = patient_ids
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
        outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/AtoIediting/tmp2024/'
        df1 = df_female_tumor
        df2 = df_female_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            df_out['dFre'] = df_out['ave1'] - df_out['ave2']
            outf = outdir + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type) 
            df_out.to_csv(outf,sep='\t')
        df1 = df_male_tumor
        df2 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            df_out['dFre'] = df_out['ave1'] - df_out['ave2']
            outf = outdir + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df2 = df_female_tumor
        df1 = df_male_tumor
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            df_out['dFre'] = df_out['ave1'] - df_out['ave2']
            outf = outdir + '{0}_Tumor_maleVSfemale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')
        df2 = df_female_normal
        df1 = df_male_normal
        if df1.shape[1]>=5 and df2.shape[1]>=5:
            df_out = PA.differential_Wilcoxontest(df1,df2)
            df_out['dFre'] = df_out['ave1'] - df_out['ave2']
            outf = outdir + '{0}_Normal_maleVSfemale_diffResult.txt'.format(cancer_type)
            df_out.to_csv(outf,sep='\t')

def Editing_Anno_1():
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/Anno/FunctionCancer.txt')
    d_func = {}
    for line in f:
        line = line.strip().split('\t')
        id2 = line[2]+':'+line[6]
        info = [line[3]]+line[7:9]
        if id2 in d_func:
            d_func[id2].append(','.join(info))
        else:
            d_func[id2] = [','.join(info)]
    f.close()
    '''
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/Anno/SequenceCancer.txt')
    d_protein = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        id2 = line[4]
        info = line[1:]
        if id2 not in d_func:
            continue
        print(line)
        print(id2,set(d_func[id2]))
        if id2 in d_protein:
            d_protein[id2].append(','.join(info))
        else:
            d_protein[id2]=[','.join(info)]
    '''
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/Anno/editing_info.txt')
    d_info = {}
    for line in f:
        line = line.strip().split()
        if line[5] != 'protein_coding':
            continue
        id2 = '_'.join(line[1:4])+':'+line[7]
        if id2 not in d_func:
            continue
        else:
            func = list(set(d_func[id2]))[0]
        info = [line[0],line[4].split('.')[0],line[6],line[7].split('.')[0],line[9]]
        newl = [id2.split(':')[0]]+info+func.split(',')
        d_info[id2] = newl
    f.close()
    r = open('/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/Anno/SexAnnoDB_RNAediting_Annotation.txt','w')
    head = ['Position','CAID','GeneID','GeneSymbol','TransID','TransSymbol','VariantType','NTchange','AAchange']
    r.write('\t'.join(head)+'\n')
    for key in d_info:
        r.write('\t'.join(d_info[key])+'\n')
    return d_info

def Editing_Anno_2():
    d_info = {}
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/Anno/SexAnnoDB_RNAediting_Annotation.txt')
    for line in f:
        line = line.strip().split('\t')
        key = line[0]
        CAID= line[1]
        
        if key in d_info:
            d_info[key][2].append(line[4])
            d_info[key][3].append(line[5])
            d_info[key][4].append(line[6])
            d_info[key][5].append(line[7])
            d_info[key][6].append(line[8])
        else:
            d_info[key] = [line[2],line[1],[line[4]],[line[5]],[line[6]],[line[7]],[line[8]]]
    return d_info

def Editing_Summary(cancer_types,genekey):
    sum_path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary2024/'.format(genekey)
    for cancer_type in cancer_types:
        pass
        Editing_TS(cancer_type,genekey)
        Editing_FMS(cancer_type,genekey)
    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    for t in keys:
        print(t,genekey)
        EditingSummary_1(sum_path,t,genekey)
        EditingSummary_2(sum_path,t,genekey)
    EditingSummary_3(sum_path,keys,genekey)

def EditingSummary_2(sum_path,t,genekey):
    d_info = Editing_Anno_2()
    sumf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    f = open(sumf)
    head = f.readline()
    d_sum = {}
    for line in f:
        line = line.strip().split()
        editingid = line[1]
        if editingid not in d_info:
            continue
        cancer_type = line[-1]
        if editingid in d_sum:
            d_sum[editingid].append(cancer_type)
        else:
            d_sum[editingid] = [cancer_type]
    sumf2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,genekey)
    r = open(sumf2,'w')
    i = 0
    head = ['Index','Gene_ID','Editing_Position','CAeditome_ID','Transcripts','Variant_Type','NTchange','AAchange','Cancer_Type']
    r.write('\t'.join(head)+'\n')
    for key in d_sum:
        i = i + 1
        cancer_types = d_sum[key]
        info = [','.join(x) for x in d_info[key][3:]]
        newl = [str(i),d_info[key][0],key,d_info[key][1]] +  info + [','.join(cancer_types)]
        r.write('\t'.join(newl)+'\n')

def EditingSummary_3(sum_path,Analysis,genekey):
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/SexAnnoDB_{0}_Summary3.txt'.format(genekey)
    d_out = {}
    for i in range(0,len(Analysis)):
        t = Analysis[i]
        sumf2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,genekey)
        f = open(sumf2)
        head = f.readline()
        #print(head.split())
        for line in f:
            line = line.strip().split('\t')
            key = '\t'.join(line[1:8])
            if key in d_out:
                d_out[key][i] = line[-1]
            else:
                d_out[key] = ['-']*len(Analysis)
                d_out[key][i] = line[-1]
        f.close()
    i = 0 
    head = ['Index', 'Gene_ID', 'Editing_Position', 'CAeditome_ID', 'Transcripts', 'Variant_Type', 'NTchange', 'AAchange','Sex_biased','Male_specific','Female_specific' ]
    r = open(outf,'w')
    r.write('\t'.join(head)+'\n')
    for key in d_out:
        i = i + 1
        newl = [str(i),key]+d_out[key]
        r.write('\t'.join(newl)+'\n')
    


def EditingSummary_1(sum_path,t,genekey):
    i = 0
    df_l = []
    for sumf in os.listdir(sum_path):
        if t not in sumf:
            continue
        #print(sumf)
        cancer_type = sumf.split('_')[-1].split('.')[0]
        sumf = sum_path + sumf
        df = pd.read_csv(sumf,sep='\t') #,index_col=0)
        df['cancer_type'] = cancer_type
        print(df)
        df_l.append(df)
    df_out = pd.concat(df_l)
    df_out.index.name='Index'
    newindex = [str(x) for x in range(1,df_out.shape[0]+1)]
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    df_out.to_csv(outf,sep='\t')

def Editing_FMS(cancer_type,genekey):
    path ='/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp2024/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary2024/'.format(genekey)
    if genekey == 'AtoIediting':
        f1 = path + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type)
        f2 = path + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            female_ids = [x for x in df1._stat_axis.values.tolist() if x not in df2._stat_axis.values.tolist()]
            male_ids = [x for x in df2._stat_axis.values.tolist() if x not in df1._stat_axis.values.tolist()]
            df1 = df1.T
            df1 = df1[female_ids]
            df1 = df1.T
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            df2 = df2.T
            df2 = df2[male_ids]
            df2 = df2.T
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
            df2.to_csv(outf2,sep='\t')
        else:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            #df1 = df1.loc[(df1['p_adj']<0.05) & (abs(df1['dFre'])>0.1),]
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
    else:
        if os.path.exists(f2) == True:
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df2 = df2.loc[df2['p_adj']<0.05,]
            #df2 = df2.loc[(df2['p_adj']<0.05) & (abs(df2['dFre'])>0.1),]
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df2.to_csv(outf2,sep='\t')


def Editing_TS(cancer_type,genekey):
    path ='/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp2024/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary2024/'.format(genekey)
    if genekey == 'AtoIediting':
        f1 = path + '{0}_Tumor_maleVSfemale_diffResult.txt'.format(cancer_type)
        f2 = path + '{0}_Normal_maleVSfemale_diffResult.txt'.format(cancer_type)
    fdr_cut = 0.05
    print(cancer_type)
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            ##Tumor specific sex biased gene
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            #df1 = df1.loc[(df1['p_adj']<0.05) & (abs(df1['dFre'])>0.1),]
            #df2 = df2.loc[(df2['p_adj']<0.05) & (abs(df2['dFre'])>0.1),]
            tumor_ids = [x for x in df1._stat_axis.values.tolist() if x not in df2._stat_axis.values.tolist()]
            normal_ids = [x for x in df2._stat_axis.values.tolist() if x not in df1._stat_axis.values.tolist()]
            df1 = df1.T
            df1 = df1[tumor_ids]
            df1 = df1.T
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
        else:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')

def Edheeatmap(cancer_types):
    Edheeatmap_Summaryfile(cancer_types)
    #Edheeatmap_makefile(cancer_types)
    #Edheeatmap_Drawfile(cancer_types)
    
def Edheeatmap_Summaryfile(cancer_types):
    d_Patientinfo = PA.PatientSubgroup()
    f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/AtoIediting/SexAnnoDB_AtoIediting_Summary3.txt')
    edtings = [line.strip().split()[2] for line in f]
    edtings = ['chr12_68851361_-']
    f.close()
    cancer_types = ['SARC']  #set([x.split('_')[0] for x in d_Patientinfo.keys()])
    path = '/data2/myang9/SexAnnoDB/data/TCGA/RNAediting/'
    colname = []
    df_out = []
    l_df = []
    for  cancer_type in cancer_types:
        editing_f = path + cancer_type + '.RNAediting'
        df = pd.read_csv(editing_f,sep='\t',index_col=0)
        patient_ids =[x[0:15] for x in df.columns.values.tolist()]
        df.columns = patient_ids
        df = df.T
        edtings = [x for x in edtings if x in df.columns.values.tolist()]
        df = df[edtings]
        df = df.T
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
        #for df_tmp in df_female_tumor:
        df_tmp = df_female_tumor.T.dropna()
        df_tmp['group'] = 'female_tumor'
        df_tmp['cancer'] = cancer_type
        print(df_tmp)
        l_df.append(df_tmp)
        df_tmp = df_male_tumor.T.dropna()
        df_tmp['group'] = 'male_tumor'
        df_tmp['cancer'] = cancer_type
        print(df_tmp)
        l_df.append(df_tmp)
        #df_mean = df_tmp.mean(skipna=True)
        #df_out.append(df_mean)
    #result = pd.concat(df_out,axis=1)
    #result.columns = colname
    #outf = '/data2/myang9/SexAnnoDB/output/Figures/Edheatmap/Edheatmap_summary.txt'
    #result.to_csv(outf)
    df_out = pd.concat(l_df)
    print(df_out)
    outf = '/data2/myang9/SexAnnoDB/output/Revise2024/Editing_candidate_SARC_CPM.txt'
    df_out.to_csv(outf,sep='\t')

def Edheeatmap_makefile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/AtoIediting/SexAnnoDB_AtoIediting_Summary3.txt')
    d_eding = {}
    for line in f:
        line = line.strip().split()
        edting = line[2]
        geneid = line[1]
        if geneid in d_eding:
            d_eding[geneid].append(edting)
        else:
            d_eding[geneid] = [edting]
    f.close()
    #for geneid in d_eding:
    #    print(geneid,len(d_eding[geneid]))
    inf = '/data2/myang9/SexAnnoDB/output/Figures/Edheatmap/Edheatmap_summary.txt'
    df = pd.read_csv(inf,sep=',',index_col=0)
    df = df.T
    for geneid in d_eding:
        outf = '/data2/myang9/SexAnnoDB/output/Figures/Edheatmap/file/{0}_heatmap.txt'.format(geneid)
        Edids = [x for x in d_eding[geneid] if x in df.columns.values.tolist()]
        df_tmp = df[Edids]
        df_tmp = df_tmp.T
        #df_tmp = df_tmp.dropna(axis=1,how='all')
        #df_tmp = df_tmp.dropna(axis=0,how='all')
        if df_tmp.shape[0] == 0:
            continue
        df_tmp.to_csv(outf,sep=',')


def Edheeatmap_Drawfile(cancer_types):
    path = '/data2/myang9/SexAnnoDB/output/Figures/Edheatmap/file/'
    r = open('/data2/myang9/SexAnnoDB/output/Figures/Edheatmap/Edheatmap_run.sh','w')
    i = 0
    for inf in os.listdir(path):
        geneid = inf.split('_')[0]
        inf = path + inf
        svgf = inf.replace('file','figure')
        svgf = svgf.replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/heatmap_Ed.r {0} {1}'.format(inf,svgf)
        r.write('nohup ' + cmd + ' &\n')
        i = i + 1
        if i%15 == 0:
            r.write('wait;\n')














