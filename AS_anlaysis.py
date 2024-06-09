import os,sys
import pandas as pd
import Proteindata_Ana as PA


def AS_diffalaysis():
    d_Patientinfo = PA.PatientSubgroup()
    cancer_types = set([x.split('_')[0] for x in d_Patientinfo.keys()])
    AS_path = '/data2/myang9/DriverRBP/data/ASevents/'
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        as_f = AS_path + cancer_type + '_exon_skip_ASmatrix.txt'
        df = pd.read_csv(as_f,sep='\t',index_col=0)
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
        #print(patient_ids[0:10])
        print(cancer_type,len(l_female_tumor),len(l_female_normal),len(l_male_tumor),len(l_male_normal))
        
        df_female_tumor = df[l_female_tumor]
        df_female_normal = df[l_female_normal]
        df_male_tumor = df[l_male_tumor]
        df_male_normal = df[l_male_normal]
        #print(df_female_tumor,df_female_normal,df_male_tumor,df_male_normal)
        outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/ASevents/'
        df1 = df_female_tumor
        df2 = df_female_normal
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
def ORF_anno():
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/ASevents/exon_skip_annotated_with_transcript_4sitematch_with_gsymbol_all_ORF.txt')
    d_anno = {}
    for line in f:
        #print(line)
        line = line.strip().split()
        region = ':'.join(line[2:4])
        d_anno[region] = line[-1]
    f.close()
    '''
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/ASevents/exon_skip_annotated_with_transcript_4sitematch_with_gsymbol_In-Frame_only_AA_loc_lost_protein_features.txt')
    d_protein = {}
    for line in f:
        line = line.strip().split()
        print(line)
    '''
    f = open('/data2/myang9/DriverRBP/data/information/exon_skip_AS_information.txt')
    d_orf = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        esid = line[0]
        region = line[-1]
        if region not in d_anno:
            continue
        anno = d_anno[region]
        #print(esid,region,anno)
        d_orf[esid]=anno
    return d_orf

def EXid_gene():
    f = open('/data2/myang9/DriverRBP/data/information/exon_skip_AS_information_4sitematch.txt')
    d_exidchange = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        #print(line)
        key = line[1]
        info = [line[6],line[0],key]+line[4:6] #+[line[8]]
        d_exidchange[key]=info
    return d_exidchange

def ASSummary(cancer_types,genekey):
    sum_path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    for cancer_type in cancer_types:
        pass
    '''
        ExonSkipping_TS(cancer_type,genekey)
        ExonSkipping_FMS(cancer_type,genekey)
    '''
    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    for t in keys:
        print(t,genekey)
        #ASSummary_1(sum_path,t,genekey)
        #ASSummary_2(sum_path,t,genekey)
        ASSummary_2_new(sum_path,t,genekey)
    #ASSummary_3(sum_path,keys,genekey)

def ASSummary_3(sum_path,keys,genekey):
    print(sum_path,keys,genekey)
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/ASevents/'
    d_out = {}
    for i in range(0,len(keys)):
        inf = path + 'SexAnnoDB_{0}_{1}_Summary2.txt'.format(keys[i],genekey)
        f = open(inf)
        head = f.readline().split()
        for line in f:
            line = line.strip().split()
            key = '|'.join(line[1:7])
            value = line[-2]
            if key in d_out:
                d_out[key][i]= value
            else:
                d_out[key] = ['-']*3
                d_out[key][i] = value
    i = 0
    head = ['Index', 'Gene_ID', 'Gene_Symbol', 'ES_ID', 'exon_skipping_position', 'Skipped_region', 'ORF_anno',  'Sex_biased_ES_events','Male_specific_DEESs','Female_specific_DEESs']
    r = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/ASevents/SexAnnoDB_diffES_Summary3.txt','w')
    r.write('\t'.join(head)+'\n')
    for key in d_out:
        i = i+1
        wl = [str(i)] + key.split('|')+d_out[key]
        r.write('\t'.join(wl)+'\n')
    

def ASSummary_1(sum_path,t,genekey):
    i = 0
    df_l = []
    for sumf in os.listdir(sum_path):
        if t not in sumf:
            continue
        print(sumf)
        cancer_type = sumf.split('_')[-1].split('.')[0]
        #print(cancer_type)
        sumf = sum_path + sumf
        df = pd.read_csv(sumf,sep='\t') #,index_col=0)
        df['cancer_type'] = cancer_type
        df_l.append(df)
    df_out = pd.concat(df_l)
    df_out.index.name='Index'
    newindex = [str(x) for x in range(1,df_out.shape[0]+1)]
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    df_out.to_csv(outf,sep='\t')

def ASSummary_2(sum_path,t,genekey):
    d_exidchange = EXid_gene()
    d_orf = ORF_anno()
    sumf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    f = open(sumf)
    head = f.readline()
    d_sum = {}
    for line in f:
        line = line.strip().split()
        exid = line[1]
        cancer_type = line[-1]
        if exid in d_sum:
            d_sum[exid].append(cancer_type)
        else:
            d_sum[exid] = [cancer_type]
    sumf2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,format(genekey))
    r = open(sumf2,'w')
    i = 0
    head = ['Index','Gene_ID','Gene_Symbol','ES_ID','exon_skipping_position','Skipped_region','ORF_anno','cancer_types','Number']
    r.write('\t'.join(head)+'\n')
    for key in d_sum:
        i = i + 1
        cancer_types = d_sum[key]
        if key not in d_exidchange:
            continue
        info = d_exidchange[key]
        if key not in d_orf:
            anno = '-'
        else:
            anno = d_orf[key]
        newl = [str(i)] + info +[anno] + [','.join(cancer_types),str(len(cancer_types))]
        r.write('\t'.join(newl)+'\n')

def ASSummary_2_new(sum_path,t,genekey):
    d_exidchange = EXid_gene()
    d_orf = ORF_anno()
    sumf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    f = open(sumf)
    sumf2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2_new.txt'.format(t,genekey)
    r = open(sumf2,'w')
    head = f.readline().strip().split()
    head = head + ['anno']
    r.write('\t'.join(head)+'\n')
    for line in f:
        line = line.strip().split()
        exid = line[1]
        cancer_type = line[-1]
        if exid not in d_exidchange:
            continue
        gene = d_exidchange[exid][1]
        #print(gene)
        if exid not in d_orf:
            anno = '-'
        else:
            anno = d_orf[exid]
        newl = line +[anno,gene]
        r.write('\t'.join(newl)+'\n')

def ExonSkipping_TS(cancer_type,genekey):
    path ='/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    if genekey == 'ASevents':
        f1 = path + '{0}_Tumor_femaleVSmale_diffResult.txt'.format(cancer_type)
        f2 = path + '{0}__Normal_femaleVSmale_diffResult.txt'.format(cancer_type)
    dPSI=0.1
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            ##Tumor specific sex biased gene
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            df1['dPSI'] = df1['ave1']-df1['ave2']
            df1 = df1.loc[df1['dPSI'].abs()>0.1,]
            df2['dPSI'] = df2['ave1']-df2['ave2']
            df2 = df2.loc[df2['dPSI'].abs()>0.1,]
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
            df1['dPSI'] = df1['ave1']-df1['ave2']
            df1 = df1.loc[df1['dPSI'].abs()>0.1,]
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')

def ExonSkipping_FMS(cancer_type,genekey):
    path ='/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    if genekey == 'ASevents':
        f1 = path + '{0}_female_tumVSnor_diffResult.txt'.format(cancer_type)
        f2 = path + '{0}_male_tumVSnor_diffResult.txt'.format(cancer_type)
    dPSI=0.1
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['p_adj']<0.05,]
            df2 = df2.loc[df2['p_adj']<0.05,]
            df1['dPSI'] = df1['ave1']-df1['ave2']
            df1 = df1.loc[df1['dPSI'].abs()>0.1,]
            df2['dPSI'] = df2['ave1']-df2['ave2']
            df2 = df2.loc[df2['dPSI'].abs()>0.1,]
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
            df1['dPSI'] = df1['ave1']-df1['ave2']
            df1 = df1.loc[df1['dPSI'].abs()>0.1,]
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
    else:
        if os.path.exists(f2) == True:
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df2 = df2.loc[df2['p_adj']<0.05,]
            df2['dPSI'] = df2['ave1']-df2['ave2']
            df2 = df2.loc[df2['dPSI'].abs()>0.1,]
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df2.to_csv(outf2,sep='\t')

def ASheeatmap(cancer_types):
    ASheeatmap_Summaryfile(cancer_types)
    #ASheeatmap_makefile(cancer_types)
    #ASheeatmap_Drawfile(cancer_types)

def ASheeatmap_Summaryfile(cancer_types):
    d_Patientinfo = PA.PatientSubgroup()
    df_out = []
    path = '/data2/myang9/DriverRBP/data/ASevents/'
    colname = []
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        exp_f = path + '{0}_exon_skip_ASmatrix.txt'.format(cancer_type)
        df = pd.read_csv(exp_f,sep='\t',index_col=0)
        patient_ids =[x[0:15] for x in df.columns.values.tolist()]
        df.columns = patient_ids
        #print(df)
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
            df_mean = df_tmp.mean(skipna=True)
            print(df_mean)
            df_out.append(df_mean)
    result = pd.concat(df_out,axis=1)
    result.columns = colname
    outf = '/data2/myang9/SexAnnoDB/output/Figures/ASheatmap/ASheatmap_summary.txt'
    result.to_csv(outf)
    

def ASheeatmap_makefile(cancer_types):
    f = open('/data2/myang9/DriverRBP/data/information/exon_skip_AS_information_4sitematch.txt')
    f.readline()
    d_exid = {}
    for line in f:
        line = line.strip().split()
        geneid = line[6]
        exid = line[1]
        if geneid in d_exid:
            d_exid[geneid].append(exid)
        else:
            d_exid[geneid] = [exid]
    f.close()
    inf = '/data2/myang9/SexAnnoDB/output/Figures/ASheatmap/ASheatmap_summary.txt'
    df = pd.read_csv(inf,sep=',',index_col=0)
    df = df.T
    for geneid in d_exid:
        outf = '/data2/myang9/SexAnnoDB/output/Figures/ASheatmap/file/{0}_heatmap.txt'.format(geneid)
        exids = [x for x in d_exid[geneid] if x in df.columns.values.tolist()]
        df_tmp = df[exids]
        df_tmp = df_tmp.T
        #df_tmp = df_tmp.dropna(axis=1,how='all')
        #df_tmp = df_tmp.dropna(axis=0,how='all')
        print(df_tmp.shape)
        if df_tmp.shape[0] == 0:
            continue
        df_tmp.to_csv(outf,sep=',')

def ASheeatmap_Drawfile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    geneids = [line.strip().split()[1] for line in f]
    f.close()
    path = '/data2/myang9/SexAnnoDB/output/Figures/ASheatmap/file/'
    r = open('/data2/myang9/SexAnnoDB/output/Figures/ASheatmap/ASheatmap_run.sh','w')
    i = 0

    for inf in os.listdir(path):
        geneid = inf.split('_')[0]
        if geneid not in geneids:
            continue
        inf = path + inf
        svgf = inf.replace('file','figure')
        svgf = svgf.replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/heatmap_AS.r {0} {1}'.format(inf,svgf)
        r.write('nohup ' + cmd + ' &\n')
        i = i + 1
        if i%15 == 0:
            r.write('wait;\n')
