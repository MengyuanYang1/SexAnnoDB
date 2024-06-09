import sys,os
import pandas as pd
import Proteindata_Ana as PA

def PCIDchange():
    '''
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    d_gene = {}
    for line in f:
        line = line.strip().split('\t')
        if len(line)<5:
            continue
        if line[2] != 'gene':
            continue
        genetype = line[8].split('"')[3]
        if genetype!='protein_coding':
            continue
        genename = line[8].split('"')[7]
        geneid = line[8].split('"')[1].split('.')[0]
        d_gene[genename] = geneid
        #print(genename,geneid)
    '''
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    f.readline()
    d_gene = {}
    for line in f:
        line = line.replace('\n','').split('\t')
        genenames = line[8].replace('"','').split('|')
        #print(genenames)
        genenames = genenames + [line[2]]
        geneid = line[1]
        for genename in genenames:
            d_gene[genename]=geneid
    return d_gene

def CodingGene_mut(cancer_types,genekey):
    sum_path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    for cancer_type in cancer_types:
        #CodingGene_TS(cancer_type,genekey)
        #CodingGene_FMS(cancer_type,genekey)
        pass
    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    #d_gene = PCIDchange()
    for t in keys:
        print(t,genekey)
        CodingGene_sum1(sum_path,t,genekey)
        CodingGene_sum2(t,genekey)
    

def CodingGene_sum1(sum_path,key,genekey):
    #r = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(key,genekey),'w')
    i = 0
    head = ['Index', 'cancer_type', 'geneid', 'genesymbol', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    df_l = []
    for sumf in os.listdir(sum_path):
        if key not in sumf:
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
    #print(newindex[0:10],newindex[-1])
    df_out.index=newindex
    outf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(key,genekey)
    df_out.to_csv(outf,sep='\t')
        

def CodingGene_sum2(t,genekey):
    sumf = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary.txt'.format(t,genekey)
    f = open(sumf)
    head = f.readline()
    d_sum = {}
    for line in f:
        line = line.strip().split()
        if genekey in ['lncRNA','CodingGene']:
            geneid = line[1]
            genesymbol = line[2]
            cancer_type = line[-1] 
            padj = line[8]       
        else:
            cancer_type = line[-1]
            geneid = line[1]
            genesymbol = '-'
            padj = line[7]
        value = cancer_type+'|'+padj
        key = geneid + '|' + genesymbol
        if key in d_sum:
            d_sum[key].append(value)
        else:
            d_sum[key] = [value]
    sumf2 = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{1}/SexAnnoDB_{0}_{1}_Summary2.txt'.format(t,format(genekey))
    r = open(sumf2,'w')
    n = 0 
    head = ['Index','Gene ID','Gene Symbol','Cancer Types','Number','Padjs']
    if genekey  == 'lncRNA':
        head = ['Index','lncRNA ID','lncRNA Symbol','Cancer Types','Number','Padjs']
    if genekey == 'CodingGene':
        head = ['Index','Gene ID','Gene Symbol','Cancer Types','Number','Padjs']
    if genekey == 'miRNA':
        head = ['Index','miRNA','Cancer Types','Number','Padjs']
    r.write('\t'.join(head)+'\n')
    for key in d_sum:
        n = n + 1
        #print(key,d_sum[key],len(d_sum[key]))
        cancer_types = [x.split('|')[0] for x in d_sum[key]]
        padj = [x.split('|')[1] for x in d_sum[key]]
        if genekey in ['lncRNA','CodingGene']:
            line = [str(n)] + key.split('|') + [','.join(cancer_types),str(len(d_sum[key])),','.join(padj)]
        if genekey == 'miRNA':
            line = [str(n),key.split('|')[0]]+ [','.join(cancer_types),str(len(d_sum[key])),','.join(padj)]
        r.write('\t'.join(line)+'\n')

def CodingGene_FMS(cancer_type,genekey):
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    if genekey=='CodingGene':
        f1 = path + 'TCGA-{0}_coding_turmol_female_normal_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_coding_turmol_male_normal_male_deg'.format(cancer_type)
    if genekey == 'lncRNA':
        f1 = path + 'TCGA-{0}_lncRNA_turmol_female_normal_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_lncRNA_turmol_male_normal_male_deg'.format(cancer_type)
    if genekey == 'miRNA':
        f1 = path + 'TCGA-{0}_miRNA_turmol_female_normal_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_miRNA_turmol_male_normal_male_deg'.format(cancer_type)
    lo2gfc_cut = 1
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['padj']<0.05,]
            df1 = df1.loc[df1['log2FoldChange'].abs()>1,]
            df2 = df2.loc[df2['padj']<0.05,]
            df2 = df2.loc[df2['log2FoldChange'].abs()>1,]
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
            df1 = df1.loc[df1['padj']<0.05,]
            df1 = df1.loc[df1['log2FoldChange'].abs()>1,]
            outf1 = outdir+'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
    else:
        if os.path.exists(f2) == True:
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df2 = df2.loc[df2['padj']<0.05,]
            df2 = df2.loc[df2['log2FoldChange'].abs()>1,]
            outf2 = outdir+'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
            df2.to_csv(outf2,sep='\t')

            
def CodingGene_TS(cancer_type,genekey):
    path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/tmp/'.format(genekey)
    outdir = '/data2/myang9/SexAnnoDB/output/Diff_analysis/{0}/Summary/'.format(genekey)
    if genekey == 'CodingGene':
        f1 = path + 'TCGA-{0}_coding_turmol_male_turmol_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_coding_normal_male_normal_female_deg'.format(cancer_type)
    if genekey == 'lncRNA':
        f1 = path + 'TCGA-{0}_lncRNA_turmol_male_turmol_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_lncRNA_normal_male_normal_female_deg'.format(cancer_type)
    if genekey == 'miRNA':
        f1 = path + 'TCGA-{0}_miRNA_turmol_male_turmol_female_deg'.format(cancer_type)
        f2 = path + 'TCGA-{0}_miRNA_normal_male_normal_female_deg'.format(cancer_type)
    lo2gfc_cut = 1
    fdr_cut = 0.05
    if  os.path.exists(f1) == True:
        if os.path.exists(f2) == True:
            ##Tumor specific sex biased gene
            #r = open(outdir + 'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type),'w')
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df2 = pd.read_csv(f2,sep='\t',index_col=0)
            df1 = df1.loc[df1['padj']<0.05,]
            df1 = df1.loc[df1['log2FoldChange'].abs()>1,]
            df2 = df2.loc[df2['padj']<0.05,]
            df2 = df2.loc[df2['log2FoldChange'].abs()>1,]
            tumor_ids = [x for x in df1._stat_axis.values.tolist() if x not in df2._stat_axis.values.tolist()]
            normal_ids = [x for x in df2._stat_axis.values.tolist() if x not in df1._stat_axis.values.tolist()]
            df1 = df1.T
            df1 = df1[tumor_ids]
            df1 = df1.T
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')
        else:
            #r = open(outdir + 'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type),'w')
            df1 = pd.read_csv(f1,sep='\t',index_col=0)
            df1 = df1.loc[df1['padj']<0.05,]
            df1 = df1.loc[df1['log2FoldChange'].abs()>1,]
            outf1 = outdir+'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
            df1.to_csv(outf1,sep='\t')

    else:
        print(f1 +  '  No result!')


def PCheeatmap(cancer_types):
    #PCheeatmap_Summaryfile(cancer_types)
    #PCheeatmap_makefile(cancer_types)
    PCheeatmap_Drawfile(cancer_types)

def PCheeatmap_Drawfile(cancer_types):
    path = '/data2/myang9/SexAnnoDB/output/Figures/PCheatmap/file/'
    r = open('/data2/myang9/SexAnnoDB/output/Figures/PCheatmap/PCheatmap_run.sh','w')
    i = 0
    for inf in os.listdir(path):
        inf = path + inf
        svgf = inf.replace('file','figure')
        svgf = svgf.replace('.txt','_heatmap.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/heatmap_gene.r {0} {1}'.format(inf,svgf)
        r.write('nohup ' + cmd + ' &\n')
        i = i + 1
        if i%15 == 0:
            r.write('wait;\n')

def PCheeatmap_makefile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    Geneids = [line.strip().split()[1] for line in f]
    f.close()
    f = open('/data2/myang9/SexAnnoDB/output/Figures/PCheatmap/PCheatmap_summary.txt')
    head = f.readline().strip().split(',')
    outdir = '/data2/myang9/SexAnnoDB/output/Figures/PCheatmap/'
    for line in f:
        line = line.strip().split(',')
        geneid = line[0]
        if geneid not in Geneids:
            continue
        #r = open(outdir +'file/{0}.txt'.format(geneid),'w')
        d_out = []
        newhead =['Cancer_type']+ [x[5:] for x in head[1:5]]
        #r.write(','.join(newhead)+'\n')
        d_out.append(newhead)
        for i in range(1,len(head),4):
            cancer_type = head[i][0:4]
            newl = [cancer_type] + line[i:i+4]
            #r.write(','.join(newl)+'\n')
            d_out.append(newl)
        df = pd.DataFrame(d_out)
        df = df.T
        outf = outdir +'file/{0}.txt'.format(geneid)
        df.to_csv(outf,header=0,index=False)


def PCheeatmap_Summaryfile(cancer_types):
    d_Patientinfo = PA.PatientSubgroup()
    df_out = []
    path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    colname = []
    for cancer_type in cancer_types:
        #if cancer_type != 'BRCA':
        #    continue
        exp_f = path + 'TCGA-{0}_expression_FPKM.csv'.format(cancer_type)
        df = pd.read_csv(exp_f,sep=',',index_col=0)
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
            df_mean = df_tmp.mean(skipna=True)
            df_out.append(df_mean)
    result = pd.concat(df_out,axis=1)
    result.columns = colname
    #print(result)
    outf = '/data2/myang9/SexAnnoDB/output/Figures/PCheatmap/PCheatmap_summary.txt'
    result.to_csv(outf)

def PatientInfo():
    f = open('/data2/myang9/SexAnnoDB/data/TCGA/cancer_type_info_27.txt')
    d_Patientinfo2 = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        cancer_type = line[3]
        sex_type = line[1][0].upper() 
        tumortype = line[2][0].upper() 
        key = cancer_type+':'+sex_type+'_'+tumortype
        p_id = line[0][0:15]
        d_Patientinfo2[p_id] = key
    return d_Patientinfo2


def Coding_boxplot(cancer_types):
    #Coding_boxplot_step1(cancer_types)
    #Coding_boxplot_step2(cancer_types)
    #Coding_boxplot_step3(cancer_types)
    Coding_boxplot_step4()

def Coding_boxplot_step1(cancer_types):
    d_Patientinfo2 = PatientInfo()
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    Geneids = [line.strip().split()[1] for line in f]
    f.close()
    path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    d_out = {}
    l_df = []
    for cancer_type in cancer_types:
        #if cancer_type != 'BRCA':
        #    continue
        exp_f = path + 'TCGA-{0}_expression_FPKM.csv'.format(cancer_type)
        df = pd.read_csv(exp_f,sep=',',index_col=0)
        patient_ids =[d_Patientinfo2[x[0:15]] for x in df.columns.values.tolist()]
        df.columns = patient_ids
        df = df.T
        Geneids = [x for x in Geneids if x in df.columns.values.tolist()]
        df = df[Geneids]
        print(df)
        #df = df.T
        l_df.append(df)
    df_sum = pd.concat(l_df)
    df_sum = df_sum.T
    outf = '/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/SexDiff_forDIffGeneBoxplot.txt'
    df_sum.to_csv(outf)
        
def Coding_boxplot_step2(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/SexDiff_forDIffGeneBoxplot.txt')
    head = f.readline().strip().split(',')
    outdir = '/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/file/'
    for line in f:
        line = line.strip().split(',')
        geneid = line[0]
        outf = outdir + '{0}.txt'.format(geneid)
        r = open(outf,'w')
        newh =['FPKM','cancer_type','group']
        r.write('\t'.join(newh)+'\n')
        for i in range(1,len(head)):
            newl = [line[i]] + head[i].split(':')
            r.write('\t'.join(newl)+'\n')
def Coding_boxplot_step3(cancer_types):
    types = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    d_match = {}
    for DEGtype in types:
        if DEGtype == 'Tumor_specific_sex_biased':
            groups = ['F_T','M_T']
        if  DEGtype == 'Male_specific':
            groups = ['M_N','M_T']
        if DEGtype == 'Female_specific':
            groups = ['F_T','F_N']
        f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/CodingGene/SexAnnoDB_{0}_PC_Summary2.txt'.format(DEGtype))
        f.readline()
        for line in f:
            line = line.strip().split()
            geneid = line[1]
            cancer_types = line[3].split(',')
            pajs = line[5].split(',')
            d_padj = {}
            pajs_sign = []
            for padj in pajs:
                if float(padj)<0.0001:
                    pajs_sign.append('****')
                elif float(padj)<0.001:
                    pajs_sign.append('***')
                elif float(padj)<0.01:
                    pajs_sign.append('**')
                elif float(padj)<0.05:
                    pajs_sign.append('*')
                else:
                    pajs_sign.append('ns')
            for i in range(0,len(cancer_types)):
                d_padj[cancer_types[i]] = pajs_sign[i]
            if os.path.exists('/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/file/{0}.txt'.format(geneid)) == False:
                continue
            f_tmp = open('/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/file/{0}.txt'.format(geneid))
            r_tmp = open('/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/figure/{0}/file/{1}.txt'.format(DEGtype,geneid),'w')
            head = f_tmp.readline().strip().split()
            head = head + ['sign']
            r_tmp.write('\t'.join(head)+'\n')
            for l in f_tmp:
                l = l.strip().split()
                if l[1] not in cancer_types:
                    continue
                if l[2] not in groups:
                    continue
                l = l + [d_padj[l[1]]]
                #print(l)
                r_tmp.write('\t'.join(l)+'\n')

def Coding_boxplot_step4():
    path = '/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/figure/'
    types = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    Rkeys = ['TS','MS','FS']
    for i in range(0,len(types)):
        keytype = types[i]
        Rkey = Rkeys[i]
        r = open('/data2/myang9/SexAnnoDB/output/Figures/DiffGeneBoxplot/{0}_draw.sh'.format(keytype),'w')
        inpath = path + '{0}/file/'.format(keytype)
        outpath = path + '{0}/figure/'.format(keytype)
        #print(inpath,outpath)
        Rscript = '/data2/myang9/SexAnnoDB/script/Draw/boxplot_diffGene_{0}.r'.format(Rkey)
        i = 0
        for inf in os.listdir(inpath):
            outf = inf.replace('.txt','.svg')
            inf = inpath+ inf
            outf = outpath+outf
            cmd = 'Rscript {0} {1} {2}'.format(Rscript,inf,outf)
            i = i + 1
            r.write('nohup '+cmd + ' &\n')
            if i%10==0:
                r.write('wait;\n')
            


def Summary_CanType(cancer_types):
    keys = ['Tumor_specific_sex_biased','Male_specific','Female_specific']
    sum_path = '/data2/myang9/SexAnnoDB/output/Diff_analysis/CodingGene/Summary/'
    for cancer_type in cancer_types:
        print(cancer_type)
        d_out = {}
        sumf1 = sum_path + 'SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
        if os.path.exists(sumf1)==False:
            continue
        print(sumf1)
        f = open(sumf1)
        for line in f:
            line = line.strip().split()
            geneid = line[0]+'\t'+line[1]
            if geneid in d_out:
                d_out[geneid][0] = 'Yes'
            else:
                d_out[geneid] = ['No']*3
                d_out[geneid][0] = 'Yes'
        f.close()
        sumf1 = sum_path + 'SexAnnoDB_Male_specific_PC_{0}.txt'.format(cancer_type)
        if os.path.exists(sumf1)==True:
            f = open(sumf1)
            for line in f:
                line = line.strip().split()
                geneid = line[0]+'\t'+line[1]
                if geneid in d_out:
                    d_out[geneid][1] = 'Yes'
                else:
                    d_out[geneid] = ['No']*3
                    d_out[geneid][1] = 'Yes'
            f.close()
        sumf1 = sum_path + 'SexAnnoDB_Female_specific_PC_{0}.txt'.format(cancer_type)
        if os.path.exists(sumf1)==True:
            f = open(sumf1)
            for line in f:
                line = line.strip().split()
                geneid = line[0]+'\t'+line[1]
                if geneid in d_out:
                    d_out[geneid][2] = 'Yes'
                else:
                    d_out[geneid] = ['No']*3
                    d_out[geneid][2] = 'Yes'
            f.close()
        outf = sum_path + 'SexAnnoDB_Index_Page_{0}.txt'.format(cancer_type)
        head = ['Index','Gene_ID','Gene_Symbol','Sex_biased_Coding_Gene','Male_specific_DEGs','Female_specific_DEGs']
        r = open(outf,'w')
        r.write('\t'.join(head)+'\n')
        for key in d_out:
            if key == 'coding_id\tcoding_name':
                continue
            neewl = [key] + d_out[key]
            r.write('\t'.join(neewl)+'\n')





