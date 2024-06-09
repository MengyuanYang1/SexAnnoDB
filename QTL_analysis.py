import os,sys 
import numpy as np
import pandas as pd
import Proteindata_Ana as PA
import warnings
warnings.filterwarnings("ignore")
import multiprocessing

def GetQTLinput_multip(cancer_types,Ana_type):
    task_l = []
    task_l2 = []
    d_Patientinfo = PA.PatientSubgroup()
    for cancer_type in cancer_types:
        d_group= GetQTLinput_group(cancer_type,Ana_type,d_Patientinfo)
        print(cancer_type)
        mol='FPKM'
        print(cancer_type,Ana_type,mol)
        GetQTLinput(cancer_type,d_group,Ana_type,mol)
        #GeteQTLinput(cancer_type,d_Patientinfo)
        '''
        for sex in ['female','male']:
            l = (cancer_type,sex)
            task_l2.append(l)
            #GeteQTLinput2(cancer_type,sex)
    print('Run GeteQTLinput_multip for {0} cancer_type!'.format(len(task_l)))
    cores = multiprocessing.cpu_count()
    cores = 10
    pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(GeteQTLinput,task_l)
    #pool.starmap(GeteQTLinput2,task_l2)
'''

def GetQTLinput_group(cancer_type,Ana_type,d_Patientinfo):
        if Ana_type == 'eQTL':
            f1 = '/data2/myang9/SexAnnoDB/data/GDCdata/TCGA-{0}_expression_FPKM.csv'.format(cancer_type)
        if Ana_type == 'sQTL':
            f1 = '/data2/myang9/DriverRBP/data/ASevents/{0}_exon_skip_ASmatrix.txt'.format(cancer_type)
        df1 = pd.read_csv(f1,sep=',',index_col=0,nrows=5)
        f2 = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.gen'.format(cancer_type)
        df2 = pd.read_csv(f2,sep='\t',index_col=0,nrows=5)
        df1_col = [x[0:15] for x in df1.columns.values.tolist()]
        df2_col = [x[0:15] for x in df2.columns.values.tolist()]
        patient_ids = [x for x in df1_col if x in df2_col]
        d_group = {}
        if cancer_type+'_turmol_female' not in d_Patientinfo:
            l_female_tumor = []
        else:
            l_female_tumor = [x for x in d_Patientinfo[cancer_type+'_turmol_female'] if x in patient_ids]
        d_group['female_tumor'] = l_female_tumor
        if cancer_type+'_normal_female' not in d_Patientinfo:
            l_female_normal = []
        else:
            l_female_normal = [x for x in d_Patientinfo[cancer_type+'_normal_female'] if x in patient_ids ]
        d_group['female_normal'] = l_female_normal 
        if cancer_type+'_turmol_male' not in d_Patientinfo:
            l_male_tumor = []
        else:
            l_male_tumor =[x for x in  d_Patientinfo[cancer_type+'_turmol_male'] if x in patient_ids]
        d_group['male_tumor'] = l_male_tumor
        if cancer_type+'_normal_male' not in d_Patientinfo:
            l_male_normal = []
        else:
            l_male_normal =[x for x in  d_Patientinfo[cancer_type+'_normal_male'] if x in patient_ids]
        d_group['male_normal'] = l_male_normal
        print(cancer_type,'eQTL',len(l_female_tumor),len(l_female_normal),len(l_male_tumor),len(l_male_normal))
        return d_group       
        
def GetQTLinput(cancer_type,d_group,Ana_type,mol):
    if mol == 'FPKM':
        f1 =  '/data2/myang9/SexAnnoDB/data/GDCdata/TCGA-{0}_expression_FPKM.csv'.format(cancer_type)
    if mol == 'SNP':
        f1 = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/UVM.hg38.gen'
    df1_sum = pd.read_csv(f1,sep=',',index_col=0,chunksize=100000)
    i = 0
    for df1 in df1_sum:
        i = i+1
        df1_col = [x[0:15] for x in df1.columns.values.tolist()]
        df1.columns = df1_col
        outdir = '/data2/myang9/SexAnnoDB/data/{0}_input/tmp/'.format(Ana_type)
        for key in ['male_tumor','female_tumor']:
            l = d_group[key]
            df1_tmp = df1[l]
            if df1_tmp.shape[1]!=0:
                outf = outdir + '{0}_{1}_{2}_{3}.txt'.format(cancer_type,key,Ana_type+'_'+mol,str(i))
                df1_tmp.to_csv(outf,sep='\t',header=0)

        if df1_female_tumor.shape[1]!=0 and True in df1_female_tumor.columns.duplicated():
            df1_female_tumor = df1_female_tumor.loc[:,~df1_female_tumor.columns.duplicated()]
        outf = outdir + '{0}_female_tumor_FPKM_{1}.txt'.format(cancer_type,str(i))
        df1_female_tumor.to_csv(outf,sep='\t',header=0)
        df1_female_normal = df1[l_female_normal]
        if df1_female_normal.shape[1]!=0 and True in df1_female_normal.columns.duplicated():
            df1_female_normal = df1_female_normal.loc[:,~df1_female_normal.columns.duplicated()]
        outf = outdir + '{0}_female_normal_FPKM_{1}.txt'.format(cancer_type,str(i))
        df1_female_normal.to_csv(outf,sep='\t')
        df1_male_tumor = df1[l_male_tumor]
        if df1_male_tumor.shape[1]!=0 and  True in df1_male_tumor.columns.duplicated():
            df1_male_tumor = df1_male_tumor.loc[:,~df1_male_tumor.columns.duplicated()]
        outf = outdir + '{0}_male_tumor_FPKM_{1}.txt'.format(cancer_type,str(i))
        df1_male_tumor.to_csv(outf,sep='\t')
        df1_male_normal = df1[l_male_normal]
        if df1_male_normal.shape[1]!=0 and  True in df1_male_normal.columns.duplicated():
            df1_male_normal = df1_male_normal.loc[:,~df1_male_normal.columns.duplicated()]
        outf = outdir + '{0}_male_normal_FPKM_{1}.txt'.format(cancer_type,str(i))
        df1_male_normal.to_csv(outf,sep='\t')

        df2_female_tumor = df2[l_female_tumor]
        if df2_female_tumor.shape[1]!=0 and  True in df2_female_tumor.columns.duplicated() :
            df2_female_tumor = df2_female_tumor.loc[:,~df2_female_tumor.columns.duplicated()]
        outf = outdir + '{0}_female_tumor_SNP.txt'.format(cancer_type)
        df2_female_tumor.to_csv(outf,sep='\t')
        df2_female_normal = df2[l_female_normal]
        if df2_female_normal.shape[1]!=0  and True in df2_female_normal.columns.duplicated():
            df2_female_normal = df2_female_normal.loc[:,~df2_female_normal.columns.duplicated()]
        outf = outdir + '{0}_female_normal_SNP.txt'.format(cancer_type)
        df2_female_normal.to_csv(outf,sep='\t')
        df2_male_tumor = df2[l_male_tumor]
        if df2_male_tumor.shape[1]!=0 and True in  df2_male_tumor.columns.duplicated():
            df2_male_tumor = df2_male_tumor.loc[:,~df2_male_tumor.columns.duplicated()]
        outf = outdir + '{0}_male_tumor_SNP.txt'.format(cancer_type)
        df2_male_tumor.to_csv(outf,sep='\t')
        df2_male_normal = df2[l_male_normal]
        if df2_male_normal.shape[1]!=0 and True in df2_male_normal.columns.duplicated():
            df2_male_normal = df2_male_normal[:,~df2_male_normal.columns.duplicated()]
        outf = outdir + '{0}_male_normal_SNP.txt'.format(cancer_type)
        df2_male_normal.to_csv(outf,sep='\t')
  
   

def eQTL_cmd():
    path = '/data2/myang9/SexAnnoDB/data/eQTL_input/'
    r = open('/data2/myang9/SexAnnoDB/output/eQTL_rn.sh','w')
    j = 0
    for inf in os.listdir(path):
        if '_FPKM.txt' not in inf:
            continue
        print(inf)
        cancer_type = inf.split('_')[0]
        group = '_'.join(inf.split('_')[1:3])
        print(cancer_type,group)
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/eQTL.r {0} {1}'.format(cancer_type,group)
        j = j +1
        r.write('nohup '+cmd+'  &\n')
        if j%10 ==0:
            r.write('wait;\n')
    

def sQTL_cmd():
    path = '/data2/myang9/SexAnnoDB/data/sQTL_input/'
    j = 0
    r = open('/data2/myang9/SexAnnoDB/output/sQTL_rn.sh','w')
    for inf in os.listdir(path):
        if '_PSI.txt' not in inf:
            continue
        print(inf)
        cancer_type = inf.split('_')[0]
        group = '_'.join(inf.split('_')[1:3])
        #print(cancer_type,group)
        outf = '/data2/myang9/SexAnnoDB/output/Network/sQTL/{0}_{1}_sQTL_cis.txt'.format(cancer_type,group)
        #if os.path.exists(outf):
        #    continue
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/sQTL.r {0} {1}'.format(cancer_type,group)
        j = j +1
        r.write('nohup '+cmd+'  &\n')
        if j%10 ==0:
            r.write('wait;\n')



def GetsQTLinput_multip(cancer_types):
    task_l = []
    task_l2= []
    d_Patientinfo = PA.PatientSubgroup()
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        #if os.path.exists('/data2/myang9/SexAnnoDB/data/sQTL_input/{0}_female_tumor_SNP.txt'.format(cancer_type)):
        #    continue
        else:
            l = (cancer_type,d_Patientinfo)
            task_l.append(l)
            #GetsQTLinput(cancer_type,d_Patientinfo)
            for sex in ['female','male']:
                l2 = (cancer_type,sex)
                task_l2.append(l2)
                #GetsQTLinput2(cancer_type,sex)
                #GetsQTLinput3(cancer_type,sex)
    #cores = multiprocessing.cpu_count()
    cores = 20
    pool = multiprocessing.Pool(processes=cores)
    pool.starmap(GetsQTLinput,task_l)
    #pool.starmap(GetsQTLinput2,task_l2)
    #pool.starmap(GetsQTLinput3,task_l2)

def GetsQTLinput(cancer_type,d_Patientinfo):
       # print(cancer_type)
        path = '/data2/myang9/SexAnnoDB/data/GDCdata/' 
        f1 = '/data2/myang9/DriverRBP/data/ASevents/{0}_exon_skip_ASmatrix.txt'.format(cancer_type)
        df1 = pd.read_csv(f1,sep='\t',index_col=0)
        f2 = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.gen'.format(cancer_type)
        df2 = pd.read_csv(f2,sep='\t',index_col=0)
        df1_col = [x[0:15] for x in df1.columns.values.tolist()]
        df2_col = [x[0:15] for x in df2.columns.values.tolist()]
        df1.columns = df1_col
        df2.columns = df2_col
        #print(df1,df2)
        patient_ids = [x for x in df1_col if x in df2_col]
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
        df1_female_tumor = df1[l_female_tumor]
        outdir = '/data2/myang9/SexAnnoDB/data/sQTL_input/'
        if df1_female_tumor.shape[1]!=0 and True in df1_female_tumor.columns.duplicated():
            df1_female_tumor = df1_female_tumor.loc[:,~df1_female_tumor.columns.duplicated()]
        outf = outdir + '{0}_female_tumor_PSI.txt'.format(cancer_type)
        df1_female_tumor.to_csv(outf,sep='\t')
        df1_female_normal = df1[l_female_normal]
        if df1_female_normal.shape[1]!=0 and True in df1_female_normal.columns.duplicated():
            df1_female_normal = df1_female_normal.loc[:,~df1_female_normal.columns.duplicated()]
        outf = outdir + '{0}_female_normal_PSI.txt'.format(cancer_type)
        df1_female_normal.to_csv(outf,sep='\t')
        df1_male_tumor = df1[l_male_tumor]
        if df1_male_tumor.shape[1]!=0 and  True in df1_male_tumor.columns.duplicated():
            df1_male_tumor = df1_male_tumor.loc[:,~df1_male_tumor.columns.duplicated()]
        outf = outdir + '{0}_male_tumor_PSI.txt'.format(cancer_type)
        df1_male_tumor.to_csv(outf,sep='\t')
        df1_male_normal = df1[l_male_normal]
        if df1_male_normal.shape[1]!=0 and  True in df1_male_normal.columns.duplicated():
            df1_male_normal = df1_male_normal.loc[:,~df1_male_normal.columns.duplicated()]
        outf = outdir + '{0}_male_normal_PSI.txt'.format(cancer_type)
        df1_male_normal.to_csv(outf,sep='\t')

        df2_female_tumor = df2[l_female_tumor]
        if df2_female_tumor.shape[1]!=0 and  True in df2_female_tumor.columns.duplicated() :
            df2_female_tumor = df2_female_tumor.loc[:,~df2_female_tumor.columns.duplicated()]
        outf = outdir + '{0}_female_tumor_SNP.txt'.format(cancer_type)
        df2_female_tumor.to_csv(outf,sep='\t')
        df2_female_normal = df2[l_female_normal]
        if df2_female_normal.shape[1]!=0  and True in df2_female_normal.columns.duplicated():
            df2_female_normal = df2_female_normal.loc[:,~df2_female_normal.columns.duplicated()]
        outf = outdir + '{0}_female_normal_SNP.txt'.format(cancer_type)
        df2_female_normal.to_csv(outf,sep='\t')
        df2_male_tumor = df2[l_male_tumor]
        if df2_male_tumor.shape[1]!=0 and True in  df2_male_tumor.columns.duplicated():
            df2_male_tumor = df2_male_tumor.loc[:,~df2_male_tumor.columns.duplicated()]
        outf = outdir + '{0}_male_tumor_SNP.txt'.format(cancer_type)
        df2_male_tumor.to_csv(outf,sep='\t')
        df2_male_normal = df2[l_male_normal]
        if df2_male_normal.shape[1]!=0 and True in df2_male_normal.columns.duplicated():
            df2_male_normal = df2_male_normal[:,~df2_male_normal.columns.duplicated()]
        outf = outdir + '{0}_male_normal_SNP.txt'.format(cancer_type)
        df2_male_normal.to_csv(outf,sep='\t')
        
def GetsQTLinput2(cancer_type,sex):
    psifile = '/data2/myang9/SexAnnoDB/data/sQTL_input/{0}_{1}_tumor_PSI.txt'.format(cancer_type,sex)
    df1 = pd.read_csv(psifile,sep='\t',index_col=0)
    df1 = df1.T
    l_index = []
    n = df1.shape[0]*0.1
    for index,row in df1.iteritems():
        tmp = df1[(df1[index]<0.9)&(df1[index]>0.1)]
        if tmp.shape[0] > n:
            l_index.append(index)
    df1 = df1[l_index]
    df1 = df1.T
    outf = '/data2/myang9/SexAnnoDB/data/sQTL_input/{0}_{1}_tumor_PSI_select.txt'.format(cancer_type,sex)
    df1.to_csv(outf,sep='\t')
    

def EX_position():
    f = open('/data2/myang9/DriverRBP/data/information/exon_skip_AS_information.txt')
    f.readline()
    r = open('/data2/myang9/SexAnnoDB/data/QTL_position/EX_position.txt','w')
    head = ['geneid','chr','s1','s2']
    r.write('\t'.join(head)+'\n')
    for line in f:
        line = line.strip().split()
        exid = line[0]
        chrom = line[3].split(':')[0]
        l = line[3].split(':')[1:]
        l = [int(x) for x in l]
        l.sort()
        newl = [exid,chrom,str(l[0]),str(l[-1])]
        #print(newl)
        r.write('\t'.join(newl)+'\n')

def SNP_matrix_idchange_multip():
    cancer_types = ['BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL']
    sexs = ['female','male']
    Types = ['tumor','normal']
    task_l = []
    for cancer_type in cancer_types:
        for sex in sexs:
            for Type in Types:
                l = (cancer_type,sex,Type)
                task_l.append(l)
    cores = multiprocessing.cpu_count()
    cores = int(cores-50)
    pool = multiprocessing.Pool(processes=cores)
    pool.starmap(SNP_matrix_idchange,task_l)

def SNP_matrix_idchange(cancer_type,sex,Type):
    path = '/data2/myang9/SexAnnoDB/data/eQTL_input/'
    snp_f = path + '{0}_{1}_{2}_SNP.txt'.format(cancer_type,sex,Type)
    print(snp_f)
    r = open(snp_f.replace('_SNP.txt','_SNP_rsid.txt'),'w')
    f = open(snp_f)
    head = f.readline()
    r.write(head)
    for line in f:
            line = line.strip().split()
            line[0] = line[0].split('|')[1].split(':')[0]
            r.write('\t'.join(line)+'\n')
            #print(line)


def SNP_hg19tohg38_multip():
    task_l = []
    cancer_types = ['BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL']
    for cancer_type in cancer_types:
        task = (cancer_type,'hg19')
        task_l.append(task)
        #print(task)
        #SNP_hg19tohg38_step3(cancer_type,'hg38')
    cores = 30
    pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(SNP_hg19tohg38_step1,task_l)
    #SNP_hg19tohg38_step2(cancer_types)
    pool.starmap(SNP_hg19tohg38_step3,task_l)

def SNP_hg19tohg38_step1(cancer_type,hg19):
    inf = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.location'.format(cancer_type)
    ouf = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg19.bed'.format(cancer_type)
    f = open(inf)
    head = f.readline()
    r = open(ouf,'w')
    for line in f:
        line = line.strip().split()
        newl = line[1:]+[str(int(line[-1])+1),line[0]]
        r.write('\t'.join(newl)+'\n')

def SNP_hg19tohg38_step2(cancer_types):
    chain = '/data2/myang9/library/human/liftOver/hg19ToHg38.over.chain'
    for cancer_type in cancer_types:
        hg19f = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg19.bed'.format(cancer_type)
        hg38f = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.bed'.format(cancer_type)
        unmap =  '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.unmap.bed'.format(cancer_type)
        cmd = 'liftOver {0} {1} {2} {3}'.format(hg19f,chain,hg38f,unmap)
        print(cmd)
        os.system('nohup ' + cmd+' &')

def SNP_hg19tohg38_step3(cancer_type,hg38loc):
    inf = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.bed'.format(cancer_type)
    outf = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.location'.format(cancer_type)
    f = open(inf)
    head = ['SNPID','Chromosome','Position']
    r = open(outf,'w')
    r.write('\t'.join(head)+'\n')
    d_loc = {}
    for line in f:
        line = line.strip().split()
        hg19id = line[-1]
        hg38id = line[0]+'|'+hg19id.split('|')[1].split(':')[0]+':'+line[1]+':'+hg19id.split('|')[1].split(':')[2]+':'+hg19id.split('|')[1].split(':')[3]
        d_loc[hg19id]=hg38id
        newl = [hg38id]+line[0:2]
        r.write('\t'.join(newl)+'\n')
    r.close()
    snpf= '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.gen'.format(cancer_type)
    f = open(snpf)
    head = f.readline().strip().split()
    newhead = [x[0:15] for x in head]
    snpf_hg38 = '/data2/myang9/SexAnnoDB/data/TCGA/SNP_impute/{0}.hg38.gen'.format(cancer_type)
    r = open(snpf_hg38,'w')
    r.write('\t'.join(newhead)+'\n')
    for line in f:
        line = line.strip().split()
        if line[0] not in d_loc:
            continue
        line[0] = d_loc[line[0]]
        r.write('\t'.join(line)+'\n')
    r.close()













