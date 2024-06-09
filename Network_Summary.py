import os,sys
import multiprocessing
from scipy import stats
import pandas as pd
import Proteindata_Ana as PA
import ASdata_anlaysis as AA
import multiprocessing
import Methyl_analysis as MA

def RBP_EX():
    d_exidchange = AA.EXid_gene()
    f = open('/data2/myang9/SexAnnoDB/output/Network/RBP_EX/RBP_exon_id_net_classification.txt')
    head = f.readline()
    r = open('/data2/myang9/SexAnnoDB/output/Network/SexAnnoDB_RBP_EX_Summary.txt','w')
    head = ['Index','GeneID','ESid','RBP','Turmol_Fale_effect_score','Turmol_Female_effect_score','Normal_male_effect_score','Normal_Female_effect_score','calssification','cancer_type']
    i = 0 
    r.write('\t'.join(head)+'\n')
    for line in f:
        if 'Other' in line:
            line = line.replace('Other','-')
        i = i + 1
        line = line.strip().split()
        exid = line[1]
        info = d_exidchange[exid]
        skipregion = line[-1]
        newl = [str(i),line[5],line[1],line[3]]+line[7:11]+[line[-2],line[-1]]
        r.write('\t'.join(newl)+'\n')

def RBP_ES_Summary(cancer_types):
    #RBP_ES_Summary_forGene_sub1()
    #RBP_ES_Summary_forGene_sub2()
    #RBP_ES_Summary_forGene_sub3()
    #RBP_ES_Summary_forGene_plot1()
    #RBP_ES_Summary_forGene_plot2()
    RBP_ES_Summary_forRBP_sub1()
    RBP_ES_Summary_forRBP_plot1()
    RBP_ES_Summary_forRBP_plot2()
    #RBP_ES_Draw_forGene()
    #RBP_ES_Summary_forCancerType()
    #RBP_ES_Draw_forCancer(cancer_types)




def RBP_ES_Summary_forGene_sub1():
    d_exidchange = AA.EXid_gene()
    path = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/zzzz_diff_edge_filter/'
    infs = os.listdir(path)
    i = 0 
    outf = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary.txt'
    r = open(outf,'w')
    head = ['Index','WebGeneID','Cancer','RBP','TargetES','Male_Edge_score','Male_Edge_threshold','Female_Edge_score','Female_Edge_threshold','Type']
    r.write('\t'.join(head)+'\n')
    for inf in infs:
        cancer = inf.split('_')[1]
        inf = path + inf
        f = open(inf)
        head = f.readline()
        for line in f:
            line = line.strip().split()
            exid = line[1]
            if exid not in d_exidchange:
                #print(exid)
                continue
            exinfo = d_exidchange[exid]
            Targetgeneid = exinfo[0]
            TargetGenesym = exinfo[0]
            i = i + 1
            wl = [str(i),Targetgeneid,cancer,line[0],exid] + line[2:4]+line[5:7]
            male_t = float(line[3])
            female_t = float(line[6])
            if male_t > 0.98 and female_t < 0.98:
                Type = 'Male-biased'
            elif male_t < 0.98 and female_t > 0.98:
                Type = 'Female-biased'
            else:
                print(wl)
            wl.append(Type)
            r.write('\t'.join(wl)+'\n')

def RBP_ES_Summary_forGene_sub2():
    inf = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary.txt'
    f = open(inf)
    d_pairs = {}
    d_SymToID = MA.GencodeV22_SymToID()
    head = f.readline()
    for line in f:
        line = line.strip().split()
        cancer = line[2]
        genesym = line[3]
        if genesym == 'BRUNOL4':
            geneid = 'ENSG00000101489'
        elif genesym == 'BRUNOL5':
            geneid = 'ENSG00000161082'
        elif genesym == 'BRUNOL6':
            geneid = 'ENSG00000140488'
        elif genesym == 'Fusip1':
            geneid = 'ENSG00000188529'
        elif genesym == 'hnRNPK':
            geneid = 'ENSG00000165119'
        elif genesym == 'hnRNPLL':
            geneid = 'ENSG00000143889'
        elif genesym == 'STAR-PAP':
            geneid = 'ENSG00000149016'
        elif genesym == 'YB-1':
            geneid = 'ENSG00000065978'
        else:
            geneid = d_SymToID[genesym]
        exid = line[4]
        pairs = geneid + '|' + exid
        if cancer not in d_pairs:
            d_pairs[cancer] = [pairs]
        else:
            d_pairs[cancer].append(pairs)
    f.close()
    exp_path = '/data2/myang9/SexAnnoDB/data/GDCdata/'
    AS_path = '/data2/myang9/DriverRBP/data/ASevents/'
    d_Patientinfo = PA.PatientSubgroup()
    l_task = []
    cores = 30
    pool = multiprocessing.Pool(processes=cores)
    for cancer in d_pairs:
        task = (cancer,exp_path,AS_path,d_Patientinfo,d_pairs)
        #RBP_ES_cor(cancer,exp_path,AS_path,d_Patientinfo,d_pairs)
        l_task.append(task)
    pool.starmap(RBP_ES_cor,l_task)

def RBP_ES_cor(cancer,exp_path,AS_path,d_Patientinfo,d_pairs):
        as_f = AS_path + cancer + '_exon_skip_ASmatrix.txt'
        exp_f = exp_path + "TCGA-{0}_expression_FPKM.csv".format(cancer)
        df_as = pd.read_csv(as_f,sep='\t',index_col=0)
        df_exp = pd.read_csv(exp_f,sep=',',index_col=0)
        patient_ids_as =[x[0:15] for x in df_as.columns.values.tolist()]
        df_as.columns = patient_ids_as
        patient_ids_exp = [x[0:15] for x in df_exp.columns.values.tolist()]
        df_exp.columns = patient_ids_exp
        patient_ids = [pid for pid in patient_ids_as if pid in patient_ids_exp]

        if cancer+'_turmol_female' not in d_Patientinfo:
            l_female_tumor = []
        else: 
            l_female_tumor = [x for x in d_Patientinfo[cancer+'_turmol_female'] if x in patient_ids]
        if cancer+'_normal_female' not in d_Patientinfo:
            l_female_normal = []
        else:
            l_female_normal = [x for x in d_Patientinfo[cancer+'_normal_female'] if x in patient_ids]
        if cancer+'_turmol_male' not in d_Patientinfo:
            l_male_tumor = []
        else:
            l_male_tumor = [x for x in  d_Patientinfo[cancer +'_turmol_male'] if x in patient_ids]
        outf = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_{0}_RBP_ES_correlation.txt'.format(cancer)
        head = ['cancer','group','RBP','ES','cor','p']
        r = open(outf,'w')
        r.write('\t'.join(head)+'\n')
        d_overlap = {}
        d_overlap['female_tumor'] = l_female_tumor
        d_overlap['male_tumor'] = l_male_tumor
        for group in d_overlap:
            p_list = d_overlap[group]
            df_exp_select = df_exp[p_list]
            df_as_select = df_as[p_list]
            for pairs in d_pairs[cancer]:
                gene = pairs.split('|')[0]
                es = pairs.split('|')[1]
                if gene not in df_exp_select.T.columns:
                    continue
                l1 = df_exp_select.T[gene].sort_index().values.tolist()
                l2 = df_as_select.T[es].sort_index().values.tolist()
                newd = pd.DataFrame([l1,l2])
                newd = newd.T
                newd.columns = ['rbp','as']
                newd = newd.dropna()
                l1 = newd['as']
                l2 = newd['rbp']
                cor,p = stats.stats.pearsonr(l1,l2)
                wl = [cancer,group,gene,es,cor,p]
                if str(cor) == 'nan':
                    continue
                if abs(cor) < 0.3 or p>0.05:
                    continue
                wl = [str(x) for x in wl]
                r.write('\t'.join(wl)+'\n')

def RBP_ES_Summary_forRBP_sub1():
    d_SymToID = MA.GencodeV22_SymToID()
    f = open('/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary.txt')
    head = f.readline().strip().split()
    r = open('/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary2.txt','w')
    #print(head)
    head = ['Index', 'WebGeneID','ESindexID', 'Cancer', 'RBP', 'TargetES', 'Male_Edge_score', 'Male_Edge_threshold', 'Female_Edge_score', 'Female_Edge_threshold', 'Type']
    r.write('\t'.join(head)+'\n')
    d_count = {}
    for line in f:
        line = line.strip().split()
        #print(line)
        genesym = line[3]
        if genesym == 'BRUNOL4':
            geneid = 'ENSG00000101489'
        elif genesym == 'BRUNOL5':
            geneid = 'ENSG00000161082'
        elif genesym == 'BRUNOL6':
            geneid = 'ENSG00000140488'
        elif genesym == 'Fusip1':
            geneid = 'ENSG00000188529'
        elif genesym == 'hnRNPK':
            geneid = 'ENSG00000165119'
        elif genesym == 'hnRNPLL':
            geneid = 'ENSG00000143889'
        elif genesym == 'STAR-PAP':
            geneid = 'ENSG00000149016'
        elif genesym == 'YB-1':
            geneid = 'ENSG00000065978'
        elif genesym == 'TARDBP':
            geneid = 'ENSG00000120948'
        else:
            geneid = d_SymToID[genesym]
        #print(line)
        wl = [line[0],geneid,line[1]]+line[2:]
        r.write('\t'.join(wl)+'\n')
        key = line[1] + line[2]
        if key in d_count:
            d_count[key].append(wl[5])
        else:
            d_count[key] = [wl[5]]
    for key in d_count:
        if len(d_count[key]) > 10:
            print(len(d_count[key]),key)


def RBP_ES_Summary_forGene_plot1():
    inf = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary.txt'
    f = open(inf)
    d_out = {}
    head = f.readline()
    #head = ['GeneSym','GeneID','Score','Threshold','NetName']
    for line in f:
        line = line.strip().split()
        #print(line)
        wl1 = line[1:7] + [line[-1]]
        wl2 = line[1:5] + line[7:]
        #print(len(wl1),wl1)
        #print(len(wl2),wl2)
        #print('\n\n')
        geneid = line[1]
        if geneid in d_out:
            d_out[geneid]= d_out[geneid] + [wl1,wl2]
        else:
            d_out[geneid] = [wl1,wl2]
    outpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFigures/files/'
    for geneid  in d_out:
        outf = outpath + geneid+'_RBP_ES.txt'
        r = open(outf,'w')
        head = ['GeneID','Cancer','RBP','ESID','Score','Threshold','NetName']
        r.write('\t'.join(head)+'\n')
        for line in d_out[geneid]:
            r.write('\t'.join(line)+'\n')
            #print(line)
        r.close()


def RBP_ES_Summary_forGene_plot2():
    r = open('RBP_ES_Draw.sh','w')
    fpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFigures/files/'
    infs = os.listdir(fpath)
    i = 0
    r.write("conda activate r4-base\n")

    for inf in infs:
        geneid = inf.split('_')[0]
        i = i + 1
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script_0717/script/Draw/RBP-ES-pairs.R {0} {1}'.format(inf,geneid)
        if i%10 == 0 :
            r.write('wait;\n')
        r.write('nohup '+cmd+'  &\n')
    

def RBP_ES_Summary_forRBP_plot1():
    d_idchange = GeneIDtoSymbol()
    inf = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary2.txt'
    f = open(inf)
    head = f.readline()
    head = ['from','to','Score','Type','Cancer']
    d_out = {}
    for line in f:
        line = line.strip().split()
        #print(line)
        esgene = line[2]
        if esgene not in d_idchange:
            continue
        esgenesym = d_idchange[esgene]
        if line[-1] == 'Male-biased':
            wl = [line[4],esgenesym+':'+line[5],line[6],line[-1],line[3]]
        else:
            wl = [line[4],esgenesym+':'+line[5],line[8],line[-1],line[3]]
        key = line[1]
        #print(wl,'\n\n\n')
        if key in d_out:
            d_out[key] = d_out[key] + [wl]
        else:
            d_out[key] = [wl]
    outpath = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFiguresRBP/files/'
    for key in d_out:
        r = open(outpath+'{0}_RBP.txt'.format(key),'w')
        head = ['from','to','Score','Type','Cancer']
        r.write('\t'.join(head)+'\n')
        for line in d_out[key]:
            r.write('\t'.join(line)+'\n')
        r.close()
def  RBP_ES_Summary_forRBP_plot2():
    path = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/WebFiguresRBP/files/'
    infs = os.listdir(path)
    i = 0
    r = open('RBP-ES-pairs_forRBP.sh','w')
    for inf in infs:
        geneid = inf.split('_')[0]
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script_0717/script/Draw/RBP-ES-pairs_forRBP.R  {0} {1}'.format(inf,geneid)
        #print(cmd)
        r.write('nohup '+cmd+' &\n')
        i = i + 1
        if i%20==0:
            r.write('wait;\n')


def RBP_ES_Draw_forCancer(cancer_types):
    path = '/data2/myang9/SexAnnoDB/output/Network/RBP_ES/'
    for cancer_type in cancer_types:
        inf = path + '{0}_RBP_ES_Summary.txt'.format(cancer_type)
        if os.path.exists(inf) == False:
            continue
        '''
        r = open(path+'{0}_RBP_ES_Summary_Top20.txt'.format(cancer_type),'w')
        df =  pd.read_csv(inf,sep = '\t')
        df_T = df.loc[df['Rtype']=='T',]
        df_F = df.loc[df['Rtype']=='F',]
        df_T['R'] = df_T['from'].str.cat(df_T['to'],sep='-')
        df_F['R'] = df_F['from'].str.cat(df_F['to'],sep='-')
        df_T = df_T.sort_values("score",ascending=False)
        df_F = df_F.sort_values("score",ascending=False)
        df_T = df_T[:20]
        df_F = df_F[:20]
        df_T = df_T[['R','score']]
        df_F = df_F[['R','score']]
        df_T.index = [str(i) for i in range(1,21)]
        df_F.index = [str(i) for i in range(1,21)]
        df_out = pd.concat([df_T,df_F],axis=1,join='inner')
        #wl = [','.join(l_T),','.join(l_F)]
        '''
        outf = '/data2/myang9/SexAnnoDB/output/Figures/RBP_ES/figure_Cancer/{0}_RBP_ES_Summary.svg'.format(cancer_type)
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/network_RBP_ES_cancer.r {0} {1} {2}'.format(inf,outf,cancer_type)
        os.system(cmd)
        print(cmd)

def RBP_ES_Summary_forGene():
    d_exidchange = AA.EXid_gene()
    f = open('/data2/myang9/SexAnnoDB/output/Network/RBP_ES/RBP_exon_id_net_classification.txt')
    head = f.readline()
    d_out = {}
    d_score_F = {}
    d_score_T = {}
    for line in f:
        line = line.strip().split()
        exid = line[1]
        geneid = line[5]
        RBPname = line[3]
        Rtype = line[11]
        cancer_type = line[-1]
        if Rtype == 'NA':
            continue
        key = geneid+'_'+exid+'_'+cancer_type
        if key in d_out:
            d_out[key][Rtype].append(RBPname)
        else:
            d_out[key]={}
            d_out[key]['F']=[]
            d_out[key]['T']=[]
            d_out[key][Rtype].append(RBPname)
        if Rtype == 'T':
            score = line[7]
            wl = [RBPname,exid,score,cancer_type]
            key2 = geneid# +'_'+exid
            if key2 in d_score_T:
                d_score_T[key2].append(wl)
            else:
                d_score_T[key2] = [wl]
        else:
            score = line[8]
            wl = [RBPname,exid,score,cancer_type]
            key2 = geneid #+'_'+exid
            if key2 in d_score_F:
                d_score_F[key2].append(wl)
            else:
                d_score_F[key2] = [wl]
    r = open('/data2/myang9/SexAnnoDB/output/Network/RBP_ES/SexAnnoDB_RBP_ES_Summary.txt','w')
    head = ['Index','GeneID','ESID','Cancer_type','Male_biased','Female_biased']
    r.write('\t'.join(head)+'\n')
    i = 0
    for key in d_out:
        i = i+1
        if len(d_out[key]['T'])==0:
            d_out[key]['T'] = ['-']
        if len(d_out[key]['F']) == 0:
            d_out[key]['F'] = ['-']
        newl =[str(i)]+[key.split('_')[0],'_'.join(key.split('_')[1:4]),key.split('_')[4]]+[','.join(d_out[key]['T'])]+[','.join(d_out[key]['F'])]
        r.write('\t'.join(newl)+'\n')
    for key2 in d_score_T:
        r = open('/data2/myang9/SexAnnoDB/output/Figures/RBP_ES/file/{0}_RBP_ES_male.txt'.format(key2),'w')
        head = ['from','to','score','cancer_type']
        r.write('\t'.join(head)+'\n')
        for line in d_score_T[key2]:
            r.write('\t'.join(line)+'\n')
        r.close()
    for key2 in d_score_F:
        r = open('/data2/myang9/SexAnnoDB/output/Figures/RBP_ES/file/{0}_RBP_ES_female.txt'.format(key2),'w')
        head = ['from','to','score','cancer_type']
        r.write('\t'.join(head)+'\n')
        for line in d_score_F[key2]:
            r.write('\t'.join(line)+'\n')
        r.close()

def RBP_ES_Draw_forGene():
    r = open('RBP_ES_Draw.sh','w')
    path = '/data2/myang9/SexAnnoDB/output/Figures/RBP_ES/file/'
    i = 0
    for inf in os.listdir(path):
        if '_male' not in inf:
            continue
        i = i + 1
        inf_male = path + inf
        inf_female = path + inf.replace('male','female')
        outf = path + inf.replace('male.txt','network.svg').replace('file','figure')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/network_RBP_ES.r {0} {1} {2}'.format(inf_male,inf_female,outf)
        r.write('nohup '+cmd+' &\n')
        if i%10 ==0:
            r.write('wait;\n')

def TF_Gene_Summary(cancer_types):
    #TF_Gene_Summary_forGene()
    #TF_Gene_draw_forGene()
    #TF_Gene_Summary_forCancer()
    #TF_Gene_Draw_forCancer(cancer_types)
    #TF_Gene_Summary_forTF_sub1()   
    TF_Gene_Summary_forTF_sub2()

def TF_Gene_Summary_forCancer():
    f = open('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/TF_gene_net_classification.txt')
    head = f.readline()
    d_net = {}
    for line in f:
        line = line.strip().split()
        genename = line[4]
        geneid = line[1]
        TFname = line[2]
        cancer_type = line[-1]
        Rtype = line[10]
        if Rtype == 'NA':
            continue
        key = geneid+'_'+cancer_type
        if Rtype == 'T':
            score = line[6]
        if Rtype == 'F':
            score = line[7]
        wl = [TFname,genename,score,Rtype]
        if cancer_type in d_net:
            d_net[cancer_type].append(wl)
        else:
            d_net[cancer_type] = [wl]
    f.close()
    for cancer_type in d_net:
        r = open('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/{0}_TF_Gene_Summary.txt'.format(cancer_type),'w')
        head = ['from','to','score','Rtype']
        r.write('\t'.join(head)+'\n')
        for line in d_net[cancer_type]:
            r.write('\t'.join(line)+'\n')
        r.close()
    
def TF_Gene_Summary_forGene():
    pass
    #TF_Gene_Summary_forGene_sub1()
    #TF_Gene_Summary_forGene_plot1()
    #TF_Gene_Summary_forGene_plot2()


def TF_Gene_Summary_forGene_sub1():
    d_idchange = GeneIDtoSymbol()
    path = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/tmp2024/'
    infs = [x for x in os.listdir(path) if '.txt' in x]
    head = ['Index','WebGeneID','Cancer','TF','TargetGene','Male_Edge_score','Male_Edge_threshold','Female_Edge_score','Female_Edge_threshold','Type']
    i = 0 
    r = open('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/SexAnnoDB_TF_Gene_Summary.txt','w')
    r.write('\t'.join(head)+'\n')
    for inf in infs:
        cancer = inf.split('_')[0]
        #print(cancer,inf)
        f = open(path+inf)
        head = f.readline()
        for line in f:
            line = line.strip().split()
            #print(line)
            i = i + 1
            if line[1] not in d_idchange:
                print(line[1])
                continue
            if float(line[3]) > 0.98 and float(line[6]) <= 0.98:
                label = 'Male-biased'
            elif float(line[3]) <= 0.98 and float(line[6]) > 0.98:
                label = 'Female-biased'
            else:
                print(line)
            wl = [str(i),line[1],cancer,line[0],d_idchange[line[1]]] + line[2:4]+line[5:7] + [label]
            r.write('\t'.join(wl)+'\n')

def TF_Gene_Summary_forGene_plot1():
    d_idchange = GeneIDtoSymbol()
    path = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/tmp2024/'
    d_out = {}
    infs = [x for x in os.listdir(path) if '.txt' in x]
    for inf in infs:
        cancer = inf.split('_')[0]
        inf = path+inf
        f = open(inf)
        head = f.readline()
        head = ['GeneSym','GeneID','Score','Threshold','NetName']
        for line in f:
            line = line.strip().split()
            wl1= line[0:5] + [cancer]
            wl2 = line[0:2] + line[5:] + [cancer]
            geneid = line[1]
            if geneid not in d_idchange:
                continue
            if geneid in d_out:
                d_out[geneid]= d_out[geneid] + [wl1,wl2]
            else:
                d_out[geneid] = [wl1,wl2]
        f.close()
    outpath = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/files/'
    for geneid  in d_out:
        #print(geneid)
        outf = outpath + geneid+'_TF-Gene.txt'
        r = open(outf,'w')
        head = ['TF','GeneID','Score','Threshold','NetName','cancer']
        r.write('\t'.join(head)+'\n')
        for line in d_out[geneid]:
            #print(line)
            r.write('\t'.join(line)+'\n')
        r.close()

def TF_Gene_Summary_forGene_plot2():
    r = open('TF_Gene_Draw.sh','w')
    r.write('conda activate r4-base')
    fpath = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFigures/files/'
    infs = os.listdir(fpath)
    i = 0 
    for inf in infs:
        geneid = inf.split('_')[0]
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script_0717/script/Draw/TF-Gene-pairs.R {0} {1}'.format(inf,geneid)
        #print(cmd)
        i = i + 1
        r.write('nohup '+ cmd + '&\n')
        if i%30 == 0:
            r.write('wait;\n')
    r.close()

def TF_Gene_Summary_forTF_sub1():
    d_SymToID = MA.GencodeV22_SymToID()
    inf = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/SexAnnoDB_TF_Gene_Summary.txt'
    f = open(inf)
    head = f.readline().strip().split()
    print(head)
    r = open('/data2/myang9/SexAnnoDB/output/Network/TF_Gene/SexAnnoDB_TF_Gene_Summary2.txt','w')
    r.write('\t'.join(head)+'\n')
    d_out  = {}
    head = ['from','to','Score','Cancer','Type']
    for line in f:
        line = line.strip().split()
        genesym = line[3]
        if genesym in d_SymToID:
            geneid = d_SymToID[line[3]]
        else:
            if genesym=='ZSCAN5':
                geneid = 'ENSG00000131848'
        #print(line)
        line[1] = geneid
        #print(line,'\n\n\n')
        r.write('\t'.join(line)+'\n')
        #print(line)
        if line[-1] == 'Male-biased':
            wl = [line[3],line[4],line[5],line[2],line[-1]]
        else:
            wl = [line[3],line[4],line[7],line[2],line[-1]]
        #print(wl,'\n\n')
        if line[1] in d_out:
            d_out[line[1]].append(wl)
        else:
            d_out[line[1]] = [wl]
    outpath = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFiguresTF/files/'
    head = ['from','to','Score','Cancer','Type']
    for key in d_out:
        r = open(outpath+'{0}_TF-Gene_forTF.txt'.format(key),'w')
        r.write('\t'.join(head)+'\n')
        for line in d_out[key]:
            r.write('\t'.join(line)+'\n')
        f.close()
def TF_Gene_Summary_forTF_sub2():
    path = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/WebFiguresTF/files/'
    infs = os.listdir(path)
    r = open('TF-Gene_forTF_plot.txt','w')
    i = 0 
    for inf in infs:
        geneid = inf.split('_')[0]
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script_0717/script/Draw/TF-Gene-pairs_forTF.R {0} {1}'.format(inf,geneid)
        i = i + 1
        r.write('nohup  '+cmd+'  &\n')
        
        if i%30 == 0:
            r.write('wait;\n')


def TF_Gene_Draw_forCancer(cancer_types):
    path = '/data2/myang9/SexAnnoDB/output/Network/TF_Gene/'
    for cancer_type in cancer_types:
        inf = path + '{0}_TF_Gene_Summary.txt'.format(cancer_type)
        if os.path.exists(inf) == False:
            continue
        outf = '/data2/myang9/SexAnnoDB/output/Figures/TF_Gene/figure_Cancer/{0}_TF_Gene_Summary.svg'.format(cancer_type)
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/network_TF_Gene_cancer.r {0} {1} {2}'.format(inf,outf,cancer_type)
        os.system(cmd)
        print(cmd)
def ceRNA():
    f = open('/data2/myang9/SexAnnoDB/output/Network/lncRNA_miRNA_coding_gene_cor_0.3_add_index.txt')
    head = f.readline().strip().split('\t')
    r = open('/data2/myang9/SexAnnoDB/output/Network/SexAnnoDB_ceRNA_Summary.txt','w')
    head = ['Index','gene_id','ceRNA(lncRNA-miRNA-mRNA)','Group','cancer_type']
    r.write('\t'.join(head)+'\n')
    for line in f:
        line = line.strip().split('\t')
        ceRNA = ','.join([line[4],line[7],line[5]])
        newl = line[0:2]+[ceRNA,line[-3],line[-4]]
        r.write('\t'.join(newl)+'\n')

def GeneIDtoSymbol():
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt')
    d_idchange = {}
    for line in f:
        line = line.strip().split()
        d_idchange[line[1]] = line[2]
    return d_idchange

def eQTL_Summary(cancer_types):
    task_l = []
    cancer_types = ['BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC']
    for cancer_type in cancer_types:
        pass
        l = (cancer_type,'eQTL')
        task_l.append(l)
        #eQTL_Summary_1(cancer_type,'eQTL')
    import multiprocessing
    cores = 20
    #pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(eQTL_Summary_1,task_l)
    eQTL_Summary_2(cancer_types)
    #eQTL_Summary_3(cancer_types)
def eQTL_Summary_3(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTL/SexAnnoDB_eQTL_summary2.txt')
    head = f.readline()
    print(head)
    d_male = {}
    d_female = {}
    for line in f:
        line = line.strip().split()
        geneid = line[1]
        genename = line[2]
        cancer_type = line[-1]
        key = '|'.join([cancer_type,geneid,genename])
        if line[6] != '-' and line[8] == '-':
            if key in d_male:
                d_male[key].append(line[3])
            else:
                d_male[key]=[line[3]]
        if line[6] == '-' and line[8] != '-':
            if key in d_female:
                d_female[key].append(line[3])
            else:
                d_female[key]=[line[3]]
    r = open('/data2/myang9/SexAnnoDB/output/Network/eQTL/SexAnnoDB_eQTL_summary2_female.txt','w')
    i = 0
    head = ['Index','Cancer_type','GeneID','GeneName','SNPID']
    r.write('\t'.join(head)+'\n')
    for key in d_female:
        i = i+1
        wl = [str(i)]+key.split('|')+d_female[key]
        r.write('\t'.join(wl)+'\n')
    r.close()
    r = open('/data2/myang9/SexAnnoDB/output/Network/eQTL/SexAnnoDB_eQTL_summary2_male.txt','w')
    i = 0
    r.write('\t'.join(head)+'\n')
    for key in d_male:
        i = i+1
        wl = [str(i)]+key.split('|')+d_male[key]
        r.write('\t'.join(wl)+'\n')
    r.close()

def eQTL_Summary_2(cancer_types):
    import Methyl_analysis as MA
    d_Stru = MA.GencodeV22_strucutre()
    path = '/data2/myang9/SexAnnoDB/output/Network/eQTL/'
    r = open(path+'SexAnnoDB_eQTL_summary.txt','w')
    r2 = open(path+'SexAnnoDB_eQTL_summary2.txt','w')
    head = ['Index','GeneID','GeneSymbol','SNP_ID','SNPinfo','Position_to_Gene','Male.effect','Male.FDR','Female.effect','Female.FDR','Sex_biased_gene_in_cancer','cancer_type']
    r.write('\t'.join(head)+'\n')
    i = 0
    r2.write('\t'.join(head)+'\n')
    for sumf in os.listdir(path):
        if '_sQTM_summary.txt' not in sumf:
            continue
        cancer_type = sumf.split('_')[0]
        sumf = path + sumf
        f = open('/data2/myang9/SexAnnoDB/output/Diff_analysis/CodingGene/Summary/SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type))
        f.readline()
        l_ts = [line.strip().split()[0] for line in f]
        f = open(sumf)
        head = f.readline()
        for line in f:
            line = line.strip().split()
            key = line[0]
            pos = int(line[3].split(':')[1])
            gene_pos= d_Stru[key]
            pos_l = []
            for g_pos in gene_pos:
                for position in gene_pos[g_pos]:
                    s = int(position.split(":")[1].split('-')[0])
                    e = int(position.split(":")[1].split('-')[1])
                    #print(pos,position)
                    if pos>=s and pos<=e:
                        pos_l.append(g_pos)
            if len(pos_l)==0:
                pos_l.append('-')
            pos_l = list(set(pos_l))
            i = i+1
            if key in l_ts:
                anno = 'Yes'
            else:
                anno = '-'
            line = [str(i)]+line[0:4]+[','.join(pos_l)]+line[4:]+[anno,cancer_type]
            if pos_l[0] != '-':
                r.write('\t'.join(line)+'\n')
                if float(line[7])<0.0001 and float(line[9])<0.0001:
                    r2.write('\t'.join(line)+'\n')
                elif float(line[7])<0.0001 and float(line[9])>0.05:
                    #if abs(float(line[8]))>0.01:
                    #    continue
                    r2.write('\t'.join(line)+'\n')
                elif float(line[9])<0.0001 and float(line[7])>0.05:
                    #if abs(float(line[6]))>0.01:
                    #    continue
                    r2.write('\t'.join(line)+'\n')

def eQTL_Summary_1(cancer_type,Analysis):
    d_idchange = GeneIDtoSymbol()
    print(cancer_type,Analysis)
    pos_f = '/data2/myang9/SexAnnoDB/data/QTL_position/{0}_female_tumor_SNP_pos_hg38.txt'.format(cancer_type)
    d_pos = {}
    f = open(pos_f)
    for line in f:
        line = line.strip().split()
        snpid = line[0]
        pos = line[1]+':'+line[2]
        d_pos[snpid] = pos
    f.close()

    path = '/data2/myang9/SexAnnoDB/output/Network/eQTL/tmp/'
    female_f = path + '{0}_female_tumor_eQTL_cis.txt'.format(cancer_type)
    male_f = path + '{0}_male_tumor_eQTL_cis.txt'.format(cancer_type)
    if os.path.exists(male_f) == False:
        return 'OK'
    f = open(male_f)
    f.readline()
    d_out = {}
    for line in f:
        line = line.strip().split()
        key = line[0]+'\t'+line[1]
        beta = line[2]
        d_out[key] = [line[2]+'p'+line[-1],'-p-']
    f.close()
    
    f = open(female_f)
    f.readline()
    for line in f:
        line = line.strip().split()
        key = line[0]+'\t'+line[1]
        if key in d_out:
            d_out[key][1] = line[2]+'p'+line[-1]
        else:
            d_out[key] = ['-p-',line[2]+'p'+line[-1]]
    f.close()
    i = 0
    Outdir = '/data2/myang9/SexAnnoDB/output/Network/eQTL/'
    outf = Outdir + '{0}_sQTM_summary.txt'.format(cancer_type)
    r = open(outf,'w')
    head = ['GeneID','GeneSymbol','SNP_ID','SNPinfo','Male.effect','Male.FDR','Female.effect','Female.FDR']
    r.write('\t'.join(head)+'\n')
    for key in d_out:
        if 'rs' not in key:
            continue
        snpinfo = key.split('\t')[0]
        geneid = key.split('\t')[1]
        if geneid not in d_idchange:
            continue
        rsid = snpinfo.split('|')[1].split(':')[0]
        ref = snpinfo.split('|')[1].split(':')[2]
        alt = snpinfo.split('|')[1].split(':')[3]        
        i = i +1
        wl = [geneid,d_idchange[geneid],rsid,d_pos[snpinfo]+':'+ref+':'+alt]+d_out[key][0].split('p')+d_out[key][1].split('p')
        if '-' in wl:
            continue
        if float(wl[5])<0.05  and float(wl[7])>0.1:
            r.write('\t'.join(wl)+'\n')
        elif float(wl[5])>0.1 and float(wl[7])<0.05:
            r.write('\t'.join(wl)+'\n')
        elif float(wl[5])<0.05 and float(wl[7])<0.05:
            if float(wl[5])>0 and float(wl[7])<0:
                r.write('\t'.join(wl)+'\n')
            elif float(wl[5])<0 and float(wl[7])>0:
                r.write('\t'.join(wl)+'\n')

def sQTL_Summary(cancer_types):
    cancer_types = ['GBM','KIRP','PAAD','PCPG','BLCA','LIHC','STAD','SARC','HNSC','COAD','LUSC','KIRC','LGG','LUAD','THCA']   
    d_exidchange = AA.EXid_gene()
    l_task = []
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        pass
        d_pos = GetSNPPos(cancer_type)
        l = (cancer_type,d_pos,d_exidchange)
        l_task.append(l)
        #sQTL_Summary_1(cancer_type,d_pos,d_exidchange)
    cores = 20
    #pool = multiprocessing.Pool(processes=cores)
    #pool.starmap(sQTL_Summary_1,l_task)
    #sQTL_Summary_2(cancer_types)
    #sQTL_Summary_3(cancer_types)
    sQTL_Summary_4(cancer_types)

def sQTL_Summary_4(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary3.txt')
    head = f.readline()
    d_male = {}
    d_female = {}
    d_opp = {}
    for line in f:
        line = line.strip().split('\t')
        #print(line)
        geneid = line[1]
        genename = line[2]
        cancer_type = line[-1]
        exid = line[3]
        key = key = '|'.join([cancer_type,geneid,genename,exid])
        if float(line[11]) <0.0001  and float(line[13]) >0.05 :
            if key in d_male:
                d_male[key].append(line[7])
            else:
                d_male[key]=[line[7]]
        elif float(line[11]) >0.05  and float(line[12]) <0.0001:
            if key  in d_female:
                d_female[key].append(line[7])
            else:
                d_female[key]=[line[7]]
        elif float(line[11])<0.0001 and float(line[12]) <0.0001:
            if key  in d_opp:
                d_opp[key].append(line[7])
            else:
                d_opp[key]= [line[7]]
    r = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_summary3_CancerPage.txt','w')
    i = 0
    head = ['Index','Cancer_type','GeneID','GeneName','ESID','SNPID','Type']
    r.write('\t'.join(head)+'\n')
    for key in d_female:
        i = i+1
        wl = [str(i)]+key.split('|')+d_female[key] + ['model1']
        r.write('\t'.join(wl)+'\n')
    for key in d_male:
        i = i+1
        wl = [str(i)]+key.split('|')+d_male[key]+['model2']
        r.write('\t'.join(wl)+'\n')
    for key in d_opp:
        i = i+1
        wl = [str(i)]+key.split('|')+d_opp[key]+['model3&4']
        r.write('\t'.join(wl)+'\n')
    r.close()

def GetSNPPos(cancer_type):
    pos_f = '/data2/myang9/SexAnnoDB/data/QTL_position/{0}_female_tumor_SNP_pos_hg38.txt'.format(cancer_type)
    d_pos = {}
    f = open(pos_f)
    for line in f:
        line = line.strip().split()
        snpid = line[0]
        pos = line[1]+':'+line[2]
        d_pos[snpid] = pos
    f.close()
    return d_pos

def sQTL_Summary_1(cancer_type,d_pos,d_exidchange):
    print(cancer_type)
    diff_f = '/data2/myang9/SexAnnoDB/output/Diff_analysis/ASevents/Summary/SexAnnoDB_Tumor_specific_sex_biased_PC_{0}.txt'.format(cancer_type)
    f = open(diff_f)
    diffES = [line.strip().split()[0] for line in f]
    f.close()
    path = '/data2/myang9/SexAnnoDB/output/Network/sQTL/tmp/'
    female_f = path + '{0}_female_tumor_sQTL_cis.txt'.format(cancer_type)
    male_f = path + '{0}_male_tumor_sQTL_cis.txt'.format(cancer_type)
    f = open(male_f)
    f.readline()
    d_out = {}
    for line in f:
        if 'Inf' in line:
            continue
        line = line.strip().split()
        key = line[0]+'\t'+line[1]
        beta = line[2]+'p'+line[-1]
        d_out[key] = [beta,'-p-']
    f.close()
    if os.path.exists(female_f) == False:
        return 'OK'
    f = open(female_f)
    f.readline()
    for line in f:
        if 'Inf' in line:
            continue
        line = line.strip().split()
        key = line[0]+'\t'+line[1]
        beta = line[2]+'p'+line[-1]
        if key in d_out:
            d_out[key][1] = beta
        else:
            d_out[key] = ['-p-',beta]
    f.close()
    i = 0
    Outdir = '/data2/myang9/SexAnnoDB/output/Network/sQTL/'
    outf = Outdir + '{0}_sQTL_summary.txt'.format(cancer_type)
    r = open(outf,'w')
    head = ['GeneID','GeneSymbol','EX_id','EX_info','Skipped_Exon','SNP_ID','SNPinfo','Male_effect','Male_FDR','Female_effect','Female_FDR','Sex_biased_ES_events']
    r.write('\t'.join(head)+'\n')
    for key in d_out:
        if 'rs' not in key:
            continue
        snpinfo = key.split('\t')[0]
        exid = key.split('\t')[1]
        if exid in diffES:
            diff_type = 'Yes'
        else:
            diff_type = '-'
        if exid not in d_exidchange:
            continue
        info = d_exidchange[exid]
        pos = d_pos[snpinfo]
        wl = info+[snpinfo.split('|')[1].split(':')[0],pos+':'+':'.join(snpinfo.split(':')[2:])]+d_out[key][0].split('p')+d_out[key][1].split('p') #+[diff_type]

        if '-' not in wl:
            if float(wl[8])<0.05 and float(wl[10])<0.05:
                if float(wl[7])>0 and float(wl[9])<0:
                    r.write('\t'.join(wl+[diff_type])+'\n')
                elif float(wl[7])<0 and float(wl[9])>0:
                    r.write('\t'.join(wl+[diff_type])+'\n')
                else:
                    pass
            elif float(wl[8])<0.05  and float(wl[10])>0.1:
                r.write('\t'.join(wl+[diff_type])+'\n')
            elif float(wl[8])>0.1 and float(wl[10])<0.05:
                r.write('\t'.join(wl+[diff_type])+'\n')

def sQTL_Summary_2(cancer_types):
    df_l = []
    for cancer_type in cancer_types:
        if cancer_type == 'LAML':
            continue
        sumf = '/data2/myang9/SexAnnoDB/output/Network/sQTL/{0}_sQTL_summary.txt'.format(cancer_type)
        df = pd.read_csv(sumf,sep = '\t')
        df['Cancer_Type'] = cancer_type
        df_l.append(df)
    result = pd.concat(df_l)
    result.index.name='Index'
    outf = '/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary.txt'
    result.to_csv(outf,sep='\t')

def sQTL_Summary_3(cancer_types):
    import ASdata_anlaysis as AA
    import Methyl_analysis as MA
    d_ORF = AA.ORF_anno()
    d_strand = MA.gene_strand()
    inf = '/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary.txt'
    f = open(inf)
    r = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary2.txt','w')
    r3 = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary3.txt','w')
    head = f.readline().strip().split('\t')
    head = head[0:6]+['ORF_anno']+head[6:8]+['SNP_position_to_ES']+head[8:]
    r.write('\t'.join(head)+'\n')
    r3.write('\t'.join(head)+'\n')
    i = 0
    for line in f:
        line = line.strip().split()
        if line[3] not in d_ORF:
            continue
        geneid = line[1]
        rspos = line[7].split(':')
        exinfo = line[4].split(':')
        strand = d_strand[geneid]
        chrom = exinfo[0]
        l = exinfo[1:]
        l.sort()
        l = [int(x) for x in l]
        s = l[2]
        e = l[3]
        pos_l = []
        rspos = int(rspos[1])
        if strand == '+':
            if rspos<=l[0]:
                pos_l.append('Distant upstream')
            elif  rspos>=l[5]:
                pos_l.append('Distant downstream')
            elif rspos>=s-2 and rspos<=s+2:
                pos_l.append('Acceptor splice site')
            elif rspos>=e-2 and rspos<=e+2:
                pos_l.append('Donor splice site')
            elif rspos > l[0] and rspos <s-2:
                pos_l.append('Proximal upstream')
            elif rspos > e+2 and rspos <l[5] :
                pos_l.append('Proximal downstream')
            elif rspos >s+2 and rspos <e-2:
                pos_l.append('Skipped exon')
            else:
                print('Error',rspos,l,strand)
        else:
            if rspos<=l[0]:
                pos_l.append('Distant downstream')
            elif  rspos>=l[5]:
                pos_l.append('Distant upstream')
            elif rspos>=s-2 and rspos<=s+2:
                pos_l.append('Donor splice site')
            elif rspos>=e-2 and rspos<=e+2:
                pos_l.append('Acceptor splice site')
            elif rspos > l[0] and rspos <s-2:
                pos_l.append('Proximal downstream')
            elif rspos > e+2 and rspos <l[5] :
                pos_l.append('Proximal upstream')
            elif rspos >s+2 and rspos <e-2:
                pos_l.append('Skipped exon')
            else:
                print('Error',rspos,l,strand)
        i = i+1
        wl = [str(i)] + line[1:6] + [d_ORF[line[3]]]+line[6:8] + pos_l+ line[8:]
        r.write('\t'.join(wl)+'\n')
        if  pos_l[0] in ['Donor splice site','Acceptor splice site','Skipped exon','Proximal upstream','Distant upstream']:
            if float(wl[13])<0.0001 and float(wl[11])>0.05:
                if abs(float(wl[10]))>0.01:
                    continue
                r3.write('\t'.join(wl)+'\n')
            elif float(wl[11])<0.0001 and float(wl[13])>0.05:
                if abs(float(wl[12]))>0.01:
                    continue
                r3.write('\t'.join(wl)+'\n')
            elif float(wl[13])<0.0001 and float(wl[11])<0.0001:
                if float(wl[10])>0 and float(wl[12])<0:
                    r3.write('\t'.join(wl)+'\n')
                elif float(wl[12])>0 and float(wl[10])<0:
                    r3.write('\t'.join(wl)+'\n')
            
        

def sQTL_boxplot(cancer_types):
    sQTL_boxplot_makefile(cancer_types)
    sQTL_boxplot_Draw()


def sQTL_boxplot_Draw():
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary3.txt')
    f.readline()
    d_score = {}
    for line in f:
        line = line.strip().split('\t')
        rsid = line[7]
        exid = line[3]
        cancer_type = line[-1]
        key = rsid + '_' + exid+'_'+cancer_type
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
            elif fdr_m < 0.05:
                sign_m = '"*"'
            else:
                sign_m = '"ns"'
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
            elif fdr_f<0.05:
                sign_f = '"*"'
            else:
                sign_f = 'ns'
        l = [sign_f,sign_m]
        #print(sign_f,sign_m,l)
        d_score[key] = l
    r = open('/data2/myang9/SexAnnoDB/output/sQTL_Draw.sh','w')
    print('Please run "bash /data2/myang9/SexAnnoDB/output/sQTL_Draw.sh"')
    i = 0
    path = '/data2/myang9/SexAnnoDB/output/Figures/sQTL_boxplot/file/'
    for inf in os.listdir(path):
        cancer_type = inf.split('-')[-2]
        key = inf.split('-')[1]+'_'+cancer_type
        sign_f = d_score[key][0]
        sign_m = d_score[key][1]
        inf = path + inf
        outf = inf.replace('file','figure').replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/boxplot_sQTL.r {0} {1} {2} {3} {4}'.format(inf,outf,sign_f,sign_m,cancer_type)
        i = i + 1
        '''
        if  'rs2911150_exon_skip_146352-KIRP' in inf: 
            print(d_score[key])
            print(sign_f,sign_m)
            print(cmd)
        '''
        r.write('nohup '+cmd+'  &\n')
        if i%10 ==0:
            r.write('wait;\n')

def sQTL_boxplot_makefile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/sQTL/SexAnnoDB_sQTL_Summary3.txt')
    d_sQTL = {}
    f.readline()
    for line in f:
        line = line.strip().split()
        
        cancer_type = line[-1]
        info = line[1:]
        if cancer_type in d_sQTL:
            d_sQTL[cancer_type].append(info)
        else:
            d_sQTL[cancer_type] = [info]
    path = '/data2/myang9/SexAnnoDB/data/sQTL_input/'
    sextype = ['female','male']
    for cancer_type in d_sQTL:
        snpf = path+'{0}_female_tumor_SNP.txt'.format(cancer_type)
        psif = path+'{0}_female_tumor_PSI.txt'.format(cancer_type)
        df_snpf = pd.read_csv(snpf,sep = '\t',index_col=0)
        df_psif = pd.read_csv(psif,sep = '\t',index_col=0)
        snpm = path+'{0}_male_tumor_SNP.txt'.format(cancer_type)
        psim = path+'{0}_male_tumor_PSI.txt'.format(cancer_type)
        df_snpm = pd.read_csv(snpm,sep = '\t',index_col=0)
        df_psim = pd.read_csv(psim,sep = '\t',index_col=0)
        df_snpf = df_snpf.T
        colnames = [x.split('|')[1].split(':')[0] for x in df_snpf.columns.values.tolist()]
        #print(colnames)
        df_snpf.columns = colnames
        df_psif = df_psif.T
        df_snpm = df_snpm.T
        colnames = [x.split('|')[1].split(':')[0] for x in df_snpm.columns.values.tolist()]
        df_snpm.columns = colnames
        df_psim = df_psim.T
        for line in d_sQTL[cancer_type]:
            ref = line[7].split(':')[2]
            alt = line[7].split(':')[3]
            rsid = line[6]
            exid = line[2]
            snpf_v = df_snpf[rsid] # .replace({0:ref+ref,1:ref+alt,2:alt+alt},inplace=True)
            psif_v = df_psif[exid]
            f_tmp = pd.concat([snpf_v,psif_v],axis=1)
            f_tmp['Group'] = 'Female'
            snpm_v = df_snpm[rsid] # .replace({0:ref+ref,1:ref+alt,2:alt+alt},inplace=True)
            psim_v = df_psim[exid]
            m_tmp = pd.concat([snpm_v,psim_v],axis=1)
            m_tmp['Group'] = 'Male'
            d_out = pd.concat([f_tmp,m_tmp])
            d_out.columns = ['x','y','Group']
            d_out.replace({'x':{0:ref+ref,1:ref+alt,2:alt+alt}},inplace=True)
            outf = '/data2/myang9/SexAnnoDB/output/Figures/sQTL_boxplot/file/{0}-{1}_{2}-{3}-sQTL.txt'.format(line[0],rsid,exid,cancer_type)
            d_out.index.name = 'Pids'
            d_out.dropna(axis=0,how='any',inplace=True)
            d_out.to_csv(outf,sep = '\t')

def eQTL_boxplot(cancer_types):
    eQTL_boxplot_makefile(cancer_types)
    eQTL_boxplot_Draw()

def eQTL_boxplot_makefile(cancer_types):
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTL/SexAnnoDB_eQTL_summary2.txt')
    d_eQTL = {}
    f.readline()
    i = 0
    for line in f:
        line = line.strip().split()
        cancer_type = line[-1]
        info = line[1:]
        if cancer_type in d_eQTL:
            d_eQTL[cancer_type].append(info)
        else:
            d_eQTL[cancer_type] = [info]
    path = '/data2/myang9/SexAnnoDB/data/eQTL_input/'
    sextype = ['female','male']
    for cancer_type in d_eQTL:
            snpf = path+'{0}_female_tumor_SNP.txt'.format(cancer_type)
            psif = path+'{0}_female_tumor_FPKM.txt'.format(cancer_type)
            df_snpf = pd.read_csv(snpf,sep = '\t',index_col=0)
            df_psif = pd.read_csv(psif,sep = '\t',index_col=0)
            snpm = path+'{0}_male_tumor_SNP.txt'.format(cancer_type)
            psim = path+'{0}_male_tumor_FPKM.txt'.format(cancer_type)
            df_snpm = pd.read_csv(snpm,sep = '\t',index_col=0)
            df_psim = pd.read_csv(psim,sep = '\t',index_col=0)
            df_snpf = df_snpf.T
            colnames = [x.split('|')[1].split(':')[0] for x in df_snpf.columns.values.tolist()]
            df_snpf.columns = colnames
            df_psif = df_psif.T
            df_snpm = df_snpm.T
            colnames = [x.split('|')[1].split(':')[0] for x in df_snpm.columns.values.tolist()]
            df_snpm.columns = colnames
            df_psim = df_psim.T
            
            for line in d_eQTL[cancer_type]:
                ref = line[3].split(':')[2]
                alt = line[3].split(':')[3]
                rsid = line[2]
                geneid = line[0]
                snpf_v = df_snpf[rsid]
                snpf_v.replace({0:ref+ref,1:ref+alt,2:alt+alt},inplace=True)
                psif_v = df_psif[geneid]
                f_tmp = pd.concat([snpf_v,psif_v],axis=1)
                f_tmp['Group'] = 'Female'
                snpm_v = df_snpm[rsid]#.replace({0:ref+ref,1:ref+alt,2:alt+alt},inplace=True)
                snpm_v.replace({0:ref+ref,1:ref+alt,2:alt+alt},inplace=True)
                psim_v = df_psim[geneid]
                m_tmp = pd.concat([snpm_v,psim_v],axis=1)
                m_tmp['Group'] = 'Male'
                d_out = pd.concat([f_tmp,m_tmp])
                #d_out.replace({0.0:ref+ref,1.0:ref+alt,2.0:alt+alt},inplace=True)
                d_out.columns = ['x','y','Group']
                outf = '/data2/myang9/SexAnnoDB/output/Figures/eQTL_boxplot/file/{0}-{1}_{2}-{3}-sQTL.txt'.format(geneid,rsid,geneid,cancer_type)
                d_out.index.name = 'Pids'
                i = i +1 
                #print(i,outf)
                d_out.dropna(axis=0,how='any',inplace=True)
                d_out.to_csv(outf,sep = '\t')

def eQTL_boxplot_Draw():
    f = open('/data2/myang9/SexAnnoDB/output/Network/eQTL/SexAnnoDB_eQTL_summary2.txt')
    f.readline()
    d_score = {}
    for line in f:
        line = line.strip().split('\t')
        rsid = line[3]
        geneid = line[1]
        cancer_type = line[-1]
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
                sign_f = 'ns'
        l = [sign_f,sign_m]
        d_score[key] = l
    r = open('eQTL_Draw.sh','w')
    i = 0
    print('Plsease run '+'nohup bash eQTL_Draw.sh  &\n')
    path = '/data2/myang9/SexAnnoDB/output/Figures/eQTL_boxplot/file/'
    for inf in os.listdir(path):
        cancer_type = inf.split('-')[-2]
        key = inf.split('-')[1]+'_'+cancer_type
        sign_f = d_score[key][0]
        sign_m = d_score[key][1]
        inf = path + inf
        outf = inf.replace('file','figure').replace('.txt','.svg')
        cmd = 'Rscript /data2/myang9/SexAnnoDB/script/Draw/boxplot_eQTL.r {0} {1} {2} {3} {4}'.format(inf,outf,sign_f,sign_m,cancer_type)
        i = i + 1
        #print('Plsease run '+'nohup '+cmd+'  &\n')
        r.write('nohup '+cmd+'  &\n')
        if i%10 ==0:
            r.write('wait;\n')















