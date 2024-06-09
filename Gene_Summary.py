import os,sys

def Gene_Summary():
    f = open('/data2/myang9/library/human/gencode/hg38/gencode.v22.annotation.gtf')
    f.readline()
    d_info = {}
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 5:
            continue
        if line[2]!='gene':
            continue
        geneid = line[8].split('"')[1].split('.')[0]
        genetype = line[8].split('"')[3]
        genename = line[8].split('"')[7]
        #print(geneid,genetype,genename)
        d_info[geneid] = [genetype,genename]
    f.close()
    
    f = open('/data2/myang9/SexAnnoDB/output/Genetable/hgnc_complete_set.txt')
    head = f.readline()
    r = open('/data2/myang9/SexAnnoDB/output/Genetable/SexDiff_GeneSummary_PC.txt','w')
    head = ['Index','Gene ID','Symbol','HGNC','Entrez ID','Gene Type','Gene Name','Cytomap','Alias Symbol','Uniprot']
    r.write('\t'.join(head)+'\n')
    i = 0
    for line in f:
        line = line.split('\t')    
        i = i +1
        geneid = line[19]
        if geneid not in d_info:
            continue
        HGNC = line[0].split(':')[1]
        full_name = line[2]
        Cytomap = line[6]
        Alias = [line[1]]+line[8].replace('"','').split('|')
        Alias = '|'.join(Alias)
        print(Alias)
        l =[str(i),geneid,d_info[geneid][1],HGNC,line[18],d_info[geneid][0],full_name,Cytomap,Alias,line[25]]
        if l[5] != 'protein_coding':
            continue
        r.write('\t'.join(l)+'\n')

    


