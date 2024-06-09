cancer_types_27 <- c('BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL')

path = '/data2/myang9/SexAnnoDB/data/TCGA/CNV/'
PreprocessingCNV <- function(cancer){
	cnvf <- paste0(path,'TCGA-',cancer,'_CNV.Rdata')
	print(cnvf)
	cnv <- load(cnvf)
	hnsc_seg <- eval(parse(text=cnv))
	hnsc_seg <- hnsc_seg[, -1]
	hnsc_seg <- hnsc_seg[,c('Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean' )]
	tumor_seg <- hnsc_seg[substr(hnsc_seg$Sample,14,15) == "01" ,]
	segf <-  paste0(path,'TCGA-',cancer,'_seg_all.txt')
	write.table(hnsc_seg, segf, sep = "\t", col.names = TRUE, row.names=F)
	snpf <- paste0(path,"snp6.na35.remap.hg38.subset.txt")
	hg_marker_file <- read.delim(snpf)
	hg_marker_file <- hg_marker_file[hg_marker_file$freqcnv =="FALSE" ,]
	hg_marker_file <- hg_marker_file[,c(1,2,3)]
	markerf <- paste0(path,cancer,'_hg_marker_file.txt')
	write.table(hg_marker_file, markerf,sep = "\t", col.names = TRUE, row.names=F)
}


for (cancer in cancer_types_27) {
	PreprocessingCNV(cancer)	
}

clinical_all <- read.delim("clinical.tsv",check.names = F)
clinical_all <- clinical_all[which(clinical_all$gender %in% c('male','female')), ]
for (cancer in cancers){
  clinical <- clinical_all[which(clinical_all$project_id == paste0('TCGA-',cancer)), ]
  print(table(clinical$gender))
}
for (cancer in cancers){
  #cancer = 'SARC'
  clinical <- clinical_all[which(clinical_all$project_id == paste0('TCGA-',cancer)), ]
  samples <- rownames(clinical)
  clinical <- clinical[,c('Tumor_Sample_Barcode','gender')]
  all_data = read.delim(paste0('CNV/',cancer,"/all_data_by_genes.txt"),row.names = 1,check.names = F)
  all_thr = read.delim(paste0('CNV/',cancer,"/all_thresholded.by_genes.txt"),row.names = 1,check.names = F)
  cnv = all_thr[,-(1:2)]
  cnv_samples <- colnames(cnv)
  cnv_tumor = c(); for (x in cnv_samples){  sub_code <- substr(x,14,15); if (as.numeric(sub_code) < 10){  cnv_tumor <- append(cnv_tumor,x) } }
  cnv <- cnv[,as.vector(cnv_tumor)]
  cnv_tumor_barcode <- c(); for (x in cnv_tumor){ sub_x <- substr(x,1,12); cnv_tumor_barcode <- append(cnv_tumor_barcode,sub_x)}
  colnames(cnv) <- cnv_tumor_barcode
  cnv <- cnv[,order(colnames(cnv))]
  clinical <- clinical[which(clinical$Tumor_Sample_Barcode %in% cnv_tumor_barcode), ]
  clinical <- clinical[ !duplicated(clinical),  ]
  clinical <- clinical[order(clinical$Tumor_Sample_Barcode),]
  print(cancer)
  cnv[cnv >= 1] = 1
  cnv[cnv <= -1] = -1
  rownames(cnv) = 1:nrow(cnv)
  diffcnv = differential_CNV(cnv,clinical$gender)
  diffcnv$gene = rownames(all_thr)
  diffcnv$Cytoband = all_thr$Cytoband
  diffcnv <- diffcnv[which(diffcnv$adj.P.Val <0.05), ]
  outf <- paste0(cancer,'_diffcnv.txt')
  print(outf)
  write.table(diffcnv, outf, sep = "\t", col.names = TRUE, row.names=T)
  
}
















