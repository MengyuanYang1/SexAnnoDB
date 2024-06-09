#system("conda activate r44")


library(TCGAbiolinks)
library(maftools)


cancers <- c('BRCA','KIRC','LUAD','LGG','THCA','HNSC','LUSC','SKCM','COAD','BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL')


inf = '/data2/myang9/SexAnnoDB/output/MutationInfo/SexAnnoDB_Sexbiased_mutation.txt'
genes =  read.table(inf,header=TRUE,sep='\t')

download <- function() {
for (cancer in cancers){
	#genelist <- subset(genes,Cancer = cancer)$genesym
	genelist <- genes[which(genes$can == cancer), ]$genesym
	project = paste0('TCGA-',cancer)
	#query <- GDCquery(project = project, data.category = "Simple Nucleotide Variation",data.type = "Masked Somatic Mutation",access = "open")
	query <- GDCquery(project = project, data.category = "Copy Number Variation",data.type = "Masked Copy Number Segment",access = "open")
	GDCdownload(query)
	#outR = paste0(project,"_SNP.Rdata")
	outR = paste0(project,"_CNV.Rdata")
	GDCprepare(query, save = T,save.filename = outR)
	#load(outR)
}}
#cancers <- c('BLCA','STAD','LIHC','KIRP','SARC','PCPG','PAAD','READ','GBM','ESCA','LAML','THYM','MESO','UVM','ACC','KICH','DLBC','CHOL')

#analysis <- {
for (cancer in cancers){
	genelist <- genes[which(genes$CancerType == cancer), ]$Hugo_Symbol
	project = paste0('TCGA-',cancer)
	outf = paste0("/data2/myang9/SexAnnoDB/data/TCGA/mutation/",project,"_SNV.modify.filter.maf")
	#write.table(data,file=outf,quote=F,col.name=T,row.names=F)
	data<-read.maf(maf=outf, clinicalData="/data2/myang9/SexAnnoDB/data/TCGA/clinical.tsv")
	Variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",   "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation",   "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
	
	clin<- data@clinical.data
	clin.female <- subset(clin, gender=="female")$Tumor_Sample_Barcode
	clin.male <- subset(clin, gender=="male")$Tumor_Sample_Barcode
	#cat(project,as.character(length(clin.female)),as.character(length(clin.male)),sep=",")
	
	data.female <- subsetMaf(maf=data, tsb=clin.female, isTCGA=TRUE)
	data.male <- subsetMaf(maf=data, tsb=clin.male, isTCGA=TRUE)
	
	#fvsm <- mafCompare(m1=data.male, m2=data.female, m1Name="Male", m2Name="Female", minMut=5)
	#result <- paste0("/data2/myang9/SexAnnoDB/output/MutationInfo/tmp",project,"_maleVSfemale.txt")
	#write.table(fvsm$results, file=result, quote=FALSE, row.names=FALSE, sep="\t")
	
	#forestPlot(mafCompareRes=fvsm, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
	#outf = paste0("/data2/myang9/SexAnnoDB/output/MutationInfo/oncoplot",project,"_oncoplot.svg")
	#p <- coOncoplot(m1=data.female, m2=data.male, m1Name="Female", m2Name="male", genes=genelist)
	#cat(length(genelist))
	for (genesym in genelist){
		cat(project,cancer,genesym)
		outf =paste0('/data2/myang9/SexAnnoDB/output/MutationInfo/lollipopPlot2/',project,'_',genesym,"_lollipopPlot.svg")
		if (genesym %in% c('TTC6','C6orf15','CRACD','BCL9','FSCB','KIAA1217','RALGAPB','KIAA1210','FAM169A','SPATA31A1','RTL9','TENT5D','CFAP221','NOL4','CATSPERD','ERMN','MUC5AC','MRE11','SPATA31D1','DRC7','EXPH5','KATNIP','UNC79','TMEM132B','TMEM132D','RELCH','FMR1NB','PRB3','MAMLD1','FAM120C','NEXMIF','ZNF761','C22orf42','C2CD6','NMS','LSR','RAI2','OBI1','CCDC181','KPRP','C14orf39','FBXO40','UNC80','MARF1','SZT2')){next;}
		#temp <- try(lollipopPlot2(m1=data.female, m2=data.male, m1_name="Female", m2_name="Male", gene=genesym, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short"),TRUE)
		if ("try-error" %in% class(temp)) {temp<-NA} 
		{
		svg(outf,10,4)
		lollipopPlot2(m1=data.female, m2=data.male, m1_name="Female", m2_name="Male", gene=genesym, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short",showDomainLabel=FALSE)
		dev.off()
		}
		
	}}
	#svg(outf,8,2+length(genelist)*0.2)
	#coOncoplot(m1=data.female, m2=data.male, m1Name="Female", m2Name="male", removeNonMutated=FALSE,gene_mar=2,titleFontSize=1,genes=genelist,legendFontSize=1)
	#dev.off()


#}
