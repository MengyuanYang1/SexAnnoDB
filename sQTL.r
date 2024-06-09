library(MethylMix)
library(MatrixEQTL)



input_dir  = '/data2/myang9/SexAnnoDB/data/sQTL_input/'
output_dir = '/data2/myang9/SexAnnoDB/output/Network/sQTL/'
useModel = modelLINEAR
args = commandArgs(T)
cancertype = args[1]
type = args[2]

pvOutputThreshold_tra = 0.0000001
pvOutputThreshold_cis = 1
pvOutputThreshold = 0.00001
errorCovariance = numeric();

snp_file_name=paste0(input_dir,cancertype,'_',type,'_SNP.txt')
PSI_file_name= paste0(input_dir,cancertype,'_',type,'_PSI_select.txt')
output_file_name_all=paste0(output_dir,cancertype,'_',type,"_all_sQTL_out.txt")
output_file_name_cis = paste0(output_dir,cancertype,'_',type,"_sQTL_cis",".txt")
output_file_name_tra = paste0(output_dir,cancertype,'_',type,"_sQTL_tra",".txt")

## Load SNP data
print(snp_file_name)
beta_value = SlicedData$new();
beta_value$fileDelimiter = "\t";      # the TAB character
beta_value$fileOmitCharacters = "NA"; # denote missing values;
beta_value$fileSkipRows = 1;          # one row of column labels
beta_value$fileSkipColumns = 1;       # one column of row labels
beta_value$fileSliceSize = 2000;      # read file in slices of 2,000 rows
beta_value$LoadFile(snp_file_name);

print(PSI_file_name)
## Load PSI data
AS_EVENT = SlicedData$new();
AS_EVENT$fileDelimiter = "\t";      # the TAB character
AS_EVENT$fileOmitCharacters = "NA"; # denote missing values;
AS_EVENT$fileSkipRows = 1;          # one row of column labels
AS_EVENT$fileSkipColumns = 1;       # one column of row labels
AS_EVENT$fileSliceSize = 2000;      # read file in slices of 2,000 rows
AS_EVENT$LoadFile(PSI_file_name);

input_dir = '/data2/myang9/SexAnnoDB/data/QTL_position/'
snp_pos_file_name = paste(input_dir,cancertype,'_female_tumor_SNP_pos_hg38.txt',sep="")
psi_pos_file_name = paste(input_dir,"EX_position.txt",sep="")

mepos = read.table(snp_pos_file_name, header = TRUE, stringsAsFactors = FALSE);
psicoord = read.table(psi_pos_file_name, header = TRUE, stringsAsFactors = FALSE);
cisDist=10e6

me = Matrix_eQTL_main(
snps = beta_value,
gene = AS_EVENT,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
snpspos = mepos,
genepos = psicoord,
cisDist = cisDist,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);


ALL <-function(){
me = Matrix_eQTL_engine(
    snps = beta_value,
    gene = AS_EVENT,
    output_file_name = output_file_name_all ,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
}




