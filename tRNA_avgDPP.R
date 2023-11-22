PPname <- as.character(commandArgs(trailingOnly=TRUE)[1]);#
chrxx <- as.character(commandArgs(trailingOnly=TRUE)[2]);#chrxx<-"chrX"
path.tss <- as.character(commandArgs(trailingOnly=TRUE)[3]);# path.tss <- "/home/uemura/tRNA/tss/GSE132660_SupplementaryData1.txt"
path.Physical_properties <- as.character(commandArgs(trailingOnly=TRUE)[4]);# path.Physical_properties <- "/home/uemura/DPP_hg38/Duplex_free_energy/chrX.txt"
path.out <- as.character(commandArgs(trailingOnly=TRUE)[5]);#

print(c(PPname, chrxx, path.tss, path.Physical_properties, path.out))

dir.create(path.out)
setwd(path.out);
if(!sum(list.files(path.out) == PPname)){
  dir.create(PPname);
}
setwd(PPname);

distance_fromTSS.U <- 50
distance_fromTSS.D <- 50

antisense_adjust_num <- 1
if(length(grep("Bruckner", path.Physical_properties))){
  antisense_adjust_num <- 2
}else if(length(grep("Packer", path.Physical_properties))){
  antisense_adjust_num <- 3
}

tss_data <- read.table(path.tss , stringsAsFactors=F);
Physical_properties <- scan(path.Physical_properties, what=1)

tss_data_ex <- tss_data[,c(1,2,3,4,10,11)]
colnames(tss_data_ex) <- c("chr", "start", "strand", "TPM","gene_type", "gene_id");

# tRNAの遺伝子を抽出する
tss_data_ex <- tss_data_ex[which(tss_data_ex[,"chr"]==chrxx),]
tss_data_ex <- tss_data_ex[which(tss_data_ex[,"gene_type"]=="tRNAscan"),]

dat_TSS_sense <- tss_data_ex[which(tss_data_ex[,"strand"]=="+"),"start"]
dat_TSS_anti <- tss_data_ex[which(tss_data_ex[,"strand"]=="-"),"start"]

dat_gene_id_sense <- tss_data_ex[which(tss_data_ex[,"strand"]=="+"),"gene_id"]
dat_gene_id_anti <- tss_data_ex[which(tss_data_ex[,"strand"]=="-"),"gene_id"]

# tpm expression
dat_gene_tpm_sense <- tss_data_ex[which(tss_data_ex[,"strand"]=="+"),"TPM"]
dat_gene_tpm_anti <- tss_data_ex[which(tss_data_ex[,"strand"]=="-"),"TPM"]


if(!(length(dat_TSS_sense) == length(dat_gene_id_sense))){
  print("the number of TSS and id incorrect")
  q("no")
}

if(!(length(dat_TSS_anti) == length(dat_gene_id_anti))){
  print("the number of TSS and id incorrect")
  q("no")
}


PP_all <- c(); gene_id.all <- c(); gene_tpm.all <- c();
for(i_tss_sense in 1:length(dat_TSS_sense)){

  if(!length(dat_TSS_sense)){break}

  i_tss_Pos <- as.numeric(dat_TSS_sense[i_tss_sense])
  if(((i_tss_Pos-distance_fromTSS.U)<0)|((i_tss_Pos+distance_fromTSS.D)>length(Physical_properties))){next}
  i_PP <- Physical_properties[(i_tss_Pos-distance_fromTSS.U):(i_tss_Pos+distance_fromTSS.D)];

  if(sum(i_PP == 55555, na.rm=T)){next}
  if(sum(is.na(i_PP))){next}

  # gene id
  gene_id.all <- c(gene_id.all, dat_gene_id_sense[i_tss_sense])

  # gene tpm
  gene_tpm.all <- c(gene_tpm.all, dat_gene_tpm_sense[i_tss_sense])

  PP_all <- rbind(PP_all, i_PP)
}
for(i_tss_anti in 1:length(dat_TSS_anti)){

  if(!length(dat_TSS_anti)){break}

  i_tss_Pos <- as.numeric(dat_TSS_anti[i_tss_anti])
  if(((i_tss_Pos-distance_fromTSS.D)<0)|((i_tss_Pos+distance_fromTSS.U)>length(Physical_properties))){next}
  i_PP <- Physical_properties[(i_tss_Pos+distance_fromTSS.U-antisense_adjust_num):(i_tss_Pos-distance_fromTSS.D-antisense_adjust_num)];

  if(sum(i_PP == 55555, na.rm=T)){next}
  if(sum(is.na(i_PP))){next}

  # gene if
  gene_id.all <- c(gene_id.all, dat_gene_id_anti[i_tss_anti])

  # gene tpm
  gene_tpm.all <- c(gene_tpm.all, dat_gene_tpm_anti[i_tss_anti])

  PP_all <- rbind(PP_all, i_PP)
}

rownames(PP_all) <- gene_id.all
PP_all <- cbind(PP_all, gene_tpm.all)
colnames(PP_all) <- c((-distance_fromTSS.U):distance_fromTSS.D, "tpm")

#pdf(sprintf("%s_%s.pdf",PPname,chrxx))
#plot(apply(PP_all, 2, mean), type="l")
#dev.off()

#PP_all_norm <- (PP_all-min(PP_all))*(1-(-1))/(max(PP_all)-min(PP_all))+(-1);
#plot(apply(PP_all_norm, 2, mean), typel="l")

write.table(PP_all, col.names=T, row.names=T, file=sprintf("%s_%s.txt",PPname,chrxx))
