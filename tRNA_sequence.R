chrxx <- as.character(commandArgs(trailingOnly=TRUE)[1]);#chrxx<-"chrX"
path.tss <- as.character(commandArgs(trailingOnly=TRUE)[2]);# path.tss <- "/home/uemura/tRNA/tss/GSE132660_SupplementaryData1.txt"
path.sequence <- as.character(commandArgs(trailingOnly=TRUE)[3]);# path.sequence <- "/home/uemura/hg38/chrX.fa"
path.out <- as.character(commandArgs(trailingOnly=TRUE)[4]);#


dir.create(path.out)
setwd(path.out);

distance_fromTSS.U <- 50
distance_fromTSS.D <- 50


tss_data <- read.table(path.tss , stringsAsFactors=F);
sequence <- scan(path.sequence, what="", skip=1)

sequence <- unlist(strsplit(sequence, ""))

tss_data_ex <- tss_data[,c(1,2,3,4,10,11)]
colnames(tss_data_ex) <- c("chr", "start", "strand", "TPM", "gene_type", "gene_id");

# tRNAの遺伝子を抽出する
tss_data_ex <- tss_data_ex[which(tss_data_ex[,"chr"]==chrxx),]
tss_data_ex <- tss_data_ex[which(tss_data_ex[,"gene_type"]=="tRNAscan"),]

dat_TSS_sense <- tss_data_ex[which(tss_data_ex[,"strand"]=="+"),"start"]
dat_TSS_anti <- tss_data_ex[which(tss_data_ex[,"strand"]=="-"),"start"]

dat_gene_id_sense <- tss_data_ex[which(tss_data_ex[,"strand"]=="+"),"gene_id"]
dat_gene_id_anti <- tss_data_ex[which(tss_data_ex[,"strand"]=="-"),"gene_id"]

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


sequence_all <- c(); gene_id.all <- c(); gene_tpm.all <- c()
for(i_tss_sense in 1:length(dat_TSS_sense)){

if(!length(dat_TSS_sense)){break}

i_tss_Pos <- as.numeric(dat_TSS_sense[i_tss_sense])
  if(((i_tss_Pos-distance_fromTSS.U)<0)|((i_tss_Pos+distance_fromTSS.D)>length(sequence))){next}
  i_sequence <- sequence[(i_tss_Pos-distance_fromTSS.U):(i_tss_Pos+distance_fromTSS.D)];

  #if(sum(i_sequence == 55555, na.rm=T)){next}
  #if(sum(is.na(i_sequence))){next}

  gene_id.all <- c(gene_id.all, dat_gene_id_sense[i_tss_sense])

  gene_tpm.all <- c(gene_tpm.all, dat_gene_tpm_sense[i_tss_sense])

  sequence_all <- rbind(sequence_all, i_sequence)
}
for(i_tss_anti in 1:length(dat_TSS_anti)){

if(!length(dat_TSS_anti)){break}

 i_tss_Pos <- as.numeric(dat_TSS_anti[i_tss_anti])
  if(((i_tss_Pos-distance_fromTSS.D)<0)|((i_tss_Pos+distance_fromTSS.U)>length(sequence))){next}
  i_sequence <- sequence[(i_tss_Pos+distance_fromTSS.U):(i_tss_Pos-distance_fromTSS.D)];

  for(j_i_sequence in 1:length(i_sequence)){

      if(i_sequence[j_i_sequence] == "A"){
        i_sequence[j_i_sequence] <- "T"
      }else if(i_sequence[j_i_sequence] == "T"){
        i_sequence[j_i_sequence] <- "A"
      }else if(i_sequence[j_i_sequence] == "G"){
        i_sequence[j_i_sequence] <- "C"
      }else if(i_sequence[j_i_sequence] == "C"){
        i_sequence[j_i_sequence] <- "G"
      }
      
    }

  #if(sum(i_sequence == 55555, na.rm=T)){next}
  #if(sum(is.na(i_sequence))){next}

  gene_id.all <- c(gene_id.all, dat_gene_id_anti[i_tss_anti])

  gene_tpm.all <- c(gene_tpm.all, dat_gene_tpm_anti[i_tss_anti])

  sequence_all <- rbind(sequence_all, i_sequence)
}

rownames(sequence_all) <- gene_id.all
sequence_all <- cbind(sequence_all, gene_tpm.all)

colnames(sequence_all) <- c((-distance_fromTSS.U):(distance_fromTSS.D), "tpm")

#sequence_all_norm <- (sequence_all-min(sequence_all))*(1-(-1))/(max(sequence_all)-min(sequence_all))+(-1);
#plot(apply(sequence_all_norm, 2, mean), typel="l")

write.table(sequence_all, col.names=T, row.names=T, file=sprintf("%s.txt",chrxx))
