setwd("/Users/uemurakohei/Desktop/ncRNA/tRNA/sequence")

seq <- c()

for(i_chr in c(1:22, "X", "Y")){

  if(!length(scan(sprintf("chr%s.txt", i_chr),""))){
    next
  }

  if(sum(list.files() == sprintf("chr%s.txt", i_chr))){


  seq <- rbind(seq, read.table(sprintf("chr%s.txt", i_chr), skip=1, stringsAsFactors=F))

  }

}

seq_reptss <-c()
for(i_tss in unique(seq[,1])){
  each_gene <- seq[which(seq[,1]==i_tss),]
  seq_reptss <- c(seq_reptss, paste(toupper(each_gene[which.max(each_gene[,ncol(each_gene)]),2:102]), collapse=""))
}

write.table(seq_reptss, file="seq_reptss_for_weblogo.txt", col.names=F, row.names=F, quote=F)
