for(para in c("DNA_bending_stiffness", "DNA_denaturation", "Duplex_disrupt_energy", "Duplex_free_energy" ,"Protein_induced_deformability_Bp" ,"Stabilizing_energy_of_Z_DNA_AS", "Stabilizing_energy_of_Z_DNA_SA", "Stacking_energy", "Bruckner" ,"Packer")){

  # para <- "DNA_bending_stiffness"

  setwd(sprintf("/Users/uemurakohei/Desktop/ncRNA/tRNA/avgDPP_a_b_box/%s",para))

  tss_data  <- read.table("/Users/uemurakohei/Desktop/ncRNA/tRNA/tss/GSE132660_SupplementaryData1.txt" , stringsAsFactors=F)
  tss_data_tRNA <- tss_data[which(tss_data[,10] == "tRNAscan"),]

  dpp <- c()

  for(i_chr in c(1:22,"X","Y","M")){
    if(!sum(list.files()==sprintf("%s_chr%s_Abox.txt", para,i_chr))){next}
    dpp <- rbind(dpp, read.table(sprintf("/Users/uemurakohei/Desktop/ncRNA/tRNA/avgDPP_a_b_box/%s/%s_chr%s_Abox.txt", para,para,i_chr), header=T))
  }

  #
  dpp <- dpp[which(rownames(dpp) %in% unique(unlist(strsplit(tss_data_tRNA[,11], ";")))), ]


  #dpp_reptss <- c()
  #for(i_tss in unique(dpp[,1])){

  #  dpp.i_gene <- dpp[which(dpp[,1] == i_tss),]

  #  dpp_reptss <- rbind(dpp_reptss, dpp.i_gene[which.max(dpp.i_gene[,ncol(dpp)]),])

  #}

  colnames(dpp) <- -50:50

  x<- 0:30
  mean_force<- apply(dpp[,as.character(x)], 2, function(x) mean(as.numeric(x)))
  #assume constant standard deviation across the
  sd<-apply(dpp[,as.character(x)], 2, function(x) sd(as.numeric(x)))
  #determine error band
  psd<-mean_force+sd
  nsd<-mean_force-sd


  #y_lim
  source("/Users/uemurakohei/Desktop/ncRNA/tRNA/tRNA_ylim.R")


  pdf(sprintf("%s_avg_A_box.pdf", para),width=3.4, height=4)
  #dev.new(width=5, height=4)
  par(mar = c(4, 5, 2, 1))
  par(family = "Helvetica")
  plot(x, mean_force, ty="l", col="black",
       ylab="",
       xlab="",
       las=1,
       cex.axis = 1.8,
       #main=property_name,
       lty=1,lwd=2,

       #ylim=c(min(nsd), max(psd))
       ylim=c(y_lim)
     )
  #draw boundary and fill
  lines(x, psd, col="gray")
  lines(x, nsd, col="gray")

  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
  #redraw line on top

  #shuffle
  #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
  #lines(x, mean_shuffle, col="red", lwd=3)
  lines(x, mean_force, col="black",lwd=2)
  #rect(x[11], y_lim[1], x[11+11], y_lim[1]+1, col="red")
  lines(x[11]:x[11+11], rep(y_lim[1], length(x[11]:x[11+11])), col="red",lwd=2, lty="dashed")

  dev.off()

}
