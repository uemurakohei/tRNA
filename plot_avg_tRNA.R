for(para in c("DNA_bending_stiffness", "DNA_denaturation", "Duplex_disrupt_energy", "Duplex_free_energy" ,"Protein_induced_deformability_Bp" ,"Stabilizing_energy_of_Z_DNA_AS", "Stabilizing_energy_of_Z_DNA_SA", "Stacking_energy", "Bruckner" ,"Packer")){

  setwd(sprintf("/Users/uemurakohei/Desktop/ncRNA/tRNA/avgDPP/%s",para))

  dpp <- c()

  for(i_chr in c(1:22,"X","Y","M")){
    if(!sum(list.files()==sprintf("%s_chr%s.txt", para,i_chr))){next}
    dpp <- rbind(dpp, read.table(sprintf("/Users/uemurakohei/Desktop/ncRNA/tRNA/avgDPP/%s/%s_chr%s.txt", para,para,i_chr), skip=1))
  }

  dpp_reptss <- c()
  for(i_tss in unique(dpp[,1])){

    dpp.i_gene <- dpp[which(dpp[,1] == i_tss),]

    dpp_reptss <- rbind(dpp_reptss, dpp.i_gene[which.max(dpp.i_gene[,ncol(dpp)]),])

  }

  x<- -50:50
  mean_force<- apply(dpp_reptss[,2:102], 2, mean)
  #assume constant standard deviation across the
  sd<-apply(dpp_reptss[,2:102], 2, sd)
  #determine error band
  psd<-mean_force+sd
  nsd<-mean_force-sd


  #y_lim
  source("/Users/uemurakohei/Desktop/ncRNA/tRNA/tRNA_ylim.R")


  pdf(sprintf("%s_average.pdf", para),width=6, height=4)
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

       #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
       ylim=c(y_lim))
  #draw boundary and fill
  lines(x, psd, col="gray")
  lines(x, nsd, col="gray")

  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
  #redraw line on top

  #shuffle
  #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
  #lines(x, mean_shuffle, col="red", lwd=3)
  lines(x, mean_force, col="black",lwd=2)
  dev.off()

  # tss　下流5まで
  x<- -50:5
  mean_force<- apply(dpp_reptss[,2:57], 2, mean)
  #assume constant standard deviation across the
  sd<-apply(dpp_reptss[,2:57], 2, sd)
  #determine error band
  psd<-mean_force+sd
  nsd<-mean_force-sd


  source("/Users/uemurakohei/Desktop/ncRNA/tRNA/tRNA_ylim.R")



  pdf(sprintf("%s_average_tss5.pdf", para),width=6, height=4)
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

       #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
       ylim=c(y_lim))
  #draw boundary and fill
  lines(x, psd, col="gray")
  lines(x, nsd, col="gray")

  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
  #redraw line on top

  #shuffle
  #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
  #lines(x, mean_shuffle, col="red", lwd=3)

  #縦線
  abline(v=0, col="red", lty="dotted", lwd=4)
  abline(v=-27, col="blue", lty="dotted", lwd=4)

  lines(x, mean_force, col="black",lwd=2)
  #points(-27, mean_force[which(x == -27)], col=rgb(1,0,0),cex = 1.5, lwd=2)
  dev.off()

  # tss　下流10まで
  x<- -50:10
  mean_force<- apply(dpp_reptss[,2:62], 2, mean)
  #assume constant standard deviation across the
  sd<-apply(dpp_reptss[,2:62], 2, sd)
  #determine error band
  psd<-mean_force+sd
  nsd<-mean_force-sd



  pdf(sprintf("%s_average_tss10.pdf", para),width=6, height=4)
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

       #ylim=c(min(c(data_mRNA_averagesd[,"value"]-data_mRNA_averagesd[,"sd"])),max(c(data_mRNA_averagesd[,"value"]+data_mRNA_averagesd[,"sd"]))))
       ylim=c(y_lim))
  #draw boundary and fill
  lines(x, psd, col="gray")
  lines(x, nsd, col="gray")

  polygon(x=c(x, rev(x)), y=c(psd, rev(nsd)), col="gray", density = 100, angle=90)
  #redraw line on top

  #shuffle
  #lines(x, mean_shuffle, col=rgb(1,0,0, alpha=0.5), lwd=1,lty=1)
  #lines(x, mean_shuffle, col="red", lwd=3)
  lines(x, mean_force, col="black",lwd=2)
  points(-27, mean_force[which(x == -27)], col=rgb(1,0,0),lwd=2)
  dev.off()
}
