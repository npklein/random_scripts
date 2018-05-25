
QQplotFromPs<- function(VctPvalues,outputfile,name=''){
  VctPvalues <- VctPvalues[!is.na(VctPvalues)]
  ## Plot function ##
  plotQQ <- function(p,color){
    #p <- 2*pnorm(-abs(z))
    
    p <- sort(p)
    expected <- c(1:length(p))
    lobs <- -(log10(p))
    lexp <- -(log10(expected / (length(expected)+1)))
    
    # plots all points with p < 1e-3
    p_sig = subset(p,p<0.001)
    points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=.4, col=color, bg=color)
    
    # samples 2500 points from p > 1e-3
    n=2501
    i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
    lobs_bottom=subset(lobs[i],lobs[i] <= 3)
    lexp_bottom=lexp[i[1:length(lobs_bottom)]]
    points(lexp_bottom, lobs_bottom, pch=23, cex=.4, col=color, bg=color)
  }
  
  
  z=qnorm(VctPvalues/2)
#  print(z)
  maxP <- max(-log10(VctPvalues)) +1
#  print(maxP)
  ## calculates lambda
  lambda = round(median(z^2)/qchisq(0.5,df=1),5)
  print(name)
  print(paste('lambda:',lambda))
  ## Plots axes and null distribution
  pdf(outputfile, width=6, height=6)
  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,7), ylim=c(0,maxP), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))
  
  ## plots data
  plotQQ(VctPvalues,"black");
  
  ## provides legend
  legend(.25,7,legend=c("Expected (null)","Observed"), pch=c((vector("numeric",2)+1)*23), cex=c((vector("numeric",2)+0.8)), pt.bg=c("red","black"))
  rm(z)
  dev.off()
  
}


noBMI <- read.table('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CRP_model1_LLD_18072017.txt',sep="\t",header=T)
bmi <- read.table('/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/CRP_model2_LLD_18072017.txt',sep="\t",header=T)

QQplotFromPs(noBMI$P_VAL, '/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/qqplotNoBMI.pdf','noBMI')
QQplotFromPs(bmi$P_VAL, '/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/qqplotBMI.pdf','BMI')
