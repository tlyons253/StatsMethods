library(nimble)
setwd("C:/Users/rerussell/Desktop/Reno_workshop/Revised_workshop_materials/Module3_Chapter4")

black.ducks<-read.csv("blackduck.csv")
cond<-read.csv("cond.csv")  ###c
####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=nrow(black.ducks), left=black.ducks$Left,
                  right=black.ducks$Right, cond=cond[,1])
###data is the response variable, things that are being estimated
bugs4_1.data=list(censored=black.ducks$Censored)
initsgen <- function(){
  temp <- list(gamma0=-5,sd=1,rho=rep(0, 64),lhrcond=0)
  return(temp)
}

simInitsList <- initsgen()

######UGS3_1 Hierarchical CAR1 model

bugs4_1<-nimbleCode({ 
  gamma0~dflat()       
  rho[1]~dflat()
  rho[2]~dflat()
  for (i in 3:64) {
    diff[i]<-2*rho[i-1]-rho[i-2]
    rho[i]~dnorm(diff[i],tau)
  }
  sd~dunif(0,10)
  tau<-1/(sd*sd)
  lhrcond~dnorm(0,0.0000001)  # log hazard ratio
  hrcond<-exp(lhrcond)  # hazard ratio
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp( gamma0+rho[k]+lhrcond * (cond[j]-4.30))
      # Centered so cond = 4.30 is the baseline 
    }
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (i in 1:64) { # baseline, not population, survival
    UCH0[i]<-exp(gamma0+rho[i])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 							
}
)

######## rung the models
bug4_1.m<-nimbleMCMC(bugs4_1, data = bugs4_1.data, inits = simInitsList, constants=constants,
               monitors = c("gamma0", "S0", "UCH0", "sd", "hrcond", "lhrcond", "SLR"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE, WAIC=TRUE)
####Make a data frame to hold the data you want to plot
ex4_1.plot<-as.data.frame(bug4_1.m$summary$all.chains)
ex4_1.plot$type<-paste(sapply(strsplit(rownames(ex4_1.plot),split= "0"),'[',1), "COV", sep="")

ex4_1.plot<-cbind(ex4_1.plot[ex4_1.plot$type=="SCOV",], ex4_1.plot[ex4_1.plot$type=="UCHCOV",])
names(ex4_1.plot)[1:6]<-paste("S", names(ex4_1.plot)[1:6], sep=".")
names(ex4_1.plot)[7:12]<-paste("UCH", names(ex4_1.plot)[7:12], sep=".")
library(ggplot2)
###############Plot the data
ggplot(ex4_1.plot, aes(x=rep(seq(1,64)), group=S.type, col=S.type))+ 
  geom_line(aes(y = S.Median))+ 
  geom_line(aes(y = UCH.Median*10, group=UCH.type), col="magenta")+
  scale_y_continuous(limits = c(0, 1), name="Survival")+ theme_bw()+
  geom_ribbon(data=ex4_1.plot,aes(ymin=ex4_1.plot$'S.95%CI_low' ,ymax=ex4_1.plot$'S.95%CI_upp',col=ex4_1.plot$S.type, linetype=NA),alpha=0.1, fill="gray50")+
  scale_y_continuous(sec.axis = sec_axis(~.*1/10, name = "UCH"))+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Survival Probability")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+theme(legend.position="none")


