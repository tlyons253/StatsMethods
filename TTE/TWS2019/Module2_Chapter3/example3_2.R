library(nimble)
setwd("C:/Users/rerussell/Desktop/Reno_workshop/Revised_workshop_materials/Module2_Chapter3")

duck.data<-read.csv("blackduck.csv", header=T)
####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right)
###data is the response variable, things that are being estimated
CAR1.data=list(censored=duck.data$Censored)


initsgen <- function(){
  temp <- list(gamma0=-5,sd=1,rho=rep(0, 64))
  return(temp)
}

simInitsList <- initsgen()
######UGS3_1 Hierarchical CAR1 model

CAR.model2<-nimbleCode({ 
  gamma0~dflat()       
  rho[1]~dflat()
  rho[2]~dflat()
  for (i in 3:64) {
    diff[i]<-2*rho[i-1]-rho[i-2]
    rho[i]~dnorm(diff[i],tau) # ICAR2 model
  }
  sd~dunif(0,10)
  tau<-1/(sd*sd)
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp( gamma0+rho[k])
    }
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (i in 1:63) {
    UCH0[i]<-exp(gamma0+rho[i])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 							
})



CARm2<-nimbleMCMC(CAR.model2, data = CAR1.data, inits = simInitsList, constants=constants,
               monitors = c("gamma0", "S0", "UCH0", "sd", "SLR"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE, WAIC=TRUE)
####Plot of the standard deviation
plot.data<-as.data.frame(CARm2$summary$all.chains)
plot.data$type<-paste(sapply(strsplit(rownames(plot.data),split= "0"),'[',1), "CAR", sep="")
plot.data<-plot.data[which(plot.data$type%in%c("SCAR", "UCHCAR")),]

example3_2.data<-cbind(plot.data[1:63,], plot.data[64:126,])
names(example3_2.data)[1:6]<-paste("S", names(example3_2.data)[1:6], sep=".")
names(example3_2.data)[7:12]<-paste("UCH", names(example3_2.data)[7:12], sep=".")
plot.ex3_2<-example3_2.data

ggplot(plot.ex3_2, aes(x=rep(seq(1,63),1)))+ 
  geom_line(aes(y = S.Median), col="cadetblue")+ 
  geom_line(aes(y = UCH.Median*10), col="chartreuse")+
  scale_y_continuous(limits = c(0, 1), name="Survival")+ theme_bw()+
  geom_ribbon(data=plot.ex3_2,aes(ymin=plot.ex3_2$'S.95%CI_low' ,ymax=plot.ex3_2$'S.95%CI_upp',col=plot.ex3_2$S.type, linetype=NA),alpha=0.1, fill="gray50")+
  scale_y_continuous(sec.axis = sec_axis(~.*1/10, name = "UCH"))+
  scale_x_continuous(limits=c(0,63), breaks=seq(0,63, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Survival Probability")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+theme(legend.position="none")


