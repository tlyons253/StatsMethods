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

CAR.model<-nimbleCode({ 
  gamma0~dflat()       
  rho[1]~dflat()
  for (i in 2:64) {
    rho[i]~dnorm(rho[i-1],tau)
    } # hyperparameter
  sd~dunif(0,10) # hyperprior
  tau<-1/(sd*sd)  # dnorm() parameterized by precision, not var.
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp( gamma0+rho[k])
    }
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (i in 1:64) {
    UCH0[i]<-exp(gamma0+rho[i])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 							
})


CAR.NimModel<- nimbleModel(code=CAR.model, constants=constants, data=CAR1.data, inits=simInitsList, check =TRUE)
CARm1<-nimbleMCMC(CAR.NimModel, data = CAR1.data, inits = simInitsList,
               monitors = c("gamma0", "S0", "UCH0", "sd", "SLR"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE, WAIC=TRUE)
####Plot of the standard deviation
plot(density(c(CARm1$samples[[1]][,"sd"],CARm1$samples[[2]][,"sd"], CARm1$samples[[3]][,"sd"])), main="sd")

ex3_1.plot<-as.data.frame(CARm1$summary$all.chains)
ex3_1.plot$type<-paste(sapply(strsplit(rownames(ex3_1.plot),split= "0"),'[',1), "CAR", sep="")
ex3_1.plot<-ex3_1.plot[which(ex3_1.plot$type%in%c("SCAR", "UCHCAR")),]

ex3_1.plot<-cbind(ex3_1.plot[1:64,], ex3_1.plot[65:128,])
names(ex3_1.plot)[1:6]<-paste("S", names(ex3_1.plot)[1:6], sep=".")
names(ex3_1.plot)[7:12]<-paste("UCH", names(ex3_1.plot)[7:12], sep=".")
plot.ex3_1<-rbind(ex3_1.plot,plot.ex2_4)

ggplot(plot.ex3_1, aes(x=rep(seq(1,64),2), group=S.type, col=S.type))+ 
  geom_line(aes(y = S.Median))+ 
  geom_line(aes(y = UCH.Median*10, group=UCH.type), col="magenta")+
  scale_y_continuous(limits = c(0, 1), name="Survival")+ theme_bw()+
  geom_ribbon(data=plot.ex3_1,aes(ymin=plot.ex3_1$'S.95%CI_low' ,ymax=plot.ex3_1$'S.95%CI_upp',col=plot.ex3_1$S.type, linetype=NA),alpha=0.1, fill="gray50")+
  scale_y_continuous(sec.axis = sec_axis(~.*1/10, name = "UCH"))+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Survival Probability")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+theme(legend.position="none")


