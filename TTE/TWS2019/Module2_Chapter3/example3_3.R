library(nimble)

new.ducks<-read.csv("new_ducks.csv")
entry<-read.csv("entry.csv")
####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=nrow(new.ducks), left=new.ducks$left,
                  right=new.ducks$right, entry=entry[,1], subjects=216)
###data is the response variable, things that are being estimated
bugs3_3.data=list(censored=new.ducks$censored)



initsgen <- function(){
  temp <- list(gamma0=-10,sd=.1,rho=rep(0, 37))
  return(temp)
}

simInitsList <- initsgen()

######UGS3_1 Hierarchical CAR1 model

bugs3_3<-nimbleCode({ 
  gamma0~dflat()       
  rho[1]~dflat()
  rho[2]~dflat()
  for (i in 3:37) {
    diff[i]<-2*rho[i-1]-rho[i-2]
    rho[i]~dnorm(diff[i],tau)
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
  
  for (i in 1:37) {
    UCH0[i]<-exp(gamma0+rho[i])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 
  for (i in 1:subjects) {
    UCHW[i,1]<-0.0
    for (j in 2:entry[i]) {UCHW[i,j]<-exp(gamma0+rho[j])} 
    CHW[i]<-sum(UCHW[i,1:entry[i]])
    SW[i]<-exp(-CHW[i]) # Survival until entry
    Virtual[i]<-1/SW[i] # Virtual subjects (includes observed)
  }   
  Total<-sum(Virtual[1:subjects]) #"Backcast" total
  
})

hist(entry[,1], breaks=seq(1,30), xlab="entry date") 
####how many ducks entered the data set after day 1
length(which(entry[,1]>1))
######## rung the models
bug3_3.m<-nimbleMCMC(bugs3_3, data = bugs3_3.data, inits = simInitsList, constants=constants,
               monitors = c("gamma0", "S0", "UCH0", "sd", "SLR"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE, WAIC=TRUE)
####Make a data frame to hold the data you want to plot
ex3_3.plot<-as.data.frame(bug3_3.m$summary$all.chains)
ex3_3.plot$type<-paste(sapply(strsplit(rownames(ex3_3.plot),split= "0"),'[',1), "CAR", sep="")

ex3_3.plot<-cbind(ex3_3.plot[which(ex3_3.plot$type=="SCAR"),], ex3_3.plot[which(ex3_3.plot$type=="UCHCAR"),])
names(ex3_3.plot)[1:6]<-paste("S", names(ex3_3.plot)[1:6], sep=".")
names(ex3_3.plot)[7:12]<-paste("UCH", names(ex3_3.plot)[7:12], sep=".")


###############Plot the data
ggplot(ex3_3.plot, aes(x=rep(seq(1,37),1)))+ 
  geom_line(aes(y = S.Median), col="thistle", lwd=1)+ 
  geom_line(aes(y = UCH.Median), col="turquoise", lwd=1)+
  scale_y_continuous(limits = c(0, 1), name="Survival")+ theme_bw()+
  geom_ribbon(data=ex3_3.plot,aes(ymin=ex3_3.plot$'S.95%CI_low' ,ymax=ex3_3.plot$'S.95%CI_upp'),alpha=0.1, fill="gray50")+
  geom_ribbon(data=ex3_3.plot,aes(ymin=ex3_3.plot$'UCH.95%CI_low' ,ymax=ex3_3.plot$'UCH.95%CI_upp'),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,37), breaks=seq(0,37, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Survival Probability")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+theme(legend.position="none")


