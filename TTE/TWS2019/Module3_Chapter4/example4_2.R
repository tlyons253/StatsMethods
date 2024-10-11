library(nimble)
setwd("C:/Users/rerussell/Desktop/Reno_workshop/Revised_workshop_materials/Module3_Chapter4")
load("bwt.RData")

####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=bwt.data$records,subjects=bwt.data$subjects, left=bwt.data$left,
                  right=bwt.data$right, droad=bwt.data$droad, entry=bwt.data$entry, xentry=bwt.data$xentry)
###data is the response variable, things that are being estimated
Cov2.data=list(censored=bwt.data$censored)



initsgen <- function(){
  temp <- list(gamma0=-5,sd=1,rho=rep(0, 380),lhrdroad=0.01)
  return(temp)
}

simInitsList <- initsgen()

######Covariate model






Cov.model2<-nimbleCode({ 
  
  gamma0~dflat()       
  rho[1]~dflat()
  rho[2]~dflat()
  for (i in 3:36) {
    diff[i]<-2*rho[i-1]-rho[i-2]
    rho[i]~dnorm(diff[i],tau)
  }
  sd~dunif(0,10)
  tau<-1/(sd*sd) 
  lhrdroad~dnorm(0,0.0000001)
  hrdroad<-exp(lhrdroad)      
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(gamma0+rho[k]+lhrdroad*(droad[j]-20.95))
    }
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (j in 1:36) {  #baseline survival function estimate
    UCH0[j]<-exp(gamma0+rho[j])
    CH0[j]<-sum(UCH0[1:j])
    S0[j]<-exp(-CH0[j])
  } 		
  
  for (i in 1:subjects) {  #compute virtual subjects and weights
    UCHW[i,1]<-0.0
    for (j in 2:entry[i]) {
      UCHW[i,j]<-
        exp(gamma0+rho[j] +lhrdroad*(xentry[i] - 20.95))
    } 
    CHW[i]<-sum(UCHW[i,1:entry[i]])
    SW[i]<-exp(-CHW[i])
    Virtual[i]<-1/SW[i]
    W[i]<-Virtual[i]/Total
    WX[i]<-W[i]*xentry[i]
  }   
  Total<-sum(Virtual[1:subjects])  
  AverageX<-subjects*mean(WX[1:subjects]) # "Back-cast" average
  #  distance
  
  for (i in 1:subjects) { #model-based survival curve estimate
    for (j in 1:37) {
      UCHi[i,j]<-exp(gamma0+rho[j] +lhrdroad*(xentry[i] - 20.95))
      CHi[i,j]<-sum(UCHi[i,1:j])
      Si[i,j]<-W[i]*exp(-CHi[i,j])	
    }			  
  }	 
  for (j in 1:37) {Smodel[j]<-sum(Si[1:subjects,j])}
  
}
)

S0.NimModel<- nimbleModel(code=Cov.model2, constants=constants, data=Cov2.data, inits=simInitsList, check = FALSE)

Covm2<-nimbleMCMC(Cov.model2, data = Cov2.data, inits = simInitsList, constants=constants,
               monitors = c("gamma0", "S0", "UCH0","Smodel", "sd", "hrdroad","lhrdroad", "SLR"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE, WAIC=TRUE)
####Plot of the standard deviation
plot.data<-as.data.frame(Covm2$summary$all.chains)
SO<-plot.data[grep("S0", row.names(plot.data)),]
UCH0<-plot.data[grep("UCH0", row.names(plot.data)),]
Smodel<-plot.data[grep("Smodel", row.names(plot.data)),]
####combine all posterior sample
hrdroad<-c(Covm2$samples[[1]][,"hrdroad"],Covm2$samples[[2]][,"hrdroad"],Covm2$samples[[3]][,"hrdroad"])
plot(density(hrdroad), main="Hazard")
lhrdroad<-c(Covm2$samples[[1]][,"lhrdroad"],Covm2$samples[[2]][,"lhrdroad"],Covm2$samples[[3]][,"lhrdroad"])
plot(density(lhrdroad))


ex4_2.data<-rbind(SO[1:36,], Smodel[1:36,])
ex4_2.data$type<-c(rep("SO",36), rep("Smodel",36))

ggplot(ex4_2.data, aes(x=rep(seq(1,36),2), y=Median, group=type, col=type))+ 
  geom_line(aes(col=type))+scale_y_continuous(limits = c(0, 1))+ theme_bw()+
  geom_ribbon(data=ex4_2.data,aes(ymin=ex4_2.data$'95%CI_low',ymax=ex4_2.data$'95%CI_upp',linetype=NA),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,36), breaks=seq(0,36, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

