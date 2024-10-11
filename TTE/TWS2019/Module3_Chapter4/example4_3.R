library(nimble)
setwd("C:/Users/rerussell/Desktop/Reno_workshop/Revised_workshop_materials/Module3_Chapter4")
black.ducks<-read.csv("blackduck.csv", header=T)
cond<-read.csv("cond.csv", header=T)

####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=nrow(black.ducks), left=black.ducks$Left,
                   right=black.ducks$Right, cond=cond[,1])
###data is the response variable, things that are being estimated
Cov3.data=list(censored=black.ducks$Censored)



initsgen <- function(){
  temp <- list(gamma0=-5,sd=1,rho=rep(0, 64),sdtv=1, beta0=0, beta1=rep(0,64))
  return(temp)
}

simInitsList <- initsgen()


######Covariate model

Cov.model3<-nimbleCode({ 
  gamma0~dflat()       
  rho[1]~dflat()
  rho[2]~dflat()
  for (i in 3:64) {
    diff[i]<-2*rho[i-1]-rho[i-2]
    rho[i]~dnorm(diff[i],tau)
  }
  sd~dunif(0,10)
  tau<-1/(sd*sd)
  
  sdtv~dunif(0,10) # create another ICAR2 prior
  tautv<-1/(sdtv*sdtv)
  beta0~dflat()       
  beta1[1]~dflat()
  beta1[2]~dflat()
  for (i in 3:64) {
    betadiff[i]<-2*beta1[i-1]-beta1[i-2]
    beta1[i]~dnorm(betadiff[i],tautv)
  }
  
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp( gamma0+rho[k]+
                         (beta0 + beta1[k]) * (cond[j]-4.30))
    }
    SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (i in 1:64) { # time-varying hazard and log-hazard
    lhr[i]<-beta0+beta1[i]
    hr[i]<-exp(lhr[i])
  } 						
}
)

##S0.NimModel<- nimbleModel(code=Cov.model3, constants=constants, data=Cov3.data, inits=inits, check = FALSE)

Covm3<-nimbleMCMC(Cov.model3, data = Cov3.data, inits = simInitsList, constants=constants,
                  monitors = c("lhr", "hr", "gamma0"), thin = 10,
                  niter = 100000, nburnin = 30000, nchains = 3,
                  summary = TRUE, WAIC=TRUE)

####combine all posterior sample
hr<-c(Covm3$samples[[1]][,"hr[1]"],Covm3$samples[[2]][,"hr[1]"],Covm3$samples[[3]][,"hr[1]"])
plot(density(hr), main="Hazard")
lhr<-c(Covm3$samples[[1]][,"lhr[1]"],Covm3$samples[[2]][,"lhr[1]"],Covm3$samples[[3]][,"lhr[1]"])
plot(density(lhr), main="log-Hazard")

hr.lhr<-data.frame(Covm3$summary$all.chains[grep("hr", rownames(Covm3$summary$all.chains)),])
hr.lhr$type<-c(rep("hr",64), rep("lhr",64))
hr<-subset(hr.lhr, hr.lhr$type=="hr")

ggplot(hr, aes(x=rep(seq(1,64)), y=Mean))+ 
  geom_line(aes())+ theme_bw()+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Hazard Ratio")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))


