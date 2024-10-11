library(nimble)
library(ggplot2)
###get the data
setwd("C:/Users/dwalsh.GS/Documents/Time2event_TWS_19/Examples/Module1_Chapter2")
duck.data<-read.csv("blackduck.csv", header=T)
####constants are elements of you data that do not change, covariates, number of observations etc.

constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right)
###data is the response variable, things that are being estimated
S0model.data=list(censored=duck.data$Censored)
###Initial Values

initsgen <- function(){
  temp <- list(gamma=runif(1,-10, -1))
  return(temp)
}


##Create initial values
simInitsList <- initsgen()

#############################Model Code

S0.model<-nimbleCode({ 
  gamma~dflat() 
  #the j loop steps through all the records
  #the next statements compute S(R | L)
  #the last statements constructs the likelihood, determining whether
  #the term should be an S(R | L) or M(R | L) 
  
  for (j in 1:records) { 
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(gamma)
    } 
    
    S[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(S[j])			
  }
  
  for (i in 1:64) { #compute the survival function S0
    UCH0[i]<-exp(gamma)
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 		
})

########Run the model
S0.NimModel<- nimbleModel(code=S0.model, constants=constants, data=S0model.data, inits=simInitsList, check = FALSE)
m<-nimbleMCMC(S0.NimModel, data = S0model.data, inits = initsgen(),
              monitors = c("S0", "UCH0", "CH0", "gamma"), thin = 1,
              niter = 30000, nburnin = 10000, nchains = 3,
              summary = TRUE, WAIC = TRUE)

###############Check model Diagnostics
##Trace plot
plot(m$samples$chain1[,"gamma"], type="l")
lines(m$samples$chain2[,"gamma"], type="l", col="violet")
lines(m$samples$chain3[,"gamma"], type="l", col="darkgreen")

####diagnostics
###use coda to plot
library(coda)
chain1<-as.mcmc(m$samples$chain1)
chain2<-as.mcmc(m$samples$chain2)
chain3<-as.mcmc(m$samples$chain3)

###creat a list use coda to check diagnostics
out.list<-mcmc.list(chain1, chain2, chain3)
gelman.diag(out.list[,"gamma"])
effectiveSize(out.list[,"gamma"]) #ESS all chains
lapply(out.list[,"gamma"], effectiveSize)  #ESS per chain

##########clean up the output and label it  
ex2_1.plot<-as.data.frame(m$summary$all.chains)
ex2_1.plot$type<-sapply(strsplit(rownames(ex2_1.plot),split= "0"),'[',1)
ex2_1.plot<-ex2_1.plot[which(ex2_1.plot$type%in%c("S", "UCH")),]


###Plot it don't plotcrastinate!
 library(ggplot2)
ggplot(ex2_1.plot, aes(x=rep(seq(1,64),2), y=Median, group=type, col=type))+ 
  geom_line(aes(col=type))+scale_y_continuous(limits = c(0, 1))+ theme_bw()+
  geom_ribbon(data=ex2_1.plot,aes(ymin=ex2_1.plot$'95%CI_low' ,ymax=ex2_1.plot$'95%CI_upp',linetype=NA),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))



###More efficient model specification
#############################Model Code

S0.model<-nimbleCode({ 
  gamma~dflat() 
  #the j loop steps through all the records
  #the next statements compute S(R | L)
  #the last statements constructs the likelihood, determining whether
  #the term should be an S(R | L) or M(R | L) 
  
  for (j in 1:records) { 
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(gamma)
    } 
      censored[j] ~ dbern(exp(-sum(UCH[j,left[j]:(right[j]-1)])))			
  }
  
  for (i in 1:64) { #compute the survival function S0
    UCH0[i]<-exp(gamma)
    S0[i]<-exp(-sum(UCH0[1:i]))
  } 		
})



########Run the model
S0.NimModel<- nimbleModel(code=S0.model, constants=constants, data=S0model.data, inits=simInitsList, check = FALSE)
m<-nimbleMCMC(S0.NimModel, data = S0model.data, inits = initsgen(),
              monitors = c("S0", "UCH0", "gamma"), thin = 1,
              niter = 30000, nburnin = 10000, nchains = 3,
              summary = TRUE, WAIC = TRUE)

