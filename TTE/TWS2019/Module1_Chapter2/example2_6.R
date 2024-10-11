library(nimble)
setwd("C:/Users/rerussell/Desktop/Reno_workshop/Revised_workshop_materials/Module1_Chapter2")
duck.data<-read.csv("blackduck.csv", header=T)
###Can do this inside the code as well or write a function.
####create adjacency matrix
days<-64
adj<-NULL
adj[1] <- 2
for(t in 2:(days-1)) {
  adj[2+(t-2)*2] <- t-1
  adj[3+(t-2)*2] <- t+1
}	
adj[(days-2)*2 + 2] <- days-1	

#################weights
wts<-NULL
wts[1] <- 1
for(t in 2:(days-1)) {
  wts[2+(t-2)*2] <- 1
  wts[3+(t-2)*2] <- 1
}	
wts[(days-2)*2 + 2] <- 1
########################
num<-NULL
num[1] <- 1
for(t in 2:(days-1)) {
  num[t] <- 2
}	
num[days] <- 1 					

################

constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right, days=64, adj=adj, wts=wts, num=num)

###data is the response variable, things that are being estimated
CARmodel.data= list(censored=duck.data$Censored)

inits <- list(gamma0 = -4.54,rho = as.vector(runif(64, -1,1)))
###rho<-dcar_normal(adj, wts, num, tau,c=2, zero_mean=1)  

CAR.model<-nimbleCode({ 
  tau<-20
  
  rho[1:days] ~ dcar_normal(adj[1:126], wts[1:126], num[1:days], tau, zero_mean=1) 
  
  gamma0~dflat()
  
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
}
  
)
CAR.nimble<- nimbleModel(code=CAR.model, constants=constants, data=CARmodel.data, inits=inits, check = FALSE)
m6<-nimbleMCMC(CAR.model, data = CARmodel.data, inits = inits, constants=constants,
               monitors = c("gamma0", "S0", "UCH0"), thin = 10,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE)
###############################################################################################
ex2_6.plot<-as.data.frame(m6$summary$all.chains)
ex2_6.plot$type<-paste(sapply(strsplit(rownames(ex2_6.plot),split= "0"),'[',1), 6, sep="")

plot.ex2_6<-cbind(ex2_6.plot[which(ex2_6.plot$type%in%c("S6")),], 
                  ex2_6.plot[which(ex2_6.plot$type%in%c("UCH6")),])
names(plot.ex2_6)[1:6]<-paste("S", names(plot.ex2_6)[1:6], sep=".")
names(plot.ex2_6)[7:12]<-paste("UCH", names(plot.ex2_6)[7:12], sep=".")

ggplot(plot.ex2_6, aes(x=rep(seq(1,64),1), ))+ 
  geom_line(aes(y=S.Median), color="steelblue")+ geom_line(aes(y = UCH.Median*10), color="turquoise")+
  scale_y_continuous(limits = c(0, 1), name="Survival")+ theme_bw()+
  scale_y_continuous(sec.axis = sec_axis(~.*1/10, name = "UCH"))+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

