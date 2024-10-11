###Step Hazard Model
###library(nimble)
###get the data
####duck.data<-read.csv("blackduck.csv", header=T)
####constants are elements of you data that do not c
constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right,
                  lookup=c(rep(c(1,2,3), each=21),3))
###data is the response variable, things that are being estimated
SHmodel.data=list(censored=duck.data$Censored)


initsgen <- function(){
  temp <- list(gamma=runif(3,-10, -1))
  return(temp)
}

#inits <- list(gamma=rep(0,3))

##Create initial values
simInitsList <- initsgen()

####BUGS2_2  3-Step hazard model
SH.model<-nimbleCode({ 
  ##priors
  for (i in 1:3) {
    gamma[i]~dflat()
    }  

  # 3-segment log hazard model  
  for (j in 1:records) {
    for (k in left[j]:(right[j]-1)) {
      UCH[j,k] <- exp(gamma[lookup[k]]) # "nested indexing"
    }
    S[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)])) 
    censored[j] ~ dbern(S[j])			
  }
  ###Summary statistics
  for (i in 1:64) {
    UCH0[i]<-exp(gamma[lookup[i]])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  }
})
###############Run the model
SH.NimModel<- nimbleModel(code=SH.model, constants=constants, data=SHmodel.data, inits=simInitsList, check = FALSE)
m2<-nimbleMCMC(SH.NimModel, data = SHmodel.data, inits = initsgen(),
              monitors = c("S0", "UCH0", "CH0", "gamma"), thin = 1,
              niter = 30000, nburnin = 10000, nchains = 3,
              summary = TRUE, WAIC = TRUE)

####################Plot the data
ex2_2.plot<-as.data.frame(m2$summary$all.chains)
ex2_2.plot$type<-paste(sapply(strsplit(rownames(ex2_2.plot),split= "0"),'[',1), 2, sep="")


####Combine results from model 1 and 2 to plot 
plot.ex2_2<-rbind(ex2_1.plot[which(ex2_1.plot$type%in%c("S","UCH")),], 
                  ex2_2.plot[which(ex2_2.plot$type%in%c("S2","UCH2")),])

ggplot(plot.ex2_2, aes(x=rep(seq(1,64),4), y=Median, group=type, col=type))+ 
  geom_line(aes(col=type))+scale_y_continuous(limits = c(0, 1))+ theme_bw()+
  geom_ribbon(data=plot.ex2_2,aes(ymin=plot.ex2_2$'95%CI_low' ,ymax=plot.ex2_2$'95%CI_upp',linetype=NA),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

 