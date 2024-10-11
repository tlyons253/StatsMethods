library(nimble)
###############################Write a function
S.CAR<- nimbleFunction(
  run = function(records = double(0),
                 left=double(1),
                 right = double(1),
                 gamma0 = double(0), rho=double(1))
  {
    SLR <- nimNumeric(records)
    UCH <-nimMatrix(value = 0,nrow = records, ncol = 64)
    UCH0<-nimNumeric(64)
    
    for (j in 1:records) {
      for (k in left[j]:(right[j]-1)) {
        UCH[j,k] <- exp(gamma0+rho[k])
      }
      # total prob of survivingy
      SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)]))
      
    }
    
    returnType(double(1))
    return(SLR[1:records])
  })
########################################Model
SH.modelCAR<-nimbleCode( { 
  gamma0~dflat()       
  rho[1]~dflat()
  
  for (i in 2:64) {
    rho[i]~dnorm(rho[i-1],30)
    } # ICAR1 model
  
  SLR[1:records]<-S.CAR(records=records,
                           left = left[1:records],
                           right = right[1:records],
                           gamma0=gamma0, rho=rho[1:64])  
  
  for (j in 1:records) {
    censored[j] ~ dbern(SLR[j])			
  }
  
  for (i in 1:64) {
    UCH0[i]<-exp(gamma0+rho[i])
    CH0[i]<-sum(UCH0[1:i])
    S0[i]<-exp(-CH0[i])
  } 
  
})
  

##########################################Run the model
constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right)

###data is the response variable, things that are being estimated
SHmodelCAR.data= list(censored=duck.data$Censored)

initsgen <- function(){
  temp <- list(gamma0=runif(1,-10, -1), rho = rep(0, 64))
  return(temp)
}

#inits <-  list(gamma0=-5, rho=c(rep(0,64)))

##Create initial values
simInitsList <- initsgen()




SH.NimModelCAR<- nimbleModel(code=SH.modelCAR, constants=constants, data=SHmodelCAR.data, inits=simInitsList, 
                             check = FALSE)

a<-proc.time()
m5<-nimbleMCMC(SH.NimModelCAR, data = SHmodelCAR.data, inits = initsgen(),
               monitors = c("gamma0", "S0", "UCH0"), thin = 1,
               niter = 3000, nburnin = 2000, nchains = 3,
               summary = TRUE)
proc.time()-a
##############################################Make Plots#######################################
ex2_5.plot<-as.data.frame(m5$summary$all.chains)
ex2_5.plot$type<-paste(sapply(strsplit(rownames(ex2_5.plot),split= "0"),'[',1), 5, sep="")

plot.ex2_5<-cbind(ex2_5.plot[which(ex2_5.plot$type%in%c("S5")),], 
                  ex2_5.plot[which(ex2_5.plot$type%in%c("UCH5")),])
names(plot.ex2_5)[1:6]<-paste("S", names(plot.ex2_5)[1:6], sep=".")
names(plot.ex2_5)[7:12]<-paste("UCH", names(plot.ex2_5)[7:12], sep=".")
##############################
ggplot(plot.ex2_5, aes(x=rep(seq(1,64),1)))+ 
  geom_line(aes(y = S.Median), colour = "black")+ geom_line(aes(y = UCH.Median*10), colour = "magenta")+
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


