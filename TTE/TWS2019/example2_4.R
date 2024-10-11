library(nimble)
###get the data
###duck.data<-read.csv("blackduck.csv", header=T)
####creating an index, each unique day of death gets an index from 1 to 16 because there are 16 of them
###censored intervals get the index 17

lookup.data<-merge(data.frame(Day=seq(1,64)), 
                  data.frame(Day=sort(unique(duck.data$Left[which(duck.data$Censored==0)])), index=seq(1,16)), by="Day", all.x=T)
lookup.data$index[which(is.na(lookup.data$index))]<-17
 
###step hazard model (KM with 16 steps)

constants <- list(records=nrow(duck.data), left=duck.data$Left,
                  right=duck.data$Right,lookup=lookup.data$index)
 
###data is the response variable, things that are being estimated
SHmodelKM.data= list(censored=duck.data$Censored)

initsgen <- function(){
  temp <- list(gamma=runif(17,-10, -1))
  return(temp)
}

#inits <- list(gamma=c(rep(-5,16), NA))

##Create initial values
simInitsList <- initsgen()
 
###NEW write a nimble function to facilitate the running of the model
survival<- nimbleFunction(
  run = function(records = double(0),
                 left = double(1),
                 right = double(1),
                 gamma = double(1), 
                 lookup=double(1))
    {
    SLR <- nimNumeric(records)
    UCH <-nimMatrix(value = 0,nrow = records, ncol = 64)
    UCH0<-nimNumeric(64)
    
    for (j in 1:records) {
      for (k in left[j]:(right[j]-1)) {
        UCH[j,k] <- exp(gamma[lookup[k]])#+ age.effect[k]+ beta1*z[j] + beta2*sex[j])#
        
      }
      # total prob of survivingy
      SLR[j] <- exp(-sum(UCH[j,left[j]:(right[j]-1)]))
     
    }
    
    returnType(double(1))
    return(SLR[1:records])
  })

####BUGS2_2  step hazard model (KM with 16 steps)

SH.modelKM<-nimbleCode({ 
                          for (i in 1:16) {
                            gamma[i]~dflat()
                            }    
                          gamma[17]<- -99.9  # constrain hazard to almost zero
                          SLR[1:records]<-survival(records=records,
                                                           left = left[1:records],
                                                           right = right[1:records],
                                                          gamma=gamma[1:16], lookup=lookup[1:64])
                         
                          for(i in 1:64){
                            UCH0[i]<-exp(gamma[lookup[i]])
                            CH0[i]<-sum(UCH0[1:i])
                            S0[i]<-exp(-CH0[i])
                          }
                          
                          for (j in 1:records) {
                            censored[j] ~ dbern(SLR[j])			
                          }
                          
}
)

SH.NimModelKM<- nimbleModel(code=SH.modelKM, constants=constants, data=SHmodelKM.data, inits=simInitsList, check = FALSE)
m4<-nimbleMCMC(SH.NimModelKM, data = SHmodel.dataKM, inits = initsgen(),
               monitors = c("gamma", "S0", "UCH0"), thin = 1,
               niter = 30000, nburnin = 20000, nchains = 3,
               summary = TRUE)

##############################################Make Plots#######################################
ex2_4.plot<-as.data.frame(m4$summary$all.chains)
ex2_4.plot$type<-paste(sapply(strsplit(rownames(ex2_4.plot),split= "0"),'[',1), 4, sep="")


plot.ex2_4<-cbind(ex2_4.plot[which(ex2_4.plot$type%in%c("S4")),], 
                   ex2_4.plot[which(ex2_4.plot$type%in%c("UCH4")),])

names(plot.ex2_4)[1:6]<-paste("S", names(plot.ex2_4)[1:6], sep=".")
names(plot.ex2_4)[7:12]<-paste("UCH", names(plot.ex2_4)[7:12], sep=".")
###################################################################################################
ggplot(plot.ex2_4, aes(x=rep(seq(1,64),1), ))+ 
  geom_line(aes(y=S.Median), color="black")+ geom_line(aes(y = UCH.Median*10), color="magenta")+
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


rm(list=ls(pattern = "str"))


 