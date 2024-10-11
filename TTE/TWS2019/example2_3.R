#################16-Step hazard model
###library(nimble)
###get the data
###duck.data<-read.csv("blackduck.csv", header=T)
####constants are elements of you data that do not c
###Step Hazard Model

constants <- list(records=nrow(duck.data), left=duck.data$Left,
                   right=duck.data$Right,
                  lookup=c(rep(seq(1,16), each=4)))
###data is the response variable, things that are being estimated
SHmodel16.data=list(censored=duck.data$Censored)


initsgen <- function(){
  temp <- list(gamma=runif(16,-10, -1))
  return(temp)
}

#inits <- list(gamma=rep(-5,16))

##Create initial values
simInitsList <- initsgen()



####BUGS2_2  3-Step hazard model
SH.model16<-nimbleCode( { 
                       for (i in 1:16) {
                         gamma[i]~dflat()
                         }  # 16-segment log hazard model  
                       for (j in 1:records) {
                         for (k in left[j]:(right[j]-1)) {
                           UCH[j,k] <- exp(gamma[lookup[k]]) 
                         }
                         censored[j] ~ dbern(exp(-sum(UCH[j,left[j]:(right[j]-1)])))			
                       }
                       
                       for (i in 1:64) {
                         UCH0[i]<-exp(gamma[lookup[i]])
                         CH0[i]<-sum(UCH0[1:i])
                         S0[i]<-exp(-CH0[i])
                       } 
                       
                     }
                     
)
###############################Run the Model
SH.NimModel16<- nimbleModel(code=SH.model16, constants=constants, data=SHmodel16.data, inits=simInitsList, check = FALSE)
m3<-nimbleMCMC(SH.NimModel16, data = SHmodel.data16, inits = initsgen(),
               monitors = c("S0", "UCH0", "CH0", "gamma"), thin = 1,
               niter = 30000, nburnin = 10000, nchains = 3,
               summary = TRUE, WAIC = TRUE)

##############################################Make Plots#######################################
ex2_3.plot<-as.data.frame(m3$summary$all.chains)
ex2_3.plot$type<-paste(sapply(strsplit(rownames(ex2_3.plot),split= "0"),'[',1), 3, sep="")
ex2_3.plot<-ex2_3.plot[-grep("gamma", rownames(ex2_3.plot)),]


#############combine data from 3 models to compare
plot.ex2_3a<-rbind(ex2_1.plot[which(ex2_1.plot$type%in%c("S")),], 
                  ex2_2.plot[which(ex2_2.plot$type%in%c("S2")),],
                  ex2_3.plot[which(ex2_3.plot$type%in%c("S3")),])

ggplot(plot.ex2_3a, aes(x=rep(seq(1,64),3), y=Median, group=type, col=type))+ 
  geom_line(aes(col=type))+scale_y_continuous(limits = c(0.4, 1))+ theme_bw()+
  geom_ribbon(data=plot.ex2_3a,aes(ymin=plot.ex2_3a$'95%CI_low' ,ymax=plot.ex2_3a$'95%CI_upp',linetype=NA),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

################compare UCH
plot.ex2_3b<-rbind(ex2_1.plot[which(ex2_1.plot$type%in%c("UCH")),], 
                  ex2_2.plot[which(ex2_2.plot$type%in%c("UCH2")),],
                  ex2_3.plot[which(ex2_3.plot$type%in%c("UCH3")),])

ggplot(plot.ex2_3b, aes(x=rep(seq(1,64),3), y=Median, group=type, col=type))+ 
  geom_line(aes(col=type))+scale_y_continuous(limits = c(0, .075))+ theme_bw()+
  geom_ribbon(data=plot.ex2_3b,aes(ymin=plot.ex2_3b$'95%CI_low' ,ymax=plot.ex2_3b$'95%CI_upp',linetype=NA),alpha=0.1, fill="gray50")+
  scale_x_continuous(limits=c(0,64), breaks=seq(0,64, by=4))+
  theme(panel.grid.major = element_line(colour = "white"))+
  theme(plot.title = element_text(lineheight=.8,vjust=-0.8, face="plain", size=10))+
  scale_fill_manual("",values="gray90")+ labs(x="days" , y = "Outcome")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

