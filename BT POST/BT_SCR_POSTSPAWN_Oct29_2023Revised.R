RStudioGD() #this will set the plot device (dev) to the plot pane in rstudio
# Code provided by B. Gardner and modifed by J. Hightower and J. Raabe for spatially 
#explicit capture-recapture model for Little River
# Similar to examples in Gardner et al. 2010, Kery et al. 2010, Sollmann et al. 2011
library(jagsUI);library(nnet);library(MCMCvis);library(rjags)
library(reshape);library(tidyverse);library(mcmcOutput)
#setwd("C:/newlaptop/FISH/RapidRiverBullTrout/R_scripts/SCR_BT/")
setwd("~/patti1/Bull Trout POST")
AS10pre <- read.table("BT_SCR_POST1Knov4.txt",header=T)
AS10 <- AS10pre%>%##all weeks represented
  arrange(TagID,Week,RKM)%>%
  mutate(RKM = RKM/1000)#change from meters to km 
#mutate(RKM = 97.7629-RKM)#reverse the order so it starts at 1
AS10 <- round(AS10,1)  
melted.rkm <- melt(AS10, id=c("TagID","RKM")) 
tagid.week <- cast(melted.rkm, TagID ~ RKM ~ value, fill=0, length)
y <- tagid.week
#y[36,2,2]
#unique(AS10$TagID)
#sort(unique(AS10$RKM))
M <- length(table(AS10$TagID))    # Number of individuals
T <- length(table(AS10$Week))     # Number of periods (weeks) 
##so the end location will be 97.7629 for the reverse of this to start at 1
reach.loc <- c(96.7629,94.76,93.40,92.35,91.298,90.28,89.25,86.1122,82.0098,
               79.9855,78.9553,64.415,63.4115,61.2025,58,57,56,55,54,53,
               52,51,50,49,48,46,41,40,38,37,36,35,34,32,31,30,29,28,
               27,26,25,21,19,16,7,6,5,4,3,2,1)  
#reach.loc <- rev(reach.loc)
nreach <- length(reach.loc)  # weir + reaches in BT case (not antenna)

# Input and format data matrix:
# y = detection history (individual (i) x reach (j) x week(t)); first = river 
#immigration week, last = river emigration week, (these are counts).

##first week fish came into study
firstc<-read.csv("firstPOSTnov4.csv",stringsAsFactors = F)
firstb <- firstc %>%
  select(firstcap)
first <- as.numeric(firstb$firstcap)
#sort(unique(AS10$Week))
##last 
last <- c(13,8,19,5,18,16,4,15,9,13,14,11,8,19,12,19,19,19,16,11,4,19,19)

##make size length variable
fishsize <- firstc %>%
  select(FL,FW)%>%
  mutate_at(c("FL","FW"), ~(scale(.) %>% as.vector))
sizeFL <- as.numeric(fishsize$FL)

# flowA = average flow for 1992-1994 weighted by prop fish/year
flowA <- c(125.07713, 100.33387,  97.83228, 103.87914, 111.33861, 112.01538, 102.61095,
           90.61472,  97.32703,  92.42417,  90.88462,  92.84326,  94.09324, 108.09876,
           136.64214, 184.56233, 244.46015,303.66906, 398.57661)
##flowB = maximum flow - minimum flow within current week
flowB <- c(38.56755,  16.99577,  15.12686,  19.39704,  11.84777,  23.32176,  15.81213,
           50.52292,  34.17843,  21.99653,  23.29344,  25.34924,  23.05558,  30.09515,
           36.51174, 88.79031, 119.45745, 107.84188,  65.60447)
##flow C = change in previous week mean flow to current week mean flow
flowC <- c(-13.28950, -24.74326,  -2.50159,   6.04686,   7.45947,   0.67677,  -9.40443,
           -11.99623,   6.71231, -4.90286,  -1.53955,   1.95864,   1.24998,  14.00552,
           28.54338,  47.92019,  59.89782,  59.20891,  94.90761)
flow=(flowB-mean(flowB))/(sd(flowB))#standardize flow variable

##Temp variables created Feb15 see paper for context
#TempA is A) weekly maximum hourly water temperature avg 1992-1994
TempA <- c(15.796387, 14.338322, 12.902372, 12.528530, 10.248000,  8.648000,
           6.318000,  4.936000, 6.146000,  4.776000,  4.175947,  4.905678,
           5.290699,  5.901468,  7.233996,  7.767503,  8.828321, 10.014755, 10.566981)
#Temp B is weekly average daily maximum water temperature avg 1992-1994
TempB <- c(13.996688, 12.802443, 11.671784, 10.438367,  8.858000,  6.894857,  5.165571,
           3.610000,  4.276857,  3.622832,  2.864495,  3.625040,  3.874573,  4.241382,
           5.860427,  6.685168,  7.544143,8.440793,  9.117784)
# ##Temp C is weekly three-day moving average of the maximum daily temperature avg 1992-1994
TempC <- c(13.935234, 12.818215, 11.738988, 10.345602,  8.844389,  6.876333,  5.269389,
           3.680944, 4.296611,  3.595110,  2.850174,  3.687394,  3.979557,  4.147522,
           5.823203,  6.616615,  7.483100, 8.479822,  9.106690)
# ##Temp D is weekly average hourly temperature avg 1992-1994
TempD <- c(12.343923, 11.537544, 10.460930,  9.593330,  8.056595,  6.273119,
           4.588060,  3.113548, 3.745988,  3.156401,  2.413285,  3.035409,  3.170718,
           3.287662,  4.679939,  5.389000,  6.218292, 6.886377,  7.424620)
##Temp E is delta current - previous week average 1992-1994
TempE <- c(-1.131609, -0.806379, -1.076614, -0.867600, -1.536735, -1.783476, -1.685059,
           -1.474512,  0.632440, -0.589587, -0.743116,  0.622124,  0.135309,  0.116944,
           1.392277,  0.709061,  0.829292,  0.668085, 0.538243)
#temp=(TempB-mean(TempB))/(sd(TempB))#standardize flow variable
flowA=(flowA-mean(flowA))/(sd(flowA))#standardize flow variable
flowB=(flowB-mean(flowB))/(sd(flowB))#standardize flow variable
flowC=(flowC-mean(flowC))/(sd(flowC))#standardize flow variable
tempA=(TempA-mean(TempA))/(sd(TempA))#standardize flow variable
tempB=(TempB-mean(TempB))/(sd(TempB))#standardize flow variable
tempC=(TempC-mean(TempC))/(sd(TempC))#standardize flow variable
tempD=(TempD-mean(TempD))/(sd(TempD))#standardize flow variable
tempE=(TempE-mean(TempE))/(sd(TempE))#standardize flow variable

# test flow and temp correlations
# flowmets <- matrix(0,nrow=18, ncol = 8)
# for (i in 1:18) {
#   flowmets[i,1]<- flowA[i]
#   flowmets[i,2]<- flowB[i]
#   flowmets[i,3]<- flowC[i]
#   flowmets[i,4]<- TempA[i]
#   flowmets[i,5]<- TempB[i]
#   flowmets[i,6]<- TempC[i]
#   flowmets[i,7]<- TempD[i]
#   flowmets[i,8]<- TempE[i]
#   }
#  colnames(flowmets) <- c("flowA","flowB","flowC","TempA","TempB","TempC","TempD","TempE")
#  flowmets<-as.data.frame(flowmets)
#   to_test<-c("flowA","flowB","flowC","TempA","TempB","TempC","TempD","TempE")
#   cor(flowmets[,to_test])
##all flow vars correlated but not with temp

#         flowA        flowB       flowC        TempA       TempB       TempC       TempD      TempE
# flowA 1.0000000  0.908030989  0.87487883  0.181329701  0.19190239  0.18916194  0.12738607  0.4876063
# flowB 0.9080310  1.000000000  0.85470468  0.003772714  0.01636113  0.01238583 -0.04880912  0.5475822
# flowC 0.8748788  0.854704680  1.00000000 -0.089003398 -0.08106945 -0.08749481 -0.14548046  0.6902300
# TempA 0.1813297  0.003772714 -0.08900340  1.000000000  0.99730601  0.99673158  0.99519765 -0.2634882
# TempB 0.1919024  0.016361128 -0.08106945  0.997306014  1.00000000  0.99983602  0.99642502 -0.2642031
# TempC 0.1891619  0.012385828 -0.08749481  0.996731583  0.99983602  1.00000000  0.99662232 -0.2688198
# TempD 0.1273861 -0.048809118 -0.14548046  0.995197653  0.99642502  0.99662232  1.00000000 -0.3255345
# TempE 0.4876063  0.547582213  0.69023003 -0.263488174 -0.26420308 -0.26881980 -0.32553445  1.0000000

#week =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)#equals biweekly!!

#bring in detect file for aerial versus radio
detect<-read.csv("detectPostoct282023.csv",stringsAsFactors = F)

##make the raw detections plot
names(detect)

detectgg <- detect%>%
  mutate(RADIO_NO=as.character(RADIO_NO))
#make an abacus plot of raw detections
# RStudioGD()
# ggformat<-(
#   ggplot(detectgg,aes(x=twoweek,y=RADIO_NO,col = RELO_METHO))+
#     geom_jitter(size=1.25, shape = 20,width = 0.17,height =0.13) +
#     theme_bw()+
#     labs(title = "", x = "BiWeeks", y = "Fish")+
#     # scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "May 26","3" = "Jun 9",
#     #                 "5" = "Jun 23","7" = "Jul 7","9" = "Jul 22", "11" = "Aug 4","13" = "Aug 18",
#     #                 "15" = "Sept 1","17" = "Sept 15"), limits = c(1, 18))+
#     scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "Aug 19", 
#            "3" = "Sept 16", "5" = "Oct 14","7" = "Nov 11","9" = "Dec 9",
#           "11" = "Jan 6", "13" = "Feb 4","15" = "Mar 4","17" = "Apr 1"))+
#     
#     
#     theme( panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
#     theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=8,color="black"))+
#     theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=8,color="black"))+
#     theme(axis.title.y = element_text(size = rel(1.6), angle = 90))+
#     theme(axis.title.x = element_text(size = rel(1.6), angle = 00))+
#     theme(plot.margin = margin(0,1,0,0, "cm"))+
#     theme(legend.position = "none")
# )
# ggsave(filename="BTPOST1992_1994abacusplot.png", plot=ggformat, device="png",
#        height=6, width=7, units="in", dpi=800)

##make season variable for fish years in study NOO, use flow and temp for year variation
# season <- firstc %>%
#   select(YEAR)%>%
#   mutate(YEAR=recode(YEAR,"1992"="1","1993"="2"))
# year <- as.numeric(season$YEAR)
###################################################################
#####  WRITE TEXT FILE WITH JAGS SPECIFICATION   ###########
###################################################################
sink("ModelCJS.txt")
cat("
model { #model
# Priors
  lam0 ~ dunif(0.1,0.2)  # constant Weekly encounter rate
  alphaphi ~ dnorm(0,0.368)      #intercept for phi
  alphasig ~ dnorm(0,0.368)  # Intercept in sigma estimate
  #betasig ~ dnorm(0,0.368)  #Coefficient for flow in sigma estimate 
  betasig2 ~ dnorm(0,0.368)
  betaT2 ~ dnorm(0,0.368) 
  betaT3 ~ dnorm(0,0.368) 
  betaphiT ~ dnorm(0,0.368)     #beta for phi for temp
    #betaphipos ~ dnorm(0,0.368)    #Coefficient for position in phi
  tauv ~ dgamma(0.1,0.1)  	    # spread for correlated S's
  tau <- 1/(tauv*tauv)  	#  spread (precision) for correlated S's
#########################
      for (t in 1:T){  #t = primary sample occasion (Weeks)
  	      # Weekly sigma = weekly movement with temp covars
  	      log(sigma[t]) <- alphasig  + betasig2*tempB[t]# + betasig*tempE[t] 
    
    	    sigma2[t] <- sigma[t] * sigma[t]  # Weekly sigma2
	          } #t
##############################  
      for (i in 1:M){ #m individuals, first are period of first tagging   
  	   	      z[i,first[i]]  ~ dbern(1)  # state=Known to be alive at entry into study area
  	      S[i,first[i]] ~ dunif(1,100)  #Possible fish reaches in tagging period
	   
        for(j in 1:nreach) {  #j locations 
	        D2[i,j,first[i]] <- pow(S[i,first[i]]-reach.loc[j], 2)   
         lam[i,j,first[i]] <- lam0*exp(-D2[i,j,first[i]]/(2*sigma2[first[i]])) 
          tmp[i,j,first[i]] <- lam[i,j,first[i]] 
          y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])
               } #j
        
        for (t in (first[i]+1):T) { #t are weeks
            S[i,t] ~ dnorm(S[i,t-1] + betaT2*tempB[t-1]+betaT3*tempE[t-1], tau) 
 
            for(j in 1:nreach) { #j locations
		            D2[i,j,t] <- pow(S[i,t]-reach.loc[j], 2)   
                lam[i,j,t] <- lam0 * exp(-D2[i,j,t]/(2*sigma2[t]))
	              tmp[i,j,t] <- z[i,t]*lam[i,j,t]
	              y[i,j,t] ~ dpois(tmp[i,j,t])
		 		   	        } #j
 	   	
 	   		      #put a logistic regression on phiUp here where it is in reach 
	       	      logit(phi[i,t])<-alphaphi +  betaphiT*tempB[t]#betaphipos*S[i,t]+
   	   	        phi2[i,t] <- max(0.0001,min(phi[i,t],0.9999))#so doesn't go outside 0 and 1 and crash
  	   		      phiUP[i,t] <- z[i,t-1]*phi2[i,t]      
 	              z[i,t] ~ dbern(phiUP[i,t])	# states alive or not
	       	  
	       	        } # t
              } # m
            } #model

", fill = TRUE)
sink()
######################################################
#####  END OF TEXT FILE   ############################
######################################################
#Set up a data input
JAGSdata<-list(y=y, first=first, M=M, T=T, nreach=nreach,tempB=tempB,
               reach.loc=reach.loc,tempE=tempE
) 

z=matrix(NA, M, T)#latent variable
for(i in 1:M){ 
  for(t in first[i]:T){ 
    z[i,t] <-1
  } 
}

#Set initial values
inits =  function() {list(z=z,tauv=runif(1,8,14),alphasig=runif(1,1,2.5),#betaphiT=runif(1,2,5),
                          alphaphi=runif(1,2,5)#betasig=runif(1,0.1,1),#lam0=runif(1,0.08,0.2),
                          #betaF=runif(1,0.5,5),betaT=runif(1,0.5,5)#,betaFL.S=runif(1,-5,5),
                          #betaFL.phi=runif(1,-5,5)
) 
}
# Parameters to follow
parameters <- c("alphaphi","alphasig","betasig",#"betasig2",
                "betaT2","betaphiT",#"betaT3",#"betaphipos",
                "S","phi","z","sigma","tauv","lam0"
)
#Call JAGS
Sys.time()

ZZ<-jags(data=JAGSdata,inits=inits,parameters.to.save=parameters,model.file="ModelCJS.txt",n.thin=1,
        
          #n.chains=3, n.burnin=1500,n.iter=8000,parallel=TRUE)#n.adapt = 10, 
        n.chains=3, n.burnin=100,n.iter=800,parallel=TRUE)#n.adapt = 10,
Sys.time()
###################################################################################

#check convergence in trace plots
MCMCtrace(ZZ, params = c("alphaphi","alphasig","betasig",#"betasig2",
                         "betaT2","betaphiT","lam0",#"betaT3",#"betaphipos",
                         "tauv"))
#get the draws
mco2 <- mcmcOutput(ZZ, params=c("alphaphi","alphasig","betasig",#"betasig2",
                                "betaT2","betaphiT",#"betaT3",#"betaphipos",
                                "S","phi","z","sigma","tauv","lam0"
))
mco2summary <- summary(mco2,CRImass=0.90)
View(summary(mco2["lam0"],CRImass=0.90)) 
# View(summary(mco2[c("alphaphi","alphasig","alphalam",
#                     "betaF", "betaT2","betaphipos","tauv", "gamma")],CRImass=0.90)) 
#                      
###############################################################################
###1. model selection for indicator variables
################################################################################
# mod <- mco2$w
# mod <- paste(mod[,1],mod[,2],mod[,3],mod[,4],mod[,5],mod[,6],
#              mod[,7],mod[,8],mod[,9],mod[,10],mod[,11],
#              mod[,12],mod[,13],mod[,14],mod[,15],
#              sep="")
# table(mod)#creates list of models
# #sort(round(table(mod)/length(mod), 5))#creates AIC table COOL!!!!
# restable <-round(table(mod)/length(mod), 5)
# restable <- as.data.frame(restable)
# weights <-  restable %>%
#   arrange(desc(Freq))
# write.csv(weights,"post_weights1992_.csv", row.names = F)
######################################################################
###make objects for the draws for Smov, z, phi, sigma, lam0
Smov <- mco2$S
zstate <- mco2$z 
sigma <- mco2$sigma
phi <- mco2$phi
lam0 <- mco2$lam0
##########get sigma estimates for manuscript with 90% CIs
(meansig <- mean(sigma,na.rm = TRUE))
(lCLsig <- quantile((sigma),probs= 0.05,na.rm = TRUE)) 
(UCLsig <- quantile((sigma),probs= 0.95,na.rm = TRUE))
############################################################
#get the draws for phi and S into a format we can work with by weeks and reaches
counter <-1
m_predict <- matrix(0,nrow=437, ncol = 7)
for (i in 1:23) {
  for (j in 1:19){
    m_predict[counter,1]<- mean(Smov[,i,j])
    m_predict[counter,2]<-j
    m_predict[counter,3]<-i
    m_predict[counter,4]<- mean(phi[,i,j])
    m_predict[counter,5]<- quantile(phi[,i,j],probs= 0.05,na.rm = TRUE) 
    m_predict[counter,6]<- quantile(phi[,i,j],probs= 0.95,na.rm = TRUE) 
    m_predict[counter,7]<- mean(zstate[,i,j])
    counter <- counter + 1
  }
}

new_frame<-as.data.frame(m_predict)
colnames(new_frame) = c("avg.S","week","fish","avg.phi","LCL.phi","UCL.phi","z")
#########################################################################
#get draws for phi to make weekly plot for phi 
#########################################################################
wk.phi <- new_frame%>%
  filter(avg.phi!="NA")%>%
  group_by(week)%>%
  summarize(weekphi=mean(avg.phi),LCLphi=mean(LCL.phi),UCLphi=mean(UCL.phi))
wk.phi$week2 <- c(seq(1:18))
##get total survival over the whole period
appStot <- prod(wk.phi$weekphi)
mean(wk.phi$weekphi)
mean(wk.phi$LCLphi)
mean(wk.phi$UCLphi)
##with 19-first temps, i.e., t-1 corresponding to the 18 phi estimates
TempBpost <- c( 12.802443, 11.671784, 10.438367,  8.858000,  6.894857,  5.165571,
                3.610000,  4.276857,  3.622832,  2.864495,  3.625040,  3.874573,  4.241382,
                5.860427,  6.685168,  7.544143, 8.440793,  9.117784)
wk.phi <- as.data.frame(wk.phi)
##########################################################################
##!!!!MAKE A PLOT FOR SURVIVAL BY WEEK with credible intervals!
#dev.off()
ggformat<-(
  ggplot(wk.phi, aes(x = week2, y=weekphi,ymin = LCLphi, ymax = UCLphi)) +
    geom_point()+
    geom_pointrange(color="dodgerblue",size=0.5, shape = 20)+
    geom_line(aes(y = TempBpost/15), color="red4")+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "Aug 19", 
                                                                  "3" = "Sept 16", "5" = "Oct 14","7" = "Nov 11","9" = "Dec 9",
                                                                  "11" = "Jan 6", "13" = "Feb 4","15" = "Mar 4","17" = "Apr 1"))+
    #limits = c(1, 18))+
    scale_y_continuous("Mean survival",sec.axis = sec_axis(~ . *15, name = " Temp ºC"),
                       breaks = seq(-5,15, by=5),limits = c(-0.24, 1.00))+
    theme(plot.margin = margin(0,1,0,0, "cm"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=9.5,color="black"))+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.5), angle = 00))+
    labs(title = "", x = "")
)
ggsave(filename="BTPOSTweeklyphiDec19.png", plot=ggformat, device="png",
       height=5, width=6, units="in", dpi=800) 
##############################################################
#make PLOT of FISH MOVEMENTS (by fish)
################################################################
movements <- new_frame%>%#filter out NAs and z<0.99
  filter(z>0.7)
##temperature B with S to plot
tempBpost <- as.data.frame(TempB)
tempBpost$week <- sort(unique(movements$week))
df1 <- left_join(tempBpost,movements)
##make tempE variable to plot
tempEpost <- as.data.frame(TempE)
tempEpost$week <- sort(unique(movements$week))
tempEpost <- as.data.frame(tempEpost)
df <- left_join(tempEpost,df1)
###############################################################################
####this is for movements for POST, 23 fishes and 19 time periods
##############################################################################
df$fish <- as.factor(df$fish)
#########LONG  MOVEMENTS#######################################################
longmov<- df%>%
  filter(fish==4|fish==5|fish==6|fish==9|fish==11|fish==13|fish==14|fish==16|fish==18|fish==19|fish==20)
longmov$fish <- as.factor(longmov$fish)
##################################################################
library(RColorBrewer)
cols <- c("4"="red1","5" = "black", "6"="red3","9" = "darkorchid1", "11" = "orange",
          "13" = "cyan2", "14" = "gold1","16" = "gold4", "18" = "olivedrab2", "19" = "magenta1","20" = "blue1")
ggformat <- (
  ggplot(longmov, aes(x = week, y=99-avg.S, color=fish)) + #96.7-
    geom_point()+
    geom_line()+
    scale_colour_manual(values = cols)+
    geom_line(aes(y = TempB*6.5), color="red",linetype = "dashed")+
    geom_line(aes(y = TempE*6.5), color="red",linetype = "dashed")+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "Aug 19", 
                                                                  "3" = "Sept 16", "5" = "Oct 14","7" = "Nov 11","9" = "Dec 9",
                                                                  "11" = "Jan 6", "13" = "Feb 4","15" = "Mar 4","17" = "Apr 1"))+
    scale_y_continuous("Reach location (RKM)",limits = c(-12, 100),
                       sec.axis = sec_axis(~ . /6.5,name = "Temp ºC",breaks = c(-5,0,5,10,15)))+ 
    theme_bw()+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11,color="black"))+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=8.5,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.3), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.3), angle = 00))+
    theme(legend.position = "none") 
)
ggsave(filename = "BT_MovementsPOSTLONGMOVDec192022.png", plot=ggformat,
       device="png",height=5,width=6,units = "in",dpi=800)
###########################################################################################
midmov<- df%>%
  filter(fish==7|fish==8|fish==12|fish==15|fish==17|fish==21)
midmov$fish <- as.factor(midmov$fish)
##################################################################
cols <- c("7" = "black", "8" = "darkorchid1", "12" = "orange",
          "15" = "cyan2",  "17" = "olivedrab2", "21" = "magenta1")
ggformat <- (
  ggplot(midmov, aes(x = week, y=99-avg.S, color=fish)) + #96.7-
    geom_point()+
    geom_line()+
    scale_colour_manual(values = cols)+
    geom_line(aes(y = TempB*6.5), color="red",linetype = "dashed")+
    geom_line(aes(y = TempE*6.5), color="red",linetype = "dashed")+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "Aug 19", 
                                                                  "3" = "Sept 16", "5" = "Oct 14","7" = "Nov 11","9" = "Dec 9",
                                                                  "11" = "Jan 6", "13" = "Feb 4","15" = "Mar 4","17" = "Apr 1"))+
    scale_y_continuous("Reach location (RKM)",limits = c(-12, 100),
                       sec.axis = sec_axis(~ . /6.5,name = "Temp ºC",breaks = c(-5,0,5,10,15)))+ 
    theme_bw()+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11,color="black"))+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=8.5,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.3), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.3), angle = 00))+
    theme(legend.position = "none") 
)
ggsave(filename = "BT_MovementsPOSTMIDMOVDec192022.png", plot=ggformat,
       device="png",height=5,width=6,units = "in",dpi=800)
#########################################################################
#############################################################################################
abmov<- df%>%
  filter(fish==1|fish==2|fish==3|fish==22|fish==23|fish==10)
abmov$fish <- as.factor(abmov$fish)
##################################################################
cols <- c("1"="red1","2" = "black", "3" = "darkorchid1", "22" = "orange",
          "23" = "cyan2", "10" = "gold1")
ggformat <- (
  ggplot(abmov, aes(x = week, y=99-avg.S, color=fish)) + #96.7-
    geom_point()+
    geom_line()+
    scale_colour_manual(values = cols)+
    geom_line(aes(y = TempB*6.5), color="red",linetype = "dashed")+
    geom_line(aes(y = TempE*6.5), color="red",linetype = "dashed")+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "Aug 19", 
                                                                  "3" = "Sept 16", "5" = "Oct 14","7" = "Nov 11","9" = "Dec 9",
                                                                  "11" = "Jan 6", "13" = "Feb 4","15" = "Mar 4","17" = "Apr 1"))+
    scale_y_continuous("Reach location (RKM)",limits = c(-12, 100),
                       sec.axis = sec_axis(~ . /6.5,name = "Temp ºC",breaks = c(-5,0,5,10,15)))+ 
    theme_bw()+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11,color="black"))+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=8.5,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.3), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.3), angle = 00))+
    theme(legend.position = "none") 
)
ggsave(filename = "BT_MovementsPOSTABMOVDec192022.png", plot=ggformat,
       device="png",height=5,width=6,units = "in",dpi=800)

# 
# ggformat <- (
#   ggplot(df, aes(x = week, y=97.7629-avg.S, color=fish)) +#subtract from the max distance
#     geom_point(size = 1)+
#     #geom_line(aes(y = TempE*6.5), color="red4")+
#     #geom_line(aes(y = TempB*6.5), color="red4")+
#     scale_colour_manual(values = cols)+
#     #scale_color_viridis_d()+
#     scale_x_continuous(breaks = c(3,6,9,12,15,18),labels=c("3" = "Sept 16",
#             "6" = "Oct 28","9" = "Dec 9","12" = "Jan 21", "15" = "Mar 4","18" = "Apr 15"))+
#     scale_y_continuous(name= "River location (RKM)",breaks = c(25,50,75,100))+
#                      #  sec.axis = sec_axis(~ ., name = "Max water temperature ºC"))+
#     labs(title = "", y = "River location (RKM)", x = "Time period (bi-week)")+
#      theme_bw()+
#     ylim(0,100)+
#     #theme(axis.text.x = element_text(angle = 0, hjust = 0, size=8,color="black"))+
#     theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11,color="black"))+
#     theme(axis.title.y = element_text(size = rel(1.3), angle = 90))+
#     theme(axis.title.x = element_text(size = rel(1.3), angle = 00))+
#     theme(legend.position = "none")
# )
# ggsave(filename = "BT_MovementsPOSTMay11.png", plot=ggformat,
#        device="png",height=5,width=5,units = "in",dpi=800)
# #############################################################################################
