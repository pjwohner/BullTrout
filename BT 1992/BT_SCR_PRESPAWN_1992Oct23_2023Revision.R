# Code provided by B. Gardner and modifed by J. Hightower and J. Raabe for spatially 
#explicit capture-recapture model for Little River 
# Similar to examples in Gardner et al. 2010, Kery et al. 2010, Sollmann et al. 2011
library(reshape);library(tidyverse);library(nnet)
library(MCMCvis);library(rjags);library(jagsUI);library(mcmcOutput)
#setwd("C:/newlaptop/FISH/RapidRiverBullTrout/R_scripts/SCR_BT/")
setwd("~/patti1/Bull Trout")
AS10pre <- read.table("BT_SCR_detect1992_1Knov132022.txt",header=T)
##add known fish weeks to be able to plot later on (becasue same distance between weeks)
TagID <- c(54,54,74,74,174,195)
Week <- c(11,13,6,8,8,4)
RKM <- c( 28044.0, 28044.0,10593.9,10593.9,28044.0,10593.9)
xASpre <- data.frame(TagID, Week,RKM ) 
AS10pre1992 <- rbind(AS10pre,xASpre ) 
AS10 <- AS10pre1992%>%
  arrange(TagID,Week,RKM)%>%
  #change from meters to km  
  mutate(RKM =RKM/1000)%>%
  filter(RKM!=40.7515)
############################################################################
#sort(unique(AS10$RKM))
melted.rkm <- melt(AS10, id=c("TagID","RKM")) 
tagid.week <- cast(melted.rkm, TagID ~ RKM ~ value, fill=0, length)
y <- tagid.week
#y[36,2,2]
#sort(unique(AS10$RKM))
M <- length(table(AS10$TagID))    # Number of individuals
T <- length(table(AS10$Week))     # Number of periods (weeks) 
#capture weir at 4000 km, need this to match in data file
reach.loc <- c(4.0000,  5.3178,  6.3187,  7.3269,  8.3751,  9.4551, 10.5939, 11.6243,
               12.6317, 13.6615, 14.6652, 15.6694,16.8046, 17.8232, 18.8351, 19.8584,
               20.8886, 21.9082, 22.9129, 23.9697, 24.9739, 25.9780, 27.0153, 28.0440,
               29.1170, 30.1574, 31.1865, 32.2015, 32.6660, 33.2486, 33.6660, 34.3071,
               34.6660, 35.3173, 35.6660, 36.3231, 36.6660, 37.6660, 38.9350, 40.0183)  

nreach <- length(reach.loc)  # weir + reaches in BT case (not antenna)

# Input and format data matrix:
# y = detection history (individual (i) x reach (j) x week(t)); first = river 
#immigration week, last = river emigration week, (these are counts).
##first week fish came into study
first.dat<-read.csv("1992firstcapnov13.csv",stringsAsFactors = F) #all in first week for 1992
#firstc<-read.csv("firstcap.csv",stringsAsFactors = F)
first <- first.dat %>%
  select(firstcap)
first <- as.numeric(first$firstcap)
##make a last to censor live individuals to when they were last seen alive
#if not seen again (know they are alive from later detections)
# last <- c(18,17,16,16,18,18,16,18,2,16,16,12,18,17,18,17,17,16,18,17,12,18,18,16,14,
#           17,18,18,12,16)
#############################################################################
##fish size standardized from prep file
fishsize.dat <-read.csv("fishsize1992dec22022.csv",stringsAsFactors = F) 
fishsize1992 <- fishsize.dat %>%
  select(FL)
sizeFL <- as.numeric(fishsize1992$FL)
################################################################################
#flowA92 is  average flow 1992 (weeks 1-18, not including week 0)
FlowA92 <- c(5.416608, 4.203029, 3.660964, 3.519380, 3.211940, 3.559832, 3.126989,
           2.486927, 2.244582, 1.980966, 1.806210, 1.694561, 1.675548, 1.614060,
           1.513333, 1.476521, 1.425551, 1.450227)
###flowB is delta maximum - minimum within the week 1992
FlowB92 <- c(1.047723345, 0.906139109, 0.141584236, 0.311485319, 0.538020097, 0.509703249, 
             0.509703249, 0.453069555, 0.155742659, 0.223703092, 0.104772335, 0.070792118,
             0.079287172, 0.158574344, 0.008495054, 0.008495054, 0.065128748, 0.110435704)
###flowC is delta among weekly maximum - minimum  week 1992
 FlowC92 <- c(-0.47329587, -1.21357916, -0.54206536, -0.14158424, -0.30744005,  
              0.34789269, -0.43284324, -0.64006188, -0.24234502, -0.26361636,
              -0.17475540, -0.11164928, -0.01901274, -0.06148801, -0.10072707,
              -0.03681190, -0.05097033,  0.02467611)
##################################################################################
###Temp variables: see paper for context
#TempA is A) weekly maximum hourly water temperature avg
 TempA92 <- c(13.5, 14.1, 14.7, 16.6, 17.0, 13.9, 14.8, 16.2, 16.9, 17.8,
              16.7, 17.1, 17.0, 14.4, 14.5, 13.2, 13.1, 13.4)
##Temp B is weekly average daily maximum water temperature
 TempB92 <- c(12.20000, 12.68571, 12.11429, 14.42857, 15.67143, 13.10000, 13.55714,
              15.46250, 15.70000, 16.98571,16.07143, 16.51429, 14.52857, 13.27143,
              12.60000, 12.15714, 11.91429, 11.95714)
###Temp C is weekly three-day moving average of the maximum daily temperature avg
TempC92 <- c(12.10667, 12.32667, 12.27333, 14.70000, 15.72000, 13.16000, 13.46667,
             15.50000, 15.64167, 17.01333, 15.98000, 16.51333, 14.72000, 13.42667,
             12.57333, 12.38000, 11.87333, 11.84000)
###Temp D is weekly average hourly temperature avg
TempD92 <- c(10.439286, 10.704167, 10.707143, 11.947024, 13.755952, 11.786310,
             11.945238, 13.265625, 13.338194, 14.333333, 13.761310, 14.494048,
             13.091667, 11.194643, 11.311310, 10.657143, 10.479762, 10.649405)
###Temp E is delta current - previous week average
TempE92 <- c(0.81607143,  0.26488095,  0.00297619,  1.23988095,  1.80892857,
             -1.96964286,  0.15892857,1.32038690,  0.07256944,  0.99513889,
             -0.57202381,  0.73273810, -1.40238095, -1.89702381,  0.11666667,
              -0.65416667, -0.17738095, 0.16964286)
##scale temp and flow vars
flowA=(FlowA92-mean(FlowA92))/(sd(FlowA92))#standardize flow variable
flowB=(FlowB92-mean(FlowB92))/(sd(FlowB92))#standardize flow variable
flowC=(FlowC92-mean(FlowC92))/(sd(FlowC92))#standardize flow variable
tempA=(TempA92-mean(TempA92))/(sd(TempA92))#standardize flow variable
tempB=(TempB92-mean(TempB92))/(sd(TempB92))#standardize flow variable
tempC=(TempC92-mean(TempC92))/(sd(TempC92))#standardize flow variable
tempD=(TempD92-mean(TempD92))/(sd(TempD92))#standardize flow variable
tempE=(TempE92-mean(TempE92))/(sd(TempE92))#standardize flow variable
##### test STANDARDIZED flow and temp correlations###########################################
#FlowA FlowB and FlowC uncorrelated with all Temp vars
#but flowA and flowB vars correlated 
#tempA,B,C,D correlated with eachother, TempE not correlated with anything
# flowmets <- matrix(0,nrow=18, ncol = 8)
# for (i in 1:18) {
#   flowmets[i,1]<- flowA[i]
#   flowmets[i,2]<- flowB[i]
#   flowmets[i,3]<- flowC[i]
#   flowmets[i,4]<- tempA[i]
#   flowmets[i,5]<- tempB[i]
#   flowmets[i,6]<- tempC[i]
#   flowmets[i,7]<- tempD[i]
#   flowmets[i,8]<- tempE[i]
# }
# colnames(flowmets) <- c("flowA","flowB","flowC","tempA","tempB","tempC","tempD","tempE")
# flowmets<-as.data.frame(flowmets)
# to_test<-c("flowA","flowB","flowC","tempA","tempB","tempC","tempD","tempE")
# cor(flowmets[,to_test])
###bring in covariates for reaches we want joined
JimCovRRa <- read.csv("JIMSCOVSnov3.csv",header=T)
#############################
##make plot of shan and reach distance
################################
# ggformat<-(
#   ggplot(JimCovRR,aes(x=CumSum/1000+59.0969,y=Shan))+
#     geom_point(size=2, shape = 20) +
#     geom_smooth(method = "glm",se = TRUE,level = 0.90,color = "black")+
#     theme_bw()+
#     labs(title = "", x = "Rapid River reach (RKM)", y = "Shan-Wien index")+
#     theme(axis.text.x = element_text(angle = 0, hjust = 0, size=12,color="black"))+
#     theme(axis.text.y = element_text(angle = 0, hjust = 0, size=12,color="black"))+
#     theme(axis.title.y = element_text(size = rel(1.6), angle = 90))+
#     theme(axis.title.x = element_text(size = rel(1.6), angle = 00))+
#     theme(plot.margin = margin(0,1,0,0, "cm"))+
#     ylim(1,2.5)+
#     xlim(60,100)
#    )
# ggsave(filename="BTPREshanindbyRKMdec18.png", plot=ggformat, device="png",
#        height=5, width=5, units="in", dpi=800)

JimCovRRa$CumSum <- as.numeric(JimCovRRa$CumSum)
##make average length variables
JimCovRR <- JimCovRRa %>%
  #filter out reaches we don't use
  filter(newRea1K.y !=1&newRea1K.y !=2&newRea1K.y!=3 & newRea1K.y!=41 & newRea1K.y!=43& newRea1K.y!=46)%>%
  arrange(CumSum)%>%
  mutate(AVG.LEN=(RITE_LEN+LEFT_LEN)/2)%>%#make average reach length
  #mutate(TotAreaha=TotArea/10000)%>%
  mutate(TOT_VGST=(LEN_VGST_L+LEN_VGST_R)/AVG.LEN)%>%##VGSTproportion
  mutate(LEN_UCUT=(LEN_UCUT_L+LEN_UCUT_R)/AVG.LEN)%>% #ratio Undercut bank
  mutate(dens.LWD_TOT=(LWD_AGREGS+LWD_WADS)/AVG.LEN)%>% #total LWD density
  # mutate(dens.LWD_AGREGS=LWD_AGREGS/TotArea)%>%#make a LWD aggregs density
  # mutate(dens.LWD_WADS=LWD_WADS/TotArea)%>% #make a LWD aggregs density
  mutate(HG.prop=(HG.area)/TotArea)%>%
  mutate(RN.prop=(RN.area)/TotArea)%>%
  mutate(CS.prop=(CS.area)/TotArea)%>%
  mutate(LG.prop=(LG.area)/TotArea)%>%
  mutate(pools=(PL.area+PW.area+SP.area+UP.area+DP.area)/TotArea)%>%
  select(-LEN_VGST_R,-LEN_VGST_L,-LEN_UCUT_L,-LEN_UCUT_R,-LEFT_LEN,-RITE_LEN,-LWD_AGREGS,
         -LWD_WADS,-newRea1K.x,-SegmentLen.x,-newRea1K.y,-SegmentLen.y,-NR1K)
##do PCA on the habitat area variables to come up with a new variable
pcajimcov <-JimCovRR%>% 
  select(TOT_VGST,LEN_UCUT,dens.LWD_TOT)
#pcajimcov<-JimCovRR[,2:25]
modelVS<-prcomp(pcajimcov,scale=TRUE)
modelVS$rotation[,1:3]
###to get the loadings for each data point for PC1 to 3
PCVS<-(modelVS$x[,1:2])
PCVS <- -PCVS
colnames(PCVS) <- c("PC1", "PC2")
#combine the PCs to the last dataset with coverpc
JimCovRR<-cbind(JimCovRR,PCVS)
##standardize covariates (but they are all small values anyway so I doin't do it)
JimCovRRst <-JimCovRR%>%
  mutate_at(c("Shan", "RN.prop", "CS.prop","LG.prop","LEN_UCUT","dens.LWD_TOT",
              "pools","PC1","TOT_VGST"
              ), ~(scale(.) %>% as.vector))
####check for correlation (nothing correlated!!)
#names(JimCovRR)
#  x <- JimCovRR[2:29]
#  y <- JimCovRR[2:29]
#  cor(x,y)
############################################################
##these are the covars from model selction
Cov1 <- JimCovRRst$Shan
#Cov2 <- JimCovRRst$HG.prop
#Cov3 <- JimCovRRst$RN.prop
#Cov4 <- JimCovRRst$CS.prop
#Cov5 <- JimCovRRst$LG.prop
Cov2 <- JimCovRRst$LEN_UCUT
Cov3 <- JimCovRRst$dens.LWD_TOT
#Cov8<- JimCovRRst$MeanDepth
#Cov9 <- JimCovRRst$TOT_VGST
#Cov4 <- JimCovRRst$PC1
#Cov5 <- JimCovRRst$pools#pools habitat combined##################
###################################################################
#####  WRITE TEXT FILE WITH JAGS SPECIFICATION   ###########
sink("ModelCJS.txt")
cat("
model { #model
# Priors

alphasig ~ dnorm(0,0.368) 		# Intercept in sigma estimate  
alphaphi ~ dnorm(0,0.368)     #intercept for phi (needs to be positive)
alphalam ~ dnorm(0,0.368)    #intercept for lam0
alphalam1 ~ dnorm(0,0.368)    #intercept for lam0 first encounters
  #betaphiFlowA  ~ dnorm(0,0.368) #beta for phi for flowA
  #betaphiFlowC  ~ dnorm(0,0.368) #beta for phi for flowC
betaF  ~ dnorm(0,0.368)      #Coefficient for flowA in S estimate 
 #betaF2 ~ dnorm(0,0.368)      #Coefficient for flowC in S estimate
 #betasig  ~ dnorm(0,0.368)     #Coefficient for flow C in sigma estimate 
 #betasig2 ~ dnorm(0,0.368)     #Coefficient for temp in sigma estimate
betaphipos~ dnorm(0,0.368)    #Coefficient for position in phi
 #betaphiTemp ~ dnorm(0,0.368)  #beta for phi for tempB
betaT2 ~ dnorm(0,0.368)       #Coefficient for tempB in S estimate 
 #betaFL.phi ~ dnorm(0,0.368)   #beta for phi for size
 #betaFL.S ~ dnorm(0,0.368)     #beta for S for size

tauv  ~ dgamma(0.1,0.1) 	    #spread for correlated S's (needs to be positive)
tau <- 1/(tauv*tauv)  	      #spread (precision) for correlated S's

#priors reach encounter betas
            for(zz in 1:3){
              gamma[zz]~dnorm(0,0.368)
            }
 # #priors flow and temp betas
            # for(xx in 1:17){
            #   tfgamma[xx]~dnorm(0,0.368)
            # }
            
# # #indicator variable prior
         # for(yy in 1:30){
         #  w[yy] ~ dbern(0.5)
         # }

##sigma is how far fish can be available to move into reach from other reaches
    for (t in 1:T){  # t = primary sample occasion (Weeks)
  	    log(sigma[t]) <- alphasig 
    	
    	       # log(sigma[t]) <- alphasig + w[10]*tfgamma[1]*flowA[t] + w[11]*tfgamma[2]*flowC[t]+
    	       #       w[12]*tfgamma[3]*tempE[t] + w[13]*tfgamma[4]*tempB[t]
    	  sigma2[t] <- sigma[t] * sigma[t]  # Weekly sigma2
	          } #t
	        
      for (i in 1:M){ # m individuals, first are period of first tagging   
  	    z[i,first[i]]  ~ dbern(1)  # state = Known to be alive at entry into study area
  	    S[i,first[i]]  ~ dunif(4,40)  # this is where fish enter study was 4
     
        for(j in 1:nreach) {  #j locations 
	        D2[i,j,first[i]] <- pow(S[i,first[i]]-reach.loc[j], 2)   
         log(lam0[i,j,first[i]]) <- alphalam1
         lam[i,j,first[i]] <- lam0[i,j,first[i]]*exp(-D2[i,j,first[i]]/
                                       (2*sigma2[first[i]])) 
         tmp[i,j,first[i]] <- lam[i,j,first[i]] 
         y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])
               } #j       
      for (t in (first[i]+1):T) {
        #for (t in (first[i]+1):T) { #t are weeks we have 
          S[i,t] ~ dnorm(S[i,t-1] + betaF*flowA[t-1]+ betaT2*tempB[t-1],tau)
          
             #S[i,t] ~ dnorm(S[i,t-1] + w[14]*tfgamma[5]*flowA[t] + w[15]*tfgamma[6]*flowC[t]+
     	         #w[16]*tfgamma[7]*tempE[t] + w[17]*tfgamma[8]*tempB[t]+ w[29]*betaFL.S*sizeFL[i] ,tau)
     
  #activity center is a funciton of where it was last time step and flow...with 
  #some error/variability (tau)...you can truncate this to make sure they stay in 
  #places it is possible to go (i.e., if they cant get through a barrier on the river). 
  #You could also add an intercept if you think there is some constant change in the 
  #activity center across time periods.
         		                          
        for(j in 1:nreach) { #j locations
		        D2[i,j,t] <- pow(S[i,t]-reach.loc[j], 2)   
            log(lam0[i,j,t]) <- alphalam + gamma[1]*Cov1[j] + gamma[2]*Cov2[j]+ gamma[3]*Cov3[j]
            
          #gamma[4]*Cov4[j] #+ gamma[5]*Cov5[j] #+ gamma[6]*flowA[t] + gamma[7]*flowC[t] 
          # log(lam0[i,j,t]) <- alphalam +  w[1]*gamma[1]*Cov1[j] + w[2]*gamma[2]*Cov2[j] +
          #   w[3]*gamma[3]*Cov3[j] + w[4]*gamma[4]*Cov4[j] + w[5]*gamma[5]*Cov5[j]+
          #   w[6]*gamma[6]*Cov6[j] + w[7]* gamma[7]*Cov7[j] +w[8]*gamma[8]*Cov8[j] +
          #   w[9]*gamma[9]*Cov9[j] + w[27]*gamma[10]*Cov10[j]+ w[28]*gamma[11]*Cov11[j]+
          #  w[23]*tfgamma[14]*flowA[t] + w[24]*tfgamma[15]*flowC[t] + 
          #  w[25]*tfgamma[16]*tempE[t] + w[26]*tfgamma[17]*tempB[t]
      
            lam[i,j,t] <- lam0[i,j,t] * exp(-D2[i,j,t]/(2*sigma2[t]))#pdf
	          tmp[i,j,t] <- z[i,t]*lam[i,j,t]
	          y[i,j,t] ~ dpois(tmp[i,j,t])
	                    } #j
  
          ###put a logistic regression on phi
          logit(phi[i,t])<- alphaphi + betaphipos*S[i,t] 
                      #+ betaphiFlowA*flowA[t]+betaphiTemp*tempB[t] 
        
                #+ betaFL.phi*sizeFL[i]#betaphiFlowC*flowC[t]
                # logit(phi[i,t])<- alphaphi + w[18]*tfgamma[9]*S[i,t] + w[19]*tfgamma[10]*flowA[t] +
                # w[20]*tfgamma[11]*flowC[t] + w[21]*tfgamma[12]*tempE[t] + w[22]*tfgamma[13]*tempB[t]+
                #  w[30]*betaFL.phi*sizeFL[i]
        
            phi2[i,t] <- max(0.0001,min(phi[i,t],0.9999))#prevents outside 0 and 1
            phiUP[i,t] <- z[i,t-1]*phi2[i,t]
            z[i,t] ~ dbern(phiUP[i,t])	# states alive or not
               	   } # t
  	           } # m
  	           
   ### Derived quantities #############################################
      #extrapolate survival over the periods
      #mean1 <- mean(phi[,1])
      mean2 <- mean(phi[,2])
      mean3 <- mean(phi[,3])
      mean4 <- mean(phi[,4])
      mean5 <- mean(phi[,5])
      mean6 <- mean(phi[,6])
      mean7 <- mean(phi[,7])
      mean8 <- mean(phi[,8])
      mean9 <- mean(phi[,9])
      mean10 <- mean(phi[,10])
      mean11 <- mean(phi[,11])
      mean12 <- mean(phi[,12])
      mean13 <- mean(phi[,13])
      mean14 <- mean(phi[,14])
      mean15 <- mean(phi[,15])
      mean16 <- mean(phi[,16])
      mean17 <- mean(phi[,17])
      mean18 <- mean(phi[,18])
     TOTALS <- mean2*mean3*mean4*mean5*mean6*mean7*mean8*mean9*mean10*mean11*mean12*mean13*
mean14*mean15*mean16*mean17*mean18
       
        # for (i in 1:M){#individual
        #            for (t in 2:T) {
        #            meanS[t] <- mean(phi[i,t])
        #            #TotalS[i,t] = TotalS[i,t]*phi[i,t]
        #                 }
        #           }
        
  	           
        } # end model
", fill = TRUE)
sink()
######################################################
#####  END OF TEXT FILE   ############################
######################################################
#Set up a data input
JAGSdata<-list(y=y,first=first,M=M,T=T,nreach=nreach,tempB=tempB,#tempE=tempE,
              reach.loc=reach.loc,flowA=flowA,#flowC=flowC,#sizeFL=sizeFL,
              Cov1=Cov1,Cov2=Cov2,Cov3=Cov3#,Cov4=Cov4#,Cov5=Cov5,Cov6=Cov6,
              #Cov7=Cov7#,Cov8=Cov8,Cov9=Cov9,Cov10=Cov10,Cov11=Cov11
               ) 

z=matrix(NA, M, T)#latent variable
for(i in 1:M){
  for(t in first[i]:T){
    z[i,t] <-1
  }
}

#Set initial values
inits =  function() {list(z=z#,#alphaphi=runif(1,0,25),alphalam1=runif(1,-15,3),
           # tauv=runif(1,0,20),alphasig=runif(1,-3,3),alphalam2=runif(1,-3.5,3.5),
           # betasig=runif(1,-1,1),betaF= runif(1,-2,2),betaFL.phi=runif(1,-1,1),
           # betaFL.S=runif(1,-2,2),betaphiFlow =runif(1,0,25),betaphiTemp =runif(1,0,22),
           # betaT2=runif(1,-2,2),betasig2=runif(1,-3,3)
                          )#,gamma=runif(11,-0.4,0.55),alphaS=runif(1,14,20
                            }
#Parameters to follow
parameters <- c("alphaphi","alphasig","alphalam","alphalam1",
                "betaF", "betaT2","betaphipos","TOTALS",
                "tauv", "gamma","S","phi","z","lam0","sigma"
              )
#cant get sigma in the same output as the rest because too large but do not need for plots

#Call JAGS
#Sys.time()
 # ZZ<-jags(data=JAGSdata,inits=inits,parameters.to.save=parameters,model.file="ModelCJS.txt",n.thin=2,
 #          n.chains=2, n.burnin=50,n.iter=100,parallel=TRUE)#n.adapt = 10,

ZZ<-jags(data=JAGSdata,inits=inits,parameters.to.save=parameters,model.file="ModelCJS.txt",n.thin=2,
           n.chains=3, n.adapt = 2000,n.burnin=6000,n.iter=20000,parallel=TRUE)#
#Sys.time()
#################################################################################
#check convergence in trace plots (must delete old one first)
MCMCtrace(ZZ, params = c("alphaphi","alphasig","alphalam","alphalam1",
                "betaF", "betaT2","betaphipos","tauv", "gamma"
                          ))
                         
##get deviance and parameter estimes with MCE
mco2 <- mcmcOutput(ZZ, params=c("TOTALS","alphaphi","alphasig","alphalam","alphalam1",
                         "betaF", "betaT2","betaphipos","tauv", "gamma",
                         "S","phi","z","lam0","sigma",
                          ))
mco2summary <- summary(mco2,CRImass=0.90)

# View(summary(mco2[c("alphaphi","alphasig","alphalam",
#                     "betaF", "betaT2","betaphipos","tauv", "gamma")],CRImass=0.90))
###############################################################################
###1. model selection for indicator variables
################################################################################
#this example has 11 covariate indicator variables w
# mod <- mco2$w
# mod <- paste(mod[,1],mod[,2],mod[,3],mod[,4],mod[,5],mod[,6],mod[,7],mod[,8],
#              mod[,9],mod[,10],mod[,11], mod[,12],mod[,13],mod[,14],mod[,15],
#              mod[,16], mod[,17],mod[,18], mod[,19],mod[,20],mod[,21],mod[,22],
#              mod[,23],mod[,24],mod[,25],mod[,26],mod[,27],mod[,28],
#             mod[,29],mod[,30],
#              sep="")
# table(mod)#creates list of models
# #sort(round(table(mod)/length(mod), 5))#creates AIC table COOL!!!!
# restable <-round(table(mod)/length(mod), 5)
# restable <- as.data.frame(restable)
# weights <-  restable %>%
#   arrange(desc(Freq))
# write.csv(weights,"post_weights1992_.csv", row.names = F)
######################################################################
###create the draws for S, z, sigma, lam0, and phi to use in plots
######################################################################
#mco2$deviance
Smov <- mco2$S
zstate <- mco2$z 
phi <- mco2$phi
sigma <- mco2$sigma
phi <- mco2$phi
lam0 <- mco2$lam0
##########get sigma estimates for manuscript with 90% CIs
(meansig <- mean(sigma,na.rm = TRUE))
(lCLsig <- quantile((sigma),probs= 0.05,na.rm = TRUE)) 
(UCLsig <- quantile((sigma),probs= 0.95,na.rm = TRUE))
############################################################
# #to backtransform negative estimates into odds ratio(did not do for table)
var1<-exp(-0.358)
##to get estimate back to be over 1
(p3<-1/var1)
##positvie estimates
exp(0.265)
#diagPlot(mco2)      # diagnostic plots
#plot(mco2)          # plot posterior distributions
#summary(mco2)       # display a summary of parameters in the Console
#get the draws into a format we can work with
counter <-1
m_predict <- matrix(0,nrow=540, ncol = 7)
for (i in 1:30) {
  for (j in 1:18){
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
###########################################################################
wk.phi <- new_frame%>%
  filter(avg.phi!="NA")%>%
  group_by(week)%>%
  summarize(weekphi=mean(avg.phi),LCLphi=mean(LCL.phi),UCLphi=mean(UCL.phi))
#get total phi for the entire time period for the MS
(appStot <- prod(wk.phi$weekphi))
mean(wk.phi$weekphi)
mean(wk.phi$LCLphi)
mean(wk.phi$UCLphi)
################################################################################
##!!!!MAKE A PLOT FOR SURVIVAL BY WEEK with credible intervals!
ggformat<-(
ggplot(wk.phi, aes(x = week, y=weekphi,ymin = LCLphi, ymax = UCLphi)) +
  geom_point()+
  geom_pointrange(color="dodgerblue",size=0.5, shape = 20)+
  scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "May 26",
            "3" = "Jun 9","5" = "Jun 23","7" = "Jul 7","9" = "Jul 22",
            "11" = "Aug 4", "13" = "Aug 18","15" = "Sept 1","17" = "Sept 15"),
          limits = c(1, 18))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=9.5,color="black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 10, size=11,color="black"))+
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.5), angle = 00))+
  #theme(plot.margin = margin(0,1.2,1,1, "cm"))+
  labs(title = "", y = "Mean survival", x = "")+
  ylim(0,1)
 )
 ggsave(filename="BTPRE1992weeklyphiDec19.png", plot=ggformat, device="png",
        height=5, width=5, units="in", dpi=800) 
 #dev.off()
#######################################################################
 ##make a plot for survival with S[i,t]
 ###########################################################################################
 #  ggformat<-(
 #    ggplot(new_frame,aes(x=avg.S+59.0969,y=avg.phi, color=fish))+
 #     geom_point(size=2, shape = 20) +
 #     theme_bw()+
 #    scale_color_viridis_b()+
 #      labs(title = "", x = "Position in river (RKM)", y = "Mean survival (??)")+
 #     theme(axis.text.x = element_text(angle = 0, hjust = 0, size=12,color="black"))+
 #     theme(axis.text.y = element_text(angle = 0, hjust = 0, size=12,color="black"))+
 #     theme(axis.title.y = element_text(size = rel(1.6), angle = 90))+
 #     theme(axis.title.x = element_text(size = rel(1.6), angle = 00))+
 #     theme(plot.margin = margin(0,1,0,0, "cm"))+
 #     xlim(50,100)+
 #      theme(legend.position = "none") 
 # )
 # ggsave(filename="BTPRE1992SITdec4.png", plot=ggformat, device="png",
 #        height=5, width=5, units="in", dpi=800)
#############################################################
##make df for movement for fish that are very likely alive
################################################################
#filter out NAs and z<0.7 (not likely alive)
movements <- new_frame%>%
  filter(z!="NA")%>%
  filter(z>0.70)
flow <- as.data.frame(FlowA92)
flow$week <- sort(unique(movements$week))
dfpre <- left_join(flow,movements)
##temperature
temp <- as.data.frame(TempB92)
temp$week <- sort(unique(movements$week))
df <- left_join(temp,dfpre)
df <- as.data.frame(df)
####################################################################################
###all 1992 fish MOVEMENTS (S) plotted together
###################################################################################
df$fish <- as.factor(df$fish)
ggformat <- (
  ggplot(df, aes(x = week, y=avg.S+56, color=fish)) +
  #ggplot(df, aes(x = week, y=avg.S+59.0969, color=fish)) +
    geom_point(size = 1)+
    geom_line(aes(y = FlowA92*6), color="blue4")+
    geom_line(aes(y = TempB92*6), color="red4")+
    scale_color_viridis_d()+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "May 26",
                "3" = "Jun 9","5" = "Jun 23","7" = "Jul 7","9" = "Jul 22",
                 "11" = "Aug 4", "13" = "Aug 18","15" = "Sept 1","17" = "Sept 15"),
                  limits = c(1, 18))+
   scale_y_continuous("Reach location (RKM)",limits = c(0, 102),
          sec.axis = sec_axis(~ . /6, name = "Flow (cms) + Temp ?C",breaks = c(0,5,10,15))
                                    )+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=8.5,color="black"))+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=10,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.2), angle = 00))+
    labs(title = "", x = "")+
    theme(plot.margin = margin(0,1,0,0, "cm"))+
    theme(legend.position = "none") 
)
ggsave(filename = "BT_Movements1992dec19.png", plot=ggformat,
       device="png",height=4,width=5,units = "in",dpi=800)
############################################################################
########   make plot lam0 by reach
############################################################################
   counter <-1
  m_covar <- matrix(0,nrow=21600, ncol = 8)
  for (i in 1:30) {
    for (j in 1:40){
    for (k in 1:18){
      m_covar[counter,1]<- mean(lam0[,i,j,k])#mean lam0
      m_covar[counter,2]<-k
      m_covar[counter,3]<-j
      m_covar[counter,4]<-i
      m_covar[counter,5]<- mean(exp(lam0[,i,j,k]))
      m_covar[counter,6]<- y[i,j,k]
      m_covar[counter,7]<- quantile((lam0[,i,j,k]),probs= 0.05,na.rm = TRUE) 
      m_covar[counter,8]<- quantile((lam0[,i,j,k]),probs= 0.95,na.rm = TRUE)
      counter <- counter + 1
    }
  }
  }
  m_covar<-as.data.frame(m_covar)
  colnames(m_covar) = c("pred.lam0","week","reach","fish","explam0","datay","LCIpredlam0","UCIpredlam0")
###############################################
  ##then do group by reach to get means of lam0
  lam0reach <- m_covar%>%
    filter(pred.lam0 != "NA")%>%
    filter(reach != 0)%>%
    group_by(reach)%>%
    summarize(lam0reach=mean(pred.lam0),LCI=mean(LCIpredlam0),UCI=mean(UCIpredlam0))
  lam0reach$dist <- round(reach.loc+56,0)
  lam0reach<-as.data.frame(lam0reach)
  mean(lam0reach$lam0reach)
  mean(lam0reach$LCI)
  mean(lam0reach$UCI)
#####################################################################
  ##!!!!MAKE A PLOT FOR lam0 BY REACH with credible intervals!
  ggformat<-(
    ggplot(lam0reach, aes(x = reach, y=lam0reach, ymin = LCI, ymax = UCI)) +
      geom_point()+
      geom_pointrange(color="dodgerblue",size=0.5, shape = 20)+
      theme_bw()+
      scale_x_continuous(breaks = c(1,5,9,13,17,21,25,29,33,37,40),
                labels=c("1" = "63","5" = "67", "9" = "72","13" = "76",
                "17" = "80","21" = "84", "25" = "88","29" = "92",
                "33" = "94","37" = "96","40" = "99"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0, size=10,color="black"))+
      theme(axis.text.y = element_text(angle = 0, hjust = 0, size=12,color="black"))+
      theme(axis.title.y = element_text(size = rel(1.4), angle = 90))+
      theme(axis.title.x = element_text(size = rel(1.4), angle = 00))+
      ylim(0,1)+
      labs(title = "", y = "Mean encounter rate", x = "Rapid River Reach (RKM)")
  )
  ggsave(filename="BTPRE1992reachlam0Dec19.png", plot=ggformat, device="png",
         height=5, width=5, units="in", dpi=800) 
