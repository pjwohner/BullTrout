# Code provided by B. Gardner and modifed by J. Hightower and J. Raabe for spatially 
#explicit capture-recapture model for Little River 
# Similar to examples in Gardner et al. 2010, Kery et al. 2010, Sollmann et al. 2011
library(reshape);library(tidyverse);library(nnet)
library(MCMCvis);library(rjags);library(jagsUI);library(mcmcOutput)
setwd("~/patti1/Bull Trout 1993")
#setwd("C:/newlaptop/FISH/RapidRiverBullTrout/R_scripts/SCR_BT/")
AS10pre <- read.table("BT_SCR_detect1993_1Knov132022.txt",header=T)
##add known fish weeks to be able to plot later on (same distance between weeks)
TagID <- c(464,36,36,36,226)
Week <- c(9,1,2,3,11)
RKM <- c(5317.8, 4000,4000,4000,33666.0)
xASpre <- data.frame(TagID, Week,RKM )
AS10pre1993 <- rbind(AS10pre,xASpre )
AS10 <- AS10pre1993%>%
  arrange(TagID,Week,RKM)%>%   
  mutate(RKM =RKM/1000)#change from meters to km 
############################################################################
##make abacus plot for 1993
#age2 fish Nradionum =9376,9356,"9316","9295","26","9136","36","26","9916"
detect <- read.csv("detect1993nov132022.csv",header=T)
#names(detect)

detectgg <- detect%>%
  mutate(Nradionum=as.character(Nradionum))
#make an abacus plot of raw detections
# ggformat<-(
#   ggplot(detectgg,aes(x=week,y=Nradionum,col = RELO_METHO))+
#     geom_jitter(size=1.25, shape = 20,width = 0.17,height =0.13) +
#     theme_bw()+
#     labs(title = "", x = "Weeks", y = "Fish")+
#     scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "May 26",
#                                                                   "3" = "Jun 9","5" = "Jun 23","7" = "Jul 7","9" = "Jul 22",
#                                                                   "11" = "Aug 4", "13" = "Aug 18","15" = "Sept 1","17" = "Sept 15"))+
#     #  limits = c(1, 18))+
#     theme( panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
#     theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=8,color="black"))+
#     theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=8,color="black"))+
#     theme(axis.title.y = element_text(size = rel(1.6), angle = 90))+
#     theme(axis.title.x = element_text(size = rel(1.6), angle = 00))+
#     theme(plot.margin = margin(0,1,0,0, "cm"))+
#     theme(legend.position = "none")
# )
# ggsave(filename="BTPRE1993abacusplot.png", plot=ggformat, device="png",
#        height=6, width=7, units="in", dpi=800)
############################################################################################
melted.rkm <- melt(AS10, id=c("TagID","RKM")) 
tagid.week <- cast(melted.rkm, TagID ~ RKM ~ value, fill=0, length)
y <- tagid.week
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
##############################################################################
##first week fish came into study
first.dat<-read.csv("1993firstcapnov13.csv",stringsAsFactors = F)
first <- first.dat %>%
  mutate(firstcap = ifelse(Nradionum ==36, "1",firstcap))%>%
  dplyr::select(firstcap)
first <- as.numeric(first$firstcap)
last <- c(16,17,18,18,16,15,16,18,15,16,15,18,18,17,18,15,15,16,18,18,18,18,15,15,16,16,15,15)
#week <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
#unique(AS10$TagID)
# firstflowA <- c(5.327107,20.420492,18.159190,18.159190,18.159190,
#                 13.122836, 10.796809,  7.868038,7.868038,7.868038,13.122836,
#                 13.122836,13.122836,13.122836,13.122836,13.122836,7.868038,7.868038,
#                 7.868038,7.868038,18.159190,18.159190,18.159190,13.122836,
#                 18.159190,18.159190,18.159190,18.159190)
# firstflowC <- c(-2.54093138,1.43202341,0.65533275,0.65533275,0.65533275,-5.03635353,
#                 -2.32602673,-2.92877105,-2.92877105,-2.92877105,-5.03635353,-5.03635353,
#                 -5.03635353,-5.03635353,-5.03635353,-2.92877105,-2.92877105,)  
#############################################################################
###AGE of fish in the study
##############################################################################
AGEpre <- read.csv("1993agenov13.csv",stringsAsFactors = F)
AGE <- AGEpre%>%
  select(age)
AGE <- as.numeric(AGE$age)
##############################################################################
##fish size standardized from prep file
fishsize.dat <-read.csv("fishsize1993dec22022.csv",stringsAsFactors = F) 
fishsize1993 <- fishsize.dat %>%
  select(FL)
sizeFL <- as.numeric(fishsize1993$FL)
##check and see if age and size are correlated an they are not
# x1 <- sizeFL
# y1 <- AGE
# cor(x1,y1)
################################################################################
#flowA92 is  average flow 1992 (weeks 1-18, not including week 0)
FlowA93 <- c(20.420492, 19.299954, 17.503857, 18.159190, 13.122836, 10.796809,  7.868038,
             5.327107,4.955448, 4.194939, 3.559832, 3.304981, 3.183623, 2.887105,
             2.630231,2.137113,2.038813, 1.996742)
###flowB is delta maximum - minimum within the week 1992
FlowB93 <- c(2.74673417, 9.88257966, 4.16257654, 1.72732768, 3.02990265, 1.98217930, 2.57683309,
             0.96277280, 0.48138640, 0.73623803, 0.59465379, 0.28316847, 0.39643586, 0.23786152,
             0.15291098, 0.50120819, 0.03681190, 0.05663369)
###flowC is delta among weekly maximum - minimum  week 1992
FlowC93 <- c(1.43202341, -1.12053809, -1.79609716,  0.65533275, -5.03635353, -2.32602673,
             -2.92877105, -2.54093138, -0.37165862, -0.76050961, -0.63510643, -0.25485162,
             -0.12135792, -0.29651784, -0.25687426,-0.4931176,-0.09829991, -0.04207074)
###Temp variables: see paper for context
#TempA is A) weekly maximum hourly water temperature avg
TempA93 <- c(9.50000,  9.50000,  9.40000, 10.50000, 10.60000, 11.40000,
             12.40000, 12.50000, 12.60000, 13.60000,13.70000, 12.00814, 12.55851,
             12.97647, 10.93544, 11.21063, 11.40655, 11.10232)
##Temp B is weekly average daily maximum water temperature avg
TempB93 <- c(9.085714,  9.000000,  8.757143, 10.000000, 10.057143, 10.514286,
             11.428571, 11.037500, 11.216667, 13.271429, 13.357143, 11.343617,
             12.272988, 12.310611, 10.429367, 10.605114, 11.139525,  9.960847)
###Temp C is weekly three-day moving average of the maximum daily temperature avg
TempC93 <- c(9.080000,  9.126667,  8.673333, 10.086667, 10.180000, 10.413333,
             11.693333, 11.105556, 11.075000, 13.226667, 13.413333, 11.190379,
             12.308137, 12.478730, 10.457723, 10.641740, 11.176986, 9.988438)
###Temp D is weekly average hourly temperature avg
TempD93 <- c(8.931765,  8.745238,  8.325000,  9.589286,  9.634524,  9.378571,
             9.870238,  9.836458, 10.240278, 11.667857, 12.329762, 10.692822,
             11.411288, 11.461764, 9.647976, 9.765491, 10.341771, 9.465176)
###Temp E is delta current - previous week average
TempE93 <- c(0.82819328, -0.18652661, -0.42023810,  1.26428571,  0.04523810,
             -0.25595238,  0.49166667, -0.03377976,  0.40381944,  1.42757937,
             0.66190476, -1.63693984,  0.71846594,  0.05047639, -1.81378874,
             0.11751567,  0.57627936, -0.87659433)
##scale temp and flow vars
flowA=(FlowA93-mean(FlowA93))/(sd(FlowA93))#standardize flow variable
flowB=(FlowB93-mean(FlowB93))/(sd(FlowB93))#standardize flow variable
flowC=(FlowC93-mean(FlowC93))/(sd(FlowC93))#standardize flow variable
tempA=(TempA93-mean(TempA93))/(sd(TempA93))#standardize flow variable
tempB=(TempB93-mean(TempB93))/(sd(TempB93))#standardize flow variable
tempC=(TempC93-mean(TempC93))/(sd(TempC93))#standardize flow variable
tempD=(TempD93-mean(TempD93))/(sd(TempD93))#standardize flow variable
tempE=(TempE93-mean(TempE93))/(sd(TempE93))#standardize flow variable
###############################################################################
##### test STANDARDIZED flow and temp correlations###########################################
#FlowB and FlowC uncorrelated with all Temp vars
#but flowA and flowB vars correlated and FlowA correlated with TempA,B,C,D
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
##################################################################################
###bring in covariates for reaches we want joined
JimCovRRa <- read.csv("JIMSCOVSnov3.csv",header=T)
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
  #mutate(slowhabs=(MP.area+LS.area+PL.area+PW.area+SP.area+UP.area+DP.area)/TotArea)%>%
  mutate(pools=(PL.area+PW.area+SP.area+UP.area+DP.area)/TotArea)%>%
  select(-LEN_VGST_R,-LEN_VGST_L,-LEN_UCUT_L,-LEN_UCUT_R,-LEFT_LEN,-RITE_LEN,-LWD_AGREGS,
         -LWD_WADS,-newRea1K.x,-SegmentLen.x,-newRea1K.y,-SegmentLen.y,-NR1K)
##do PCA on the habitat area variables to come up with a new variable
pcajimcov <-JimCovRR%>% 
  select(TOT_VGST,LEN_UCUT,dens.LWD_TOT)
modelVS<-prcomp(pcajimcov,scale=TRUE)
#	model.st<-prcomp(st.spdata,scale=TRUE)
modelVS$rotation[,1:3]
###to get the loadings for each data point for PC1 to 3
PCVS<-(modelVS$x[,1:2])
colnames(PCVS) <- c("PC1", "PC2")
#combine the PCs to the last dataset with coverpc
JimCovRR<-cbind(JimCovRR,PCVS)
JimCovRRst <-JimCovRR%>%
  mutate_at(c("Shan", "HG.prop", "RN.prop","CS.prop","LG.prop","TOT_VGST",
              "LEN_UCUT","PC1","dens.LWD_TOT","MeanDepth","pools"), ~(scale(.) %>% as.vector))
##################################################################################
Cov1 <- JimCovRRst$Shan
#Cov2 <- JimCovRRst$HG.prop
#Cov3 <- JimCovRRst$RN.prop
Cov2 <- JimCovRRst$CS.prop
#Cov3 <- JimCovRRst$LG.prop
#Cov6 <- JimCovRRst$LEN_UCUT
Cov3 <- JimCovRRst$dens.LWD_TOT
#Cov8<- JimCovRRst$MeanDepth
Cov4 <- JimCovRRst$TOT_VGST 
#Cov4 <- JimCovRRst$PC1 #PC1 for overhead cover
#Cov11 <- JimCovRRst$pools#pool
###############################################################
#####  WRITE TEXT FILE WITH JAGS SPECIFICATION   ###########
###################################################################
sink("ModelCJS.txt")
cat("
model { #model
# Priors
alphasig ~ dnorm(0,0.368) 		# Intercept in sigma estimate  
alphaphi ~ dnorm(0,0.368)     #intercept for phi (needs to be positive)
alphalam ~ dnorm(0,0.368)    #intercept for lam0
alphalam1 ~ dnorm(0,0.368)    #intercept for lam0 first encounters
betaT2 ~ dnorm(0,0.368)       #Coefficient for temp in S estimate 
  # betaF  ~ dnorm(0,0.368)      #Coefficient for flow in S estimate 
betaphiFlow  ~ dnorm(0,0.368) #beta for phi for flow
  #betaphipos  ~ dnorm(0,0.368) #beta for phi for S
#betaphiTemp ~ dnorm(0,0.368)  #beta for phi for temp
  #betaphiTemp2 ~ dnorm(0,0.368)  #beta for phi for temp
  #betaagephi ~ dnorm(0,0.368)  #age on survival
betaageS ~ dnorm(0,0.368)      # age on movements
 #betaFL.phi ~ dnorm(0,0.368)   #beta for phi for size
 #betaFL.S ~ dnorm(0,0.368)     #beta for S for size
tauv  ~ dgamma(0.1,0.1) 	    #spread for correlated S's (needs to be positive)
tau <- 1/(tauv*tauv)  	      #spread (precision) for correlated S's

#priors reach encounter betas
            for(zz in 1:4){
              gamma[zz]~dnorm(0,0.368)
            }
# # #priors flow and temp betas
            # for(xx in 1:17){
            #   tfgamma[xx] ~ dnorm(0,0.368)
            # }

# # ##indicator variable prior
#          for(yy in 1:32){
#           w[yy] ~ dbern(0.5)
#          }
##sigma is how far fish can be available to move into reach from other reaches
      for (t in 1:T){  # t = primary sample occasion (Weeks)
  	      log(sigma[t]) <- alphasig
  	          #log(sigma[t]) <- alphasig + betasig * flowC[t]#*temp[t]+ betasig2 *temp[t]
    	     	  ##below for model selection
    	        # log(sigma[t]) <- alphasig + w[12]*tfgamma[1]*flowA[t] + w[13]*tfgamma[2]*flowC[t]+
    	        #       w[14]*tfgamma[3]*tempE[t] + w[15]*tfgamma[4]*tempB[t]
    	    sigma2[t] <- sigma[t] * sigma[t]  # Weekly sigma2
	            } #t
####################################	        
      for (i in 1:M){ # m individuals, first are period of first tagging   
  	 for(t in 1:(first[i]-1)) { #t   periods prior to river entry
	                S[i,t] <- 0  # Fix S to 0 (not in river), needed to follow node in jags
	                z[i,t] <- 0  # Fit z to 0 (not in river), needed to follow node in jags
                          } #t
  #############################
            for(t in (last[i]+1):T) {#t  Periods after documented censor
	                    S[i,t] <- 0  # Fix S to 0 in river and known alive but not seen again this timeframe
	                    z[i,t] <- 0  # Fix z to 0 same as above needed to follow node in JAGS
               	}#t
  	 
  	    z[i,first[i]]  ~ dbern(1)  # state = Known to be alive at entry into study area
  	    S[i,first[i]]   ~ dunif(4,40)  # this is where fish enter study was reach RKM 4
################################
              for(j in 1:nreach) {  #j locations 
	             D2[i,j,first[i]] <- pow(S[i,first[i]]-reach.loc[j], 2)   
               log(lam0[i,j,first[i]]) <- alphalam1
               lam[i,j,first[i]] <- lam0[i,j,first[i]]*exp(-D2[i,j,first[i]]/
                                       (2*sigma2[first[i]])) 
               tmp[i,j,first[i]] <- lam[i,j,first[i]] 
               y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])
               } #j     
#######################################
        for (t in (first[i]+1):last[i]) { #t are weeks
            S[i,t] ~ dnorm(S[i,t-1] + betaT2*tempB[t-1]+betaageS*age[i],tau)T(4,40)#truncated normal
              ##below for model slection  
              #   S[i,t] ~ dnorm(S[i,t-1] + w[16]*tfgamma[5]*flowA[t] + w[17]*tfgamma[6]*flowC[t] +
              # 	   w[18]*tfgamma[7]*tempE[t] + w[19]*tfgamma[8]*tempB[t]+w[29]*betaFL.S*sizeFL[i]+
              # 	    w[31]*betaageS*age[i],tau) #
   	    
  #activity center is a funciton of where it was last time step and flow...with 
  #some error/variability (tau)...you can truncate this to make sure they stay in 
  #places it is possible to go (i.e., if they cant get through a barrier on the river). 
  #You could also add an intercept if you think there is some constant change in the 
  #activity center across time periods.
         		                          
          for(j in 1:nreach) { #j locations
		          D2[i,j,t] <- pow(S[i,t]-reach.loc[j], 2)   
                ##below for model slection
                # log(lam0[i,j,t]) <- alphalam +  w[1]*gamma[1]*Cov1[j] + w[2]*gamma[2]*Cov2[j] +
                #        w[3]*gamma[3]*Cov3[j] + w[4]*gamma[4]*Cov4[j] + w[5]*gamma[5]*Cov5[j]+
                #        w[6]*gamma[6]*Cov6[j] #+ w[7]* gamma[7]*Cov7[j] + w[8]*gamma[8]*Cov8[j] +
                #        #w[9]*gamma[9]*Cov9[j] + w[10]*gamma[10]*Cov10[j]+ w[11]*gamma[11]*Cov11[j]+
                #        w[20]*tfgamma[14]*flowA[t] + w[21]*tfgamma[15]*flowC[t] + 
                #       w[27]*tfgamma[16]*tempE[t] + w[28]*tfgamma[17]*tempB[t]

                log(lam0[i,j,t]) <- alphalam + gamma[1]*Cov1[j] + gamma[2]*Cov2[j] +  gamma[3]*Cov3[j]+
                    gamma[4]*Cov4[j] #+ gamma[5]*Cov5[j]
                lam[i,j,t] <- lam0[i,j,t] * exp(-D2[i,j,t]/(2*sigma2[t]))#pdf
	              tmp[i,j,t] <- z[i,t]*lam[i,j,t]
	              y[i,j,t] ~ dpois(tmp[i,j,t])
	                    } #j
    
               ###put a logistic regression on phi here   ########################################## 
              logit(phi[i,t])<- alphaphi + betaphiFlow*flowC[t] #+betaphiTemp*tempB[t]
                      ##Use below for model selection with w indicators
                      #   logit(phi[i,t])<- alphaphi + w[22]*tfgamma[9]*S[i,t] + w[23]*tfgamma[10]*flowA[t] +
                      #      w[24]*tfgamma[11]*flowC[t] + w[25]*tfgamma[12]*tempE[t] + w[26]*tfgamma[13]*tempB[t]+
                      #    w[30]*betaFL.phi*sizeFL[i]+ w[32]*betaagephi*age[i]
              phi2[i,t] <- max(0.0001,min(phi[i,t],0.9999))#prevents outside 0 and 1
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
JAGSdata<-list(y=y,first=first,M=M,T=T,nreach=nreach,tempB=tempB,tempE=tempE,
               reach.loc=reach.loc,flowC=flowC,age=AGE,last=last,#sizeFL=sizeFL,flowA=flowA,
               Cov1=Cov1,Cov2=Cov2,Cov3=Cov3,Cov4=Cov4#,Cov5=Cov5#,Cov6=Cov6#,
               #Cov7=Cov7,Cov8=Cov8,Cov9=Cov9,Cov10=Cov10,Cov11=Cov11
) 

z=matrix(NA, M, T)#latent variable
for(i in 1:M){
  for(t in first[i]:last[i]){
    z[i,t] <-1
  }
}

#Set initial values
inits =  function() {list(z=z,alphaphi=runif(1,0,1)#,alphalam1=runif(1,-15,3),
                          # tauv=runif(1,0,20),alphasig=runif(1,-3,3),alphalam2=runif(1,-3.5,3.5),
                          # betasig=runif(1,-1,1),betaF= runif(1,-2,2),betaFL.phi=runif(1,-1,1),
                          # betaFL.S=runif(1,-2,2),betaphiFlow =runif(1,0,25),betaphiTemp =runif(1,0,22),
                          # betaT2=runif(1,-2,2),betasig2=runif(1,-3,3)
)#,gamma=runif(11,-0.4,0.55),alphaS=runif(1,14,20
}
#Parameters to follow, REMEMBER to add S, phi, and z for final run
parameters <- c("alphaphi","alphalam","alphasig","betaageS",
                "betaphiFlow","betaT2","alphalam1",
                "tauv", "gamma","S","phi","z","lam0","sigma"
)#parameters  
#cant get sigma in the same output as the rest because too large  

#Call JAGS
Sys.time()

ZZ<-jags(data=JAGSdata,inits=inits,parameters.to.save=parameters,model.file="ModelCJS.txt",n.thin=2,
         n.chains=3, n.adapt = 2000,n.burnin=6000,n.iter=20000,parallel=TRUE)#
Sys.time()
###############################################################################
#check convergence in trace plots
save(ZZ,file="ZZ_1993.Rdata")##save as robject
MCMCtrace(ZZ, params = c("alphaphi","alphalam","alphasig","alphalam1",
                         "betaphiFlow","betaT2", "betaageS",
                         "tauv", "gamma"
))

##get deviance and parameter estimes with MCE
mco2 <- mcmcOutput(ZZ, params=c("alphaphi","alphalam","alphasig","alphalam1",
                                "betaphiFlow","betaT2","betaageS",
                                "tauv", "gamma[1:3]",
                                "S","phi","z","lam0",'sigma'))

mco2summary <- summary(mco2,CRImass=0.90)#summary with rhats, lcl, ucl for updating suppl table estimates

load("ZZ_1993.Rdata")
library(MCMCvis)#load bayes plotting package
#make caterpillar plot from MCMCvis package
##!!READ THIS!!! To save these plots you must go to the plot viewer,
##choose export: save as pdf, PDF size = 5x7, choose landscape and use cairo pdf device
##then open the pdf go to file:export to:image:png
##this will give a high quality and consistent image....
RStudioGD() #this will set the plot device (dev) to the plot pane in rstudio
MCMCplot(ZZ,params = c("alphasig","betaageS", "betaT2","alphalam","alphalam1",
                       'gamma','alphaphi', "betaphiFlow","tauv"),ci = c(50, 90), guide_lines = TRUE,HPD = TRUE,
         xlab = 'Parameter estimate',sz_med = 1.1,sz_labels = 1.0,sz_ax = 1.1,
         pos_tick=c(-3,-2,-1,0,1,2,3,4,5,6),sz_tick_txt=0.8,sz_ax_txt=1.1,
         labels = c("Scale int","Movement flow", "Movement maturity",
                    "Detection int 1","Detection int 2",'Habitat diversity index', 
                    'Prop casc, falls,riff','Prop veg stable bank', 
                    "Density LWD","Survival int","Survival flow","Precision")) 

##for plotting 1992 and 1993 together, could be nice
# MCMCplot(object = MCMC_data, 
#          object2 = MCMC_data2, 
#          params = 'beta', 
#          offset = 0.1,
#          ref_ovl = TRUE)
###############################################################################
###1. model selection for indicator variables
################################################################################
#this example has 11 covariate indicator variables w 
# mod <- mco2$w
# mod <- paste(mod[,1],mod[,2],mod[,3],mod[,4],mod[,5],mod[,6],
#              mod[,7],mod[,8],mod[,9],mod[,10],mod[,11],
#              mod[,12],mod[,13],mod[,14],mod[,15],mod[,16],
#              mod[,17],mod[,18],mod[,19],mod[,20],mod[,21],mod[,22],
#              mod[,23],mod[,24],mod[,25],mod[,26],mod[,27],
#              sep="")
# table(mod)#creates list of models
# #sort(round(table(mod)/length(mod), 5))#creates AIC table COOL!!!!
# restable <-round(table(mod)/length(mod), 5)
# restable <- as.data.frame(restable)
# weights <-  restable %>%
#   arrange(desc(Freq))
# write.csv(weights,"post_weights1993_.csv", row.names = F)
#######################################################################
######################################################################
###create the draws for S, z, and phiUP to use in plots
######################################################################
Smov <- mco2$S
zstate <- mco2$z 
sigma <- mco2$sigma
phi <- mco2$phi
lam0 <- mco2$lam0
##########get sigma estimates for manuscript with 90% CIs
(meansig <- mean(sigma,na.rm = TRUE))
(lCLsig <- quantile((sigma),probs= 0.10,na.rm = TRUE)) 
(UCLsig <- quantile((sigma),probs= 0.90,na.rm = TRUE))
############################################################
#extrapolate total phi for the entire time period for the MS with CI
#(phi[,x,y])
appStot= NULL
fish = NULL
for (j in 1:21000){
  for(f in 1:28){
    fish[f] <-  prod(phi[j,f,1:18], na.rm=TRUE)
  }
  appStot[j] <- mean(fish, na.rm=T)
  fish <-NULL
}

mean(appStot)
quantile(appStot,probs = c(0.1,0.9))


#get the draws into a format we can work with
counter <-1
m_predict <- matrix(0,nrow=504, ncol = 7)
for (i in 1:28) {
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
  ##take out weeks that only had one fish
  filter(week != "1"&week != "2" & week !="3"&week!="4"&week!="18")%>%
  group_by(week)%>%
  ###calculate the total survival for the entire period  
  summarize(weekphi=mean(avg.phi),LCLphi=mean(LCL.phi),UCLphi=mean(UCL.phi))
mean(wk.phi$weekphi)
mean(wk.phi$LCLphi)
mean(wk.phi$UCLphi)
wk.phi$week2 <- wk.phi$week-4

FlowC93_17w <- c( -5.03635353, -2.32602673,-2.92877105, -2.54093138, -0.37165862, 
                  -0.76050961, -0.63510643, -0.25485162,
                  -0.12135792, -0.29651784, -0.25687426,-0.4931176, -0.09829991)
# TempB93_17w <- c(10.057143, 10.514286,
#                  11.428571, 11.037500, 11.216667, 13.271429, 13.357143, 11.343617,
#                  12.272988, 12.310611, 10.429367, 10.605114, 11.139525)

wk.phi$FlowC93 <- FlowC93_17w ##this is because no estimate fo survival for week one.
#wk.phi$TempB93 <- TempB93_17w

RStudioGD() 
##!!!!MAKE A PLOT FOR SURVIVAL BY WEEK with credible intervals!
ggformat<-(
  ggplot(wk.phi, aes(x = week2, y=weekphi,ymin = LCLphi, ymax = UCLphi)) +
    geom_point()+
    geom_pointrange(color="dodgerblue",size=0.5, shape = 20)+
    geom_line(aes(y = FlowC93_17w/7.5), color="blue4")+
    theme_bw()+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13),labels=c("1" = "Jun 23","3" = "Jul 7","5" = "Jul 22",
               "7" = "Aug 4", "9" = "Aug 18","11" = "Sept 1","13" = "Sept 15"),
             limits = c(1, 13))+
    
    scale_y_continuous("Mean survival",#limits = c(-0.75, 1),
                       sec.axis = sec_axis(~ . *7.5,name = "Flow (cms)"),
                       )+#breaks = seq(-5,0, by=5)
    labs(title = "",  x = "")+
    
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=9.5,color="black"))+
    theme(axis.text.y = element_text(angle = 0, hjust = 0.2, size=11,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.5), angle = 00))+
    theme(plot.margin = margin(0,1,0,0, "cm"))
)
ggsave(filename="BTPRE1993weeklyphiOct292023.png", plot=ggformat, device="png",
       height=4, width=6, units="in", dpi=800) 
##############################################################
##Plot for movement for fish that are very likely alive
#make PLOT of FISH MOVEMENTS (by fish)
################################################################
#filter out NAs and z<0.99
movements <- new_frame%>%
  filter(z!="NA")%>%
  filter(z>0.7)
#make dataframe to plot movements with temp NOT flow
  flowC93p <- FlowC93[1:17]
  temp93p <- TempB93[1:17]
  flow <- as.data.frame(flowC93p)
  flow$week <- sort(unique(movements$week))
  dfpre <- left_join(flow,movements)
  temp <- as.data.frame(temp93p)
  temp$week <- sort(unique(movements$week))
  df <- left_join(temp,dfpre)
####################################################################################
###all fish together
###################################################################################
df$fish <- as.factor(df$fish)
ggformat <- (
  ggplot(df, aes(x = week, y=avg.S+56, color=fish)) +
    geom_point( size = 1)+
    # geom_line(aes(y = FlowC93*7.5), color="blue4")+
    geom_line(aes(y = temp93p*2.5), color="red4")+
    #scale_color_viridis_d()+
    scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15,17),labels=c("1" = "May 26",
          "3" = "Jun 9","5" = "Jun 23","7" = "Jul 7","9" = "Jul 22",
          "11" = "Aug 4", "13" = "Aug 18","15" = "Sept 1","17" = "Sept 15"),
             limits = c(1, 18))+
    scale_y_continuous("Reach location (RKM)",limits = c(0, 102),
                       sec.axis = sec_axis(~ . /6.5,name = "Temp (°C)",
                          breaks = c(0,5,10,15)))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size=8.5,color="black"))+
    theme(axis.text.y = element_text(angle = 0, hjust = 0, size=10,color="black"))+
    theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
    theme(axis.title.x = element_text(size = rel(1.2), angle = 00))+
    labs(title = "", x = "")+
    theme(plot.margin = margin(0,1,0,0, "cm"))+
    theme(legend.position = "none")
)
ggsave(filename = "BT_Movements1993oct28_1993.png", plot=ggformat,
       device="png",height=4,width=6,units = "in",dpi=800)
############################################################################
##make plot lam0 by reach
############################################################################
counter <-1
m_covar <- matrix(0,nrow=20160, ncol = 8)
for (i in 1:28) {
  for (j in 1:40){
    for (k in 1:18){
      m_covar[counter,1]<- mean(lam0[,i,j,k])#mean lam0
      m_covar[counter,2]<-k
      m_covar[counter,3]<-j
      m_covar[counter,4]<-i
      m_covar[counter,5]<- mean((lam0[,i,j,k]))
      m_covar[counter,6]<- y[i,j,k]
      m_covar[counter,7]<- quantile((lam0[,i,j,k]),probs= 0.05,na.rm = TRUE) 
      m_covar[counter,8]<- quantile((lam0[,i,j,k]),probs= 0.95,na.rm = TRUE)
      counter <- counter + 1
    }
  }
}
m_covar<-as.data.frame(m_covar)
colnames(m_covar) = c("pred.lam0","week","reach","fish","explam0","datay","LCIpredlam0","UCIpredlam0")
##then do group by reach to get means of lam0
lam0reach <- m_covar%>%
  filter(pred.lam0 != "NA")%>%
  filter(reach != "0")%>%
  group_by(reach)%>%
  summarize(lam0reach=mean(pred.lam0),LCI=mean(LCIpredlam0),UCI=mean(UCIpredlam0))
lam0reach$dist <- round(reach.loc+56,0)
lam0reach<-as.data.frame(lam0reach)

mean(lam0reach$lam0reach)
mean(lam0reach$LCI)
mean(lam0reach$UCI)
#####################################################################
##MAKE A PLOT FOR lam0 BY REACH with credible intervals!
RStudioGD()

ggformat<-(
  ggplot(lam0reach, aes(x = reach, y=lam0reach, ymin = LCI, ymax = UCI)) +
    geom_point()+
    geom_pointrange(color="black",size=0.5, shape = 20)+
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
    labs(title = "", y = "Mean detection rate", x = "Rapid River Reach (RKM)")
)
ggsave(filename="BTPRE1993reachlam0Oct28_2023.png", plot=ggformat, device="png",
       height=5, width=5, units="in", dpi=800) 
