#  data arrangement
library(tidyverse)
library(INLA)
library(BayesX) # for bnd2sp() & read.bnd()
library(spdep) # for poly2nb()
#dat.surv <- read.csv("Cleaned_data_new_m.csv") # Load the sample dataset (dat.surv.Rdata) from this repository.
dat.surv %>% names


#%%%%%%%%%%%%%%%%%%%%%
# Data preparation
#%%%%%%%%%%%%%%%%%%%%%

# Convert death day to the same scale 

###################
# Note (https://dhsprogram.com/pubs/pdf/DHSG4/Recode5DHS_23August2012.pdf)
#Age at death (AaD) of the child as reported in the questionnaire. The first digit of the age at death
#gives the units in which it was reported: 1 - Days, 2 - Months, 3 - Years, 9 - Special
#responses. The last two digits give the age at death in those units. Age at death is usually
#reported in days if it was less than one month, in months if it was less than two years and
#otherwise in years. 
##################
dat.surv$AaD %>%unique
##########
# Convert Age at death (AaD) to Day of death (dayofD)
#########

dayofD <- delta <- vector("numeric",length=length(dat.surv$AaD))

for(i in 1:length(dat.surv$AaD)){
  
  if(!is.na(dat.surv$AaD[i])){
    
    if(dat.surv$AaD[i] < 200){
      
      dayofD[i] = dat.surv$AaD[i]-100 # If AaD is recorded in days
      delta[i]  = 1
      
    }else if(dat.surv$AaD[i] > 200 &dat.surv$AaD[i] <300){
      
      if (!is.na(dat.surv$DoD[i])){
      
      dayofD[i] = ((dat.surv$AaD[i]-200)-1)*30+dat.surv$DoD[i] # If AaD is recorded in months
      delta[i]  = 1
      }else{
        dayofD[i] = ((dat.surv$AaD[i]-200)-1)*30
        delta[i]  = 1
      }
      
    }else if(dat.surv$AaD[i] > 300 & dat.surv$AaD[i] <400){
      
      if(!is.na(dat.surv$DoD[i])&!is.na(dat.surv$MoD[i])){
      dayofD[i]=(dat.surv$AaD[i]-300-1)*365+dat.surv$MoD[i]*30+dat.surv$DoD[i] # If AaD is recorded in years.
      delta[i]  = 1
      }else {
        dayofD[i]=(dat.surv$AaD[i]-300-1)*365 
        delta[i]  = 1
      }
      
    }
  }else {
    dayofD[i]=NA
    delta[i]  = 0
  }
}
# Merge data set  
dat.surv1  = dat.surv %>%as_tibble()%>% mutate(dayofD=(dayofD)/30,delta=delta)
# Find children alive
id.na2 <-   is.na(dat.surv1$dayofD) 
# Assign those children alive their current age (in months) at the survey date
dat.surv1$dayofD[id.na2]<- dat.surv1$C_age[id.na2]
# Assign the follow-up time (59) months for children who were still alive and older than 59
id.na2 <-   is.na(dat.surv1$dayofD) 
dat.surv1$dayofD[id.na2]<- 59
# Drop any abnormal data point
dat.surv1 <- dat.surv1[-which(dat.surv1$dayofD>59),]

dat.surv2  = dat.surv1 %>%as_tibble()%>%select(-c("AaD","DoD","MoD" ,
                                                  "YoD",
                                                  "DoB","C_age" 
                                                  ))

#%%%%%%%%%%%%%%%%%%%
# Modeling begins
#%%%%%%%%%%%%%%%%%%%
dat.surv2%>%names

#----Map for Spatial component

map=read.bnd("nigeria.bnd")  # Download from repository
plot(map)
map_shp=bnd2sp(map)

temp <- poly2nb(map_shp)
nb2INLA("LDN.graph", temp)

H=inla.read.graph(filename="LDN.graph") # or download from repository
A <- inla.graph2matrix(H)
image(A,xlab="",ylab="")


# Using INLA framework

d.add = 1E-4
shyp = list(prec = list(prior = "pcprec", param = c(1,0.01)))


formula = sinla.vet ~ -1 +sex_female+ edu_pri +edu_sec +edu_high+wealth_poorer+wealth_mid+
  wealth_rich+wealth_richest+water_improve+toilet_improve+
  resid_rural+mas_yes+sex_female+antenat_no+antenat_fiveplus+postnatal_yes+breast_nofeed+
  place_deliv_hm+mother_work+married+seperated+catholic+islam+trad+other_rel+multiple_birth+
  firsthbirth+fivebirthplus+caesarean+tetanus_yes+bednet+sleepundernet+wifehusband+spouse+smelse+access_health+
  +                  f(M_age,model = "rw2")+
                     f(birthweight,model = "rw2")+
                     f(birthweight2,model = "rw2")+
                     f(state, model = "bym",
                       graph = "LDN.graph",
                       scale.model=T, diagonal = d.add, hyper = "shyp")


sinla.vet <- inla.surv((dat.surv2$dayofD+0.1),dat.surv2$delta)
  
#####%%%%%
# Using the Exponential hazard model
####%%%%%%
exp.vet <- inla(formula, data = dat.surv2, family = "exponentialsurv",
                control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))
#####%%%%%
# Using the Cox hazard model
####%%%%%%
cox.vet <- inla(formula, data = dat.surv2, family = "coxph",
                    control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))

#####%%%%%
# Using Gamma hazard model
####%%%%%%
n <- nrow(dat.surv2)
prec.scale <- runif(n, min = 0.001, max = 8)

gamma.vet <- inla(formula, data = dat.surv2, family = "gammasurv ",
                  scale = prec.scale,
                    control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))
#####%%%%%
# Using the Weibull hazard model
####%%%%%%
weibull.vet <- inla(formula, data = dat.surv2, family = "weibullsurv",
                  control.family = list(variant=0),
                  control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE))
#####%%%%%
# Using the lognormal hazard model
####%%%%%%
lognorm.vet <- inla(formula, data = dat.surv2, family = "lognormalsurv",
                    control.family = list(variant=0),
                    control.compute = list(dic=TRUE, cpo=TRUE, waic=TRUE),verbose=TRUE)
