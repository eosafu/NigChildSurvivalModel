# load exp_cox_gamma_weibull_lognorm_model_new.RData into workspace

###3 Kapler Meier
library(ggplot2)
library(survival)
surv.vet <- Surv(dat.surv2$dayofD,dat.surv2$delta)
km.vet <- survfit(surv.vet ~ 1)
plot(km.vet)


## plot kapler meier

dat.kapler <- data.frame(time=km.vet$time, lw=km.vet$lower,surv=km.vet$surv,upp=km.vet$upper)
dat.kapler %>% ggplot(aes(x=time))+
  geom_line(aes(y=surv))+
  geom_line(aes(y=lw),colour="red",linetype = "dashed")+
  geom_line(aes(y=upp),colour="red",linetype = "dashed")+theme_minimal()+
  labs(x="Months",y="Survival probability")

dat.kapler <- data.frame(time1=km.vet$time[1:598], lw1=km.vet$lower[1:598],surv1=km.vet$surv[1:598],upp1=km.vet$upper[1:598],
                         time2=c(km.vet$time[599:1172],rep(NA,24)), lw2=c(km.vet$lower[599:1172],rep(NA,24)),surv2=c(km.vet$surv[599:1172],rep(NA,24)),upp2=c(km.vet$upper[599:1172],rep(NA,24)))
dat.kapler %>% ggplot()+
  geom_line(aes(x=time1,y=surv1,colour="Male"))+
  geom_line(aes(x=time1,y=lw1),colour="red",linetype = "dashed")+
  geom_line(aes(x=time1,y=upp1),colour="red",linetype = "dashed")+
  geom_line(aes(x=time2,y=surv2,colour="Female"))+
  geom_line(aes(x=time2,y=lw2),colour="orange",linetype = "dashed")+
  geom_line(aes(x=time2,y=upp2),colour="orange",linetype = "dashed")+theme_minimal()+
  labs(x="Months",y="Survival probability")+

scale_colour_manual("", 
                    breaks = c("Female", "Male"),
                    values = c("black", "blue")) +
  labs(title="")

#%%%%%%%%%%%%%%%
# Plot model results
#%%%%%%%%%%%%%%
require(rgdal)
fn <- file.path(tempdir(), "nigeria-lgas.zip", fsep = "\\")
download.file("https://s3-eu-west-1.amazonaws.com/cfa-openafrica/resources/372a616a-66cc-41f7-ac91-d8af8f23bc2b/nigeria-lgas.zip?AWSAccessKeyId=AKIAJZDSKF5QFOBBR4RQ&Expires=1618247755&Signature=dcSAuxTqhtkTzlfr4Tpn5g5eyJI%", fn)
utils::unzip(fn, exdir = tempdir())
shp <- readOGR(dsn ="Nigeria_surv/nigeria-lgas/new_lga_nigeria_2003.shp", stringsAsFactors = F)

shp.map <- fortify(shp, region = "STATE")

shp.map$count=0

state2 <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
            "Benue.1",       "Borno.1" ,      "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
            "Enugu.1"   ,    "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
            "Kebbi.1",       "Kogi.1"  ,      "Kwara.1"  ,     "Lagos.1" ,    "Nassarawa.1" , 
            "Niger.1",       "Ogun.1" ,       "Ondo.1" ,       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
            "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")
cloc2 <- c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,
           28,18,19,20,21,29,30,36)
spat.lognorm <- lognorm.vet$summary.random$state[1:37,]


est_data <- data.frame(region=state2,est=spat.lognorm$mean[cloc2],est.0.025quant=spat.lognorm$`0.025quant`[cloc2],est.0.975quant=spat.lognorm$`0.975quant`[cloc2])


s <- function(region,est){
  shp.map$count[which(shp.map$group==as.character(region))] <- est
  return(shp.map)
}

# Mean 

for (i in 1:37) {
  shp.map <-s(est_data[i,1],est_data[i,2]) 
}
shp.map$count[which(shp.map$group=="Lagos.2")] <- est_data[25,2]
shp.map$count[which(shp.map$group=="Lake.1")] <- est_data[25,2]

ggplot(shp.map, aes(x = long, y = lat, group = group, fill = count)) +
  geom_polygon(colour = "black", size = 0.5, aes(group = group)) +
  theme_minimal()+ scale_fill_gradient(name="Risk",high='green', low='red')

### Probability 

shp.map$prob  <- 0
est_data$prob <- 0

Probab <- function(k){
  
  m <-  lognorm.vet$marginals.random$state[[k]]
  m <- inla.pmarginal(0,m)
  return(m)
  
}
probb <- unlist(lapply(1:37,Probab)) 

est.t_data <- data.frame(region=state2, prob4=probb[cloc2])

#est_data$prob <- 0

s <- function(region,est){
  shp.map$count[which(shp.map$group==as.character(region))] <- est
  return(shp.map)
}

# Mean 

for (i in 1:37) {
  shp.map <-s(est.t_data[i,1],est.t_data[i,2]) 
}
shp.map$count[which(shp.map$group=="Lagos.2")] <- est_data[25,5]
shp.map$count[which(shp.map$group=="Lake.1")] <- est_data[25,5]

ggplot(shp.map, aes(x = long, y = lat, group = group, fill = count)) +
  geom_polygon(colour = "black", size = 0.5, aes(group = group)) +
  theme_minimal()+ scale_fill_gradient(name="Prob.",high='darkred', low='orange')


#### Non Linear 
nolinear.lognorm <- lognorm.vet$summary.random$M_age
nolinear.lognorm <- data.frame( id=nolinear.lognorm$ID,mean=nolinear.lognorm$mean,
                                lw= nolinear.lognorm$`0.025quant`,
                                upp=nolinear.lognorm$`0.975quant`)
nolinear.lognorm  %>% ggplot(aes(x=id))+
  geom_line(aes(y=mean))+
  geom_line(aes(y=lw),colour="red")+
  geom_line(aes(y=upp),colour="red")+theme_minimal()+
  labs(x="Mother's age")


nolinear.lognorm <- lognorm.vet$summary.random$birthweight
nolinear.lognorm <- data.frame( id=nolinear.lognorm$ID,mean=nolinear.lognorm$mean,
                                lw= nolinear.lognorm$`0.025quant`,
                                upp=nolinear.lognorm$`0.975quant`)
nolinear.lognorm  %>% ggplot(aes(x=id))+
  geom_line(aes(y=mean))+
  geom_line(aes(y=lw),colour="red")+
  geom_line(aes(y=upp),colour="red")+theme_minimal()+
  labs(x="Child's weight")

#### Fixed Effect
library(xtable)
xtable(lognorm.vet$summary.fixed[c("mean","sd","0.025quant", "0.975quant")],digits=3)

xtable(data.frame(
  cox.vet$summary.fixed[c("mean","0.025quant", "0.975quant")],
  exp.vet$summary.fixed[c("mean","0.025quant", "0.975quant")],
  gamma.vet$summary.fixed[c("mean","0.025quant", "0.975quant")],
  weibull.vet$summary.fixed[c("mean","0.025quant", "0.975quant")]),
  digits=3
)

#### density plot of time to death
dat.dens <- data.frame(x=dat.surv2$dayofD[dat.surv2$delta==1])
dat.dens%>% ggplot(aes(x=x,fill="xgreen")) +
  geom_density(alpha=0.4)+labs(x="Months")+theme_minimal()+
  theme(legend.position = "none")

###########
#DHS.component probabilities
A=read.dta("C:\\Users\\11494792\\Downloads\\NG_2018_DHS_10122021_214_117034\\NGBR7BDT\\NGBR7BFL.DTA")
B=chmortp(A)
age= row.names(B)
BB=data.frame(agegroup=age,surv.cumm=cumsum(B$PROBABILITY),B)
p<-ggplot(data=BB, aes(x=reorder(agegroup, surv.cumm), y=surv.cumm)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+ theme(axis.text.x=element_text(angle=45, hjust=1))
p=p+geom_line(data=BB, aes(x=as.numeric(reorder(agegroup, surv.cumm)), y=surv.cumm), colour = "red")
p+labs(x="Age group (Months)",y="Cummulative component probabilities")

ggplot(data=data.frame(y=ch,x=dat.surv2$dayofD))+
  geom_line(aes(x=x, y=y),size=1,color="red")+theme_minimal()+
  labs(x="Months",y="Cummulative hazard rates")
