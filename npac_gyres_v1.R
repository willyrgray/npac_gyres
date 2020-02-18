#from Gray et al 2020 GRL
#analysis of north pacific deglacial planktonic foram d18O to trace gyre boundary position
require(mgcv) #the mgcv package needs to be installed

#import data
dat<-read.csv("npac_d18o_data.csv")

#add identifier 
id<-matrix(, nrow = nrow(dat), ncol = 1)
for (i in 1:nrow(dat)){
id[i]<-paste(dat$core[i],"_",dat$species[i],sep="")
}

dat<-cbind(dat,id)

#select geographical domain 
#to run code for all basin or just east/est of 180 degrees E set domain to 'ALL', 'WEST', or 'EAST 
DOMAIN<- 'ALL'
#DOMAIN<- 'WEST'
#DOMAIN<- 'EAST'

dat<- 	if(DOMAIN == 'ALL') {dat} else 
		if(DOMAIN == 'WEST') {subset(dat, lon_e > 0)} else 
		if(DOMAIN == 'EAST') {subset(dat, lon_e < 0)}


id_sum<-as.vector(unique(dat$id))

#import sea level/sst data
#1 per mil ice volume plus global SST change from PMIP3 modelscaled to SL
esl_dat<-read.csv("lambeck2014_esl.csv") #sea level data
global_sst_dat<-read.csv("PMIP3_global_mean_DSST.csv") #pimp3 global mean DSST 

# age settings
age_holo<- 10.5 #sets Holocene age (too few sites younger than 10.5 ka)
start_age<- if(DOMAIN == 'WEST') {18} else {18.5} #age at which to start analysis
end_age<- age_holo
int<- -0.1 #time steps in ka at which to perform the analysis
int_age<- seq(start_age,end_age,by=int)

#lat settings
delta_lat_min<- -2 
delta_lat_max<- 10
lat_res <- 0.1
lat_s<- if(DOMAIN == 'EAST') {31.1} else {33.6} #5 or 10 degree lat window around steepest part of Holocene curve
lat_n<- if(DOMAIN == 'EAST') {41.1} else {38.6} 
kay<- if(DOMAIN == 'EAST') {3} else {8} #used to fit lat gam


### define shift function, which will be used later
shift <- function(x, n, invert=FALSE, default=NA){
  stopifnot(length(x)>=n)
  if(n==0){
    return(x)
  }
  n <- ifelse(invert, n*(-1), n)
  if(n<0){
    n <- abs(n)
    forward=FALSE
  }else{
    forward=TRUE
  }
  if(forward){
    return(c(rep(default, n), x[seq_len(length(x)-n)]))
  }
  if(!forward){
    return(c(x[seq_len(length(x)-n)+n], rep(default, n)))
  }
}


#########
#montecarlo loop for error propogation
mcmc_iterations<- 999 #999 iterations takes about 30 mins on my machine
s1<-matrix(, nrow = length(int_age), ncol = mcmc_iterations)

for (it in 1: mcmc_iterations) {
	
#resampling
dat_i<-dat[sample(nrow(dat), nrow(dat),replace=TRUE), ] #bootstrap resampling
dat_i$age<- rnorm(length(dat_i$age), dat_i$age, 0.5) #500 yr age uncertainties 
dat_i$d18o<- rnorm(length(dat_i$d18o), dat_i$d18o, 0.04) #0.04 per mil d18o uncertainty

#ice volume/global sst correction
Desl<- max(esl_dat$esl)-min(esl_dat$esl) #total sea level change
global_ivc_tot_per_mil<- rnorm(1,1,0.1) #lgm d18osw from Schrag et al 2002
global_sst_tot_per_mil<- rnorm(1, mean(global_sst_dat$DT), sd(global_sst_dat$DT))*-0.22 #from PMIP3/kimoneil97
Dd18Oesl_tot = global_ivc_tot_per_mil+global_sst_tot_per_mil
sl<- approx(x=esl_dat$age, y=esl_dat$esl, xout=int_age)$y

#model and predict deglacial d18O at int_age with gam
id_sum<-as.vector(unique(dat_i$id))
r1<-matrix(ncol=1,nrow=length(id_sum)) #vector of lats
r3<-matrix(ncol=length(int_age),nrow=length(id_sum))

for (i in 1:length(id_sum)){
	a <- as.character(id_sum[i])
	d<-subset(dat_i, id == a)
	lat<- dat_i$lat[match(a,dat_i$id)]
	r1[i,]<- lat
	if (age_holo >= min(d$age, na.rm=TRUE) & start_age <= max(d$age, na.rm=TRUE) & length(d$d18o) >= 11){
		k <- ifelse(length(d$d18o) > 30, 20, 10)
		b<-gam(d18o~s(age, k=k),data=d, method='REML')
		int_d18o<-predict(b,newdata=data.frame(age=int_age), se.fit=FALSE)
		Dd18Oesl<- Dd18Oesl_tot *(sl/Desl)
		r3[i,]<- 	int_d18o + Dd18Oesl } 
		else {r3[i,]<- NA}
	}

r4<-matrix(ncol=length(int_age),nrow=length(seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res)))

#fit lat gam to each time step
for (j in 1:length(int_age)){
	temp<- cbind(r1,r3[,j])
	temp<- na.omit(temp) 
	lat_int <- temp[,1]
	d18o_int <- temp[,2]
	b_int<-gam(d18o_int~s(lat_int, k=kay))
	pdat2<- data.frame(lat_int=seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res))
	p_int<-as.numeric(predict(b_int,newdata=pdat2, se.fit=TRUE)$fit)
	r4[,j]<- p_int
	}


#delta_lat calculation
pdat<- data.frame(lat=seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res))
q_holo<- r4[match(lat_s,pdat$lat):match(lat_n,pdat$lat),ncol(r4)]
r5<-matrix(ncol=1,nrow=length(int_age))

for (j in 1:length(int_age)){
	p_int<- r4[,j]
	
	delta_lat<- seq(delta_lat_min,delta_lat_max, by=0.1)
	y<-matrix(ncol=1,nrow=length(delta_lat))
	
		for (i in 1:length(delta_lat)){
		nshift<- delta_lat[i]*10
		p_shift <- shift(p_int,nshift)
		q<-p_shift[match(lat_s,pdat$lat):match(lat_n,pdat$lat)]
		y[i,] <- sum((q_holo - q)^2)
		}
	
	r5[j,]<- delta_lat[match(min(y),y)]
}

s1[,it]<- r5
}

#######
#loop over!

#find quantiles 
DLat<- apply(s1, 1, quantile, probs=0.5, na.rm=TRUE) #median value
DLat_upr95<- apply(s1, 1, quantile, probs=0.975, na.rm=TRUE) #95% range
DLat_lwr95<- apply(s1, 1, quantile, probs=0.025, na.rm=TRUE) #95% range
DLat_upr68<- apply(s1, 1, quantile, probs=0.84, na.rm=TRUE) #68% range
DLat_lwr68<- apply(s1, 1, quantile, probs=0.16, na.rm=TRUE) #68% range

#combine and export results
results<-data.frame(age=int_age,DLat=-DLat,DLat_upr95=-DLat_upr95,DLat_lwr95=-DLat_lwr95,DLat_upr68=-DLat_upr68,DLat_lwr68=-DLat_lwr68)
write.table(results,file=paste(DOMAIN,'_results.csv', sep=''),sep=',',row.names=FALSE, col.names=TRUE)


#quick wee plot!
dev.new(width=5.7, height=4); par(mar=c(3.5,3.5,0.5,0.5)); par(mgp=c(1.5,0.5,0)); par(las=1)
plot(-99, -99, col=adjustcolor("grey99",alpha=0.01),type='o',xlim=c(20,10),ylim=c(-5,1), pch=1, bty='n', axes=TRUE, xlab='Age (ka)', ylab='∆Lat (°N)')

polygon(x=c(results$age, rev(results$age)),y=c(results$DLat_lwr95, rev(results$DLat_upr95)),col = adjustcolor("grey67",alpha=0.3),border=NA)
polygon(x=c(results$age, rev(results$age)),y=c(results$DLat_upr68, rev(results$DLat_lwr68)),col = adjustcolor("grey47",alpha=0.3),border=NA)  
points(results$age, results$DLat, col=adjustcolor("grey17",alpha=0.9),type='l', pch=1, lwd=1.25)





