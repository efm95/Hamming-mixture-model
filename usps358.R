library(MBCbook)
library(AntMAN)
library(mcclust.ext)

source('../main/gibbs_sampler.R', echo=TRUE)

#=========================================================================================
# Loading data
#=========================================================================================

data('usps358')
X   = as.matrix(usps358[,-1])
cls = usps358[,1]

#=========================================================================================
# Data cleaning
#=========================================================================================

breaks  = c(-0.001,0.001,0.75,1,1.25,1.5,2)
m       = length(breaks)-1
classes = paste(breaks[1:m],(breaks[2:(m+1)]),sep="--|")
dati.m  = apply(X,2,cut,breaks = breaks,labels = classes)

#Dropping columns with less than m attributes
nclas<- rep(0,256)
for(i in 1:256){
  nclas[i]= length(table(dati.m[,i]))
}

for(i in 1:length(classes)){
  dati.m[dati.m==classes[i]]=i
}
dati.m = apply(dati.m, 2,as.numeric)
to_omit=which(nclas<m)
data.clean=as.data.frame(dati.m[,-to_omit])

#=========================================================================================
# Gibbs sampler
#=========================================================================================

gamma = AntMAN::AM_find_gamma_Pois(n=nrow(dati.omit),Lambda = 3,Kstar = 3) #gamma=0.1514657
u     = rep(3,ncol(data.clean))
v     = rep(0.5,ncol(data.clean))

set.seed(35)
gibbs = gibbs_mix_con(G=50000,
                      burnin = 10000,
                      data=data.clean,
                      u=u,
                      v=v,
                      Lambda = 3,
                      gam = gamma)


psm = comp.psm(gibbs$C[burnin:end,])

pred   = minVI(psm)$cl
pred_2 = minbinder(psm)$cl

