library(AntMAN)
library(mcclust.ext)

source('main/gibbs_sampler.R', echo=TRUE)

#=========================================================================================
# Loading data
#=========================================================================================

zoo=read.table("zoo.data",h=F,sep=",")
nam = zoo$V1
groundTruth = zoo$V18
classes = factor(groundTruth,labels=c("mammals", "birds", "reptiles", "fish", 
                                      "amphibians", "insects", "mollusks"))
names(groundTruth)<-classes

#=========================================================================================
# Data cleaning
#=========================================================================================

zoo = as.matrix(zoo[,-c(1,18)]+1)

attriList = function(mat){
  ans = list()
  for (i in 1:ncol(mat)) {
    ans[[i]]=sort(unique(mat[,i]))
  }
  return(ans)
}

attriList(zoo)
zoo[,13] = ifelse(zoo[,13]==3,2,
                  ifelse(zoo[,13]==5,3,
                         ifelse(zoo[,13]==6,4,
                                ifelse(zoo[,13]==7,5,
                                       ifelse(zoo[,13]==9,6,
                                              1)))))


#=========================================================================================
# Gibbs sampler
#=========================================================================================

Kstar  = 7
Lambda = 7
gam    = AntMAN::AM_find_gamma_Pois(n=nrow(zoo),Lambda=Lambda,Kstar=Kstar)

u = c(rep(6,12),3,rep(6,3))
v = c(rep(0.25,12),0.5,rep(0.25,3))

sim_zoo = gibbs_mix_con(G=25000,
                        burnin = 5000,
                        data=zoo,
                        u=u,v=v,
                        Lambda = Lambda,
                        gam = gam)

psm = comp.psm(sim_zoo$C)
Bin = minbinder(psm)
VI  = minVI(psm)

table(Bin$cl)
table(VI$cl)

arandi(Bin$cl,groundTruth)
arandi(VI$cl,groundTruth) 
