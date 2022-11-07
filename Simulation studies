library(mcclust)
library(mcclust.ext)

cat('Loading functions \n');cat('\n')
source('../main/data_generation.R', echo=TRUE)
source('../main/gibbs_sampler.R', echo=TRUE)
cat('Functions correctly loaded \n');cat('\n')

#=========================================================================================
# SIMULATION 1
#=========================================================================================

cat('Starting simulation set 1')

seed_set = 1:50

data_list_1 = list()
sim_1       = list()
psm_1       = list()
rand_VI_1   = list()
rand_bin_1  = list()


M.na.init = 1
k.init = 2

for (seed in seed_set) {
  
  cat('Simulation 1.',seed,'\n')
  
  data_list_1[[seed]] = ham_mix_gen(M=c(3,4,5),
                                    k = 3,
                                    p = 15,
                                    n=c(150,150,150),
                                    s=matrix(rep(0.2,3*15),ncol = 15),
                                    seed = seed)
  
  set.seed(seed)
  C.init = sample(2,(150*3),replace = T)
  sigma.init = matrix(runif((k.init+M.na.init)*15),nrow = (k.init+M.na.init), ncol=15)
  cent.init =  matrix(sample(2,(k.init+M.na.init)*15,replace = T),nrow=(k.init+M.na.init),ncol = 15)
  
  
  m = apply(data_list_1[[seed]]$data,2,function(x){length(table(x))})
  u=v=vector(length = 15)
  
  u[m==3]=5.00
  v[m==3]=0.25
  
  u[m==4]=4.50
  v[m==4]=0.25
  
  u[m==5]=4.25
  v[m==5]=0.25
  
  set.seed(seed)
  sim_1[[seed]] = gibbs_mix_con(G=10000,
                                burnin = 5000,
                                data=data_list_1[[seed]]$data,
                                C.init = C.init,
                                k.init=k.init,
                                M.na.init = M.na.init,
                                cent.init = cent.init,
                                sigma.init = sigma.init,
                                u=u,
                                v=v, 
                                Lambda = 3,
                                gam = 0.1905278)
  
  psm_1[[seed]] = comp.psm(sim_1[[seed]]$C)
  
  estim_1_VI_tmp  = mcclust.ext::minVI(psm_1[[seed]],method = 'all',cls.draw = sim_1[[seed]]$C)
  estim_1_VI = estim_1_VI_tmp$cl['best',]
  
  estim_1_bin_tmp  = mcclust.ext::minbinder.ext(psm_1[[seed]],method = 'all',cls.draw = sim_1[[seed]]$C)
  estim_1_bin = estim_1_bin_tmp$cl['best',]
  
  
  rand_VI_1[[seed]] = mcclust::arandi(estim_1_VI,data_list_1[[seed]]$groundTruth)
  rand_bin_1[[seed]] = mcclust::arandi(estim_1_bin,data_list_1[[seed]]$groundTruth)
}
cat('Saving outputs \n')
save.image("Server_output_1.RData")
