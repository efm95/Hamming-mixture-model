library(mcclust)
library(mcclust.ext)

cat('Loading functions \n');cat('\n')
source('../main/data_generation.R', echo=TRUE)
source('../main/gibbs_sampler.R', echo=TRUE)
cat('Functions correctly loaded \n');cat('\n')

seed_set = 1:50

#=========================================================================================
# SIMULATION 3
#=========================================================================================

cat('Starting simulation set 3')

data_list_3 = list()
sim_3       = list()
psm_3       = list()
rand_VI_3   = list()
rand_bin_3  = list()

M.na.init = 1
k.init = 2

for (seed in seed_set) {
  
  cat('Simulation 3.',seed,'\n')
  
  p=10
  
  data_list_3[[seed]] = ham_mix_gen(M=c(3,4,5),
                                    k = 4,
                                    p = 10,
                                    n=c(50,50,50,50),
                                    s=matrix(rep(0.7,4*p),ncol = p),
                                    seed = seed)
  
  set.seed(seed)
  C.init = sample(2,(50*4),replace = T)
  sigma.init = matrix(runif((k.init+M.na.init)*p),nrow = (k.init+M.na.init), ncol=p)
  cent.init =  matrix(sample(2,(k.init+M.na.init)*p,replace = T),nrow=(k.init+M.na.init),ncol = p)
  
  m = apply(data_list_3[[seed]]$data,2,function(x){length(table(x))})
  u=v=vector(length = p)
  
  u[m==3]=5.00
  v[m==3]=0.25
  
  u[m==4]=4.50
  v[m==4]=0.25
  
  u[m==5]=4.25
  v[m==5]=0.25
  
  sim_3[[seed]] = gibbs_mix_con(G=10000,
                                burnin = 5000,
                                data=data_list_3[[seed]]$data,
                                u=u,
                                v=v,
                                C.init = C.init,
                                k.init = k.init,
                                M.na.init = M.na.init,
                                cent.init = cent.init,
                                sigma.init = sigma.init,
                                Lambda = 4,
                                gam = 0.3028313)
  
  psm_3[[seed]] = comp.psm(sim_3[[seed]]$C)
  
  estim_3_VI_tmp  = mcclust.ext::minVI(psm_3[[seed]],method = 'all',cls.draw = sim_3[[seed]]$C)
  estim_3_VI = estim_3_VI_tmp$cl['best',]
  
  estim_3_bin_tmp  = mcclust.ext::minbinder.ext(psm_3[[seed]],method = 'all',cls.draw = sim_3[[seed]]$C)
  estim_3_bin = estim_3_bin_tmp$cl['best',]
  
  rand_VI_3[[seed]] = mcclust::arandi(estim_3_VI,data_list_3[[seed]]$groundTruth)
  rand_bin_3[[seed]] = mcclust::arandi(estim_3_bin,data_list_3[[seed]]$groundTruth)
}

cat('Saving outputs \n')
save.image("Server_output_3.RData")
