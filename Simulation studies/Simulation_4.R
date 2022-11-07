library(mcclust)
library(mcclust.ext)

cat('Loading functions \n');cat('\n')
source('../main/data_generation.R', echo=TRUE)
source('../main/gibbs_sampler.R', echo=TRUE)
cat('Functions correctly loaded \n');cat('\n')

seed_set = 1:50

#=========================================================================================
# SIMULATION 4
#=========================================================================================

cat('Starting simulation set 4')

data_list_4 = list()
sim_4       = list()
psm_4       = list()
rand_VI_4   = list()
rand_bin_4  = list()

M.na.init = 1
k.init = 2

for (seed in seed_set) {
  
  cat('Simulation 4.',seed,'\n')
  
  p=15
  
  data_list_4[[seed]] = ham_mix_gen(M=c(3,4,5),
                                    k = 3,
                                    p = 15,
                                    n=c(25,25,25),
                                    s=matrix(rep(0.5,3*p),ncol = p),
                                    seed = seed)
  
  set.seed(seed)
  C.init = sample(2,(25*3),replace = T)
  sigma.init = matrix(runif((k.init+M.na.init)*p),nrow = (k.init+M.na.init), ncol=p)
  cent.init =  matrix(sample(2,(k.init+M.na.init)*p,replace = T),nrow=(k.init+M.na.init),ncol = p)
  
  m = apply(data_list_4[[seed]]$data,2,function(x){length(table(x))})
  u=v=vector(length = p)
  
  u[m==3]=5.00
  v[m==3]=0.25
  
  u[m==4]=4.50
  v[m==4]=0.25
  
  u[m==5]=4.25
  v[m==5]=0.25
  
  sim_4[[seed]] = gibbs_mix_con(G=10000,
                                burnin=5000,
                                data=data_list_4[[seed]]$data,
                                u=u,
                                v=v,
                                C.init = C.init,
                                k.init = k.init,
                                M.na.init = M.na.init,
                                cent.init = cent.init,
                                sigma.init = sigma.init,
                                Lambda = 3,
                                gam = 0.3028313)
  
  psm_4[[seed]] = comp.psm(sim_4[[seed]]$C)
  
  estim_4_VI_tmp  = mcclust.ext::minVI(psm_4[[seed]],method = 'all',cls.draw = sim_4[[seed]]$C)
  estim_4_VI = estim_4_VI_tmp$cl['best',]
  
  estim_4_bin_tmp  = mcclust.ext::minbinder.ext(psm_4[[seed]],method = 'all',cls.draw = sim_4[[seed]]$C)
  estim_4_bin = estim_4_bin_tmp$cl['best',]
  
  rand_VI_4[[seed]] = mcclust::arandi(estim_4_VI,data_list_4[[seed]]$groundTruth)
  rand_bin_4[[seed]] = mcclust::arandi(estim_4_bin,data_list_4[[seed]]$groundTruth)
}

#HD = list()
#k_modes = list()
#for (seed in seed_set) {
#  cat('Iteration', seed,'\n');cat('\n')
#  HD[[seed]]= CategorialCluster(data_list_4[[seed]]$data)
#  k_modes[[seed]]=klaR::kmodes(data=data_list_4[[seed]]$data,
#                               modes = 3)
#}

cat('Saving outputs \n')
save.image("Server_output_4.RData")

#HD_rand = c()
#kmodes_rand = c()
#for (seed in seed_set) {
#  HD_rand[seed]=arandi(HD[[seed]][[1]],data_list_4[[seed]]$groundTruth)
#  kmodes_rand = arandi(k_modes[[seed]]$cluster,data_list_4[[seed]]$groundTruth)
#}

#boxplot(cbind(HD_rand,kmodes_rand,unlist(rand_VI_4),unlist(rand_bin_4)),ylim = c(0,1))
