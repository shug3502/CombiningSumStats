#run death process model to get sample from true posterior and compare with abc 
  library(rstan)
  library(mvtnorm)
  library(dplyr)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  set.seed(345)
##################################################  
  # exp_decay <- function(N=10,k=0.1,times){
  #   ttt = c(0)
  #   output = rep(NA,length(times))
  #   for (j in seq_along(times)){
  #     while ((ttt < times[j]) && (N>0)) {
  #       tau = rexp(1,rate=k*N)
  #       ttt = ttt + tau
  #       N = N - 1 #only reaction that can occur
  #     }       
  #     output[j] = N
  #   }
  #   return(output)
  # }
##################################################  
   identifier = 'death_processV213' #run identifier
   run_mcmc <- TRUE
   show_diagnostic_plots = TRUE
   parametersToPlot = c('k','log_k','sigma','log_sigma')
   n <- 7 #ignore first time point, not informative, just initial condition
   times <- seq(from=0,to=20,length=(n+2))
   times = times[-(n+2)]
   observations <- read.csv('../code_for_figures/death_process_v213_data.csv',header=FALSE)
   if (ncol(observations)!=n+2) warning('check using correct data, and n value, number of columns not n+2')
#   observations <- observations[1:2,]
   # observations <- exp_decay(N=N,k=0.1,times[-1])   #c(6,4,3,2,1,1,1,1,1,1,1,rep(0,20))
   N <- observations[1,1]
   y <- observations[,c(-1,-(n+2))]
   z <- observations[,(n+2)]
   repeats <- nrow(observations)
   
     ##########################
  if (run_mcmc) {
    estimates <- stan(file = 'death_process_posterior.stan',
                      data = list (
                      N = N, n = n, repeats = repeats, y = y, z=z, t = times[-1]
                      ),
                      seed = 42,
                      chains = 4,
                      warmup = 1000,
                      iter = 2000
    )
    
    tryCatch({
      saveRDS(estimates, file = paste('death_process_',identifier,'.rds',sep='')) # to load use readRDS("fit.rds")
    }, warning = function(war) {  
      # warning handler picks up where error was generated
      print(paste("WARNING:  ",war))
    }, error = function(err) {  
      # error handler picks up where error was generated
      print(paste("NOT SAVED ERROR:  ",err))
    }) # END tryCatch
  } else {
    estimates = readRDS(paste('death_process',identifier,'.rds',sep=''))
  }
  
 print(estimates, pars = parametersToPlot)
  
  #######################
  #visualisation
  ###############################
  # #cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
  # pairs(estimates, pars = parametersToPlot)
  # #dev.off()
  # #ggsave(paste('plots/pairs',identifier, '.eps',sep=''),device=cairo_ps)
  # 
  # #look at posterior predictive distn
 
  if (show_diagnostic_plots) {
    # source('hist_treedepth.R')
    # hist_treedepth(estimates)
    # source('mcmcDensity.R')
    # mcmcDensity(estimates, parametersToPlot, byChain = TRUE)
    # ggsave(paste('plots/denisty',identifier, '.eps',sep=''),device=cairo_ps)
    
    library(bayesplot)
    draws <- as.array(estimates, pars=parametersToPlot)
    mcmc_trace(draws)
#    ggsave(paste('plots/trace',identifier, '.eps',sep=''),device=cairo_ps)
    mcmc_intervals(draws,pars=parametersToPlot)
#    ggsave(paste('plots/intervals',identifier, '.eps',sep=''),device=cairo_ps)
    color_scheme_set("purple")
    mcmc_areas(draws,pars=c('log_k','log_sigma'))
    ggsave(paste('areas',identifier, '.eps',sep=''),device=cairo_ps)
    color_scheme_set("brightblue")
#    cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
#    pairs(estimates, pars = parametersToPlot)
#    dev.off()
    # mcmc_scatter(draws,pars=parametersToPlot)
    # ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)
  }
  
  ######################
  #cos saving wasn't working
 require(tidyr)
 require(ggplot2)
 pred <- as.data.frame(estimates, pars = 'y_sim') %>%
   gather(factor_key = TRUE) %>%
   group_by(key) %>%
   summarize(lb = quantile(value, probs = 0.05),
             median = quantile(value, probs = 0.5),
             ub = quantile(value, probs = 0.95))
 xdata <- data.frame(times = times[-1], y = t(y))
 xdata <- xdata %>% gather(sim,molecules,-times)
 # pred <- pred %>% bind_cols(xdata) %>%
 #   mutate(split = if_else(time %in% ts_test,'train','test')) 
 pred["times"] = times[-1]
 p1 <- ggplot(pred, aes(x = times, y = median))
 p1 <- p1 + geom_line() +
   geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
   labs(x = "Time", y = "Molecules") +
   theme(text = element_text(size = 12), axis.text = element_text(size = 12),
         legend.position = "none", strip.text = element_text(size = 8))
   p1 <- p1 + geom_line(data = xdata, aes(x = times, y = molecules, colour=sim))
  print(p1)
  ggsave(paste(identifier, '_pred.eps',sep=''),device=cairo_ps)
###################################################
  #write output to read into matlab
  e <- rstan::extract(estimates,pars=c('log_k','log_sigma'),permuted=TRUE)
  write.csv(e,file='death_process_stan_estimates.csv',row.names=FALSE)
