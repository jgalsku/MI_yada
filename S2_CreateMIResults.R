################################################################################
#                                                                              #
#                         Script 2: Mutual Information Calculation             #
#                                                                              #
################################################################################
#
# The following script imports DATA and VARIABLE INFORMATION. The script exports
# mutual information values into a .RDS file. DATA and VARIABLE INFORMATION must 
# be formatted correctly. Information to do so can be found at the following link
# https://rpubs.com/elainechu/mcp_vignette.
################################################################################

## Step 1: Package Dependencies
### See Script 1 to ensure all requisite packages are installed.

library(yada)
library(tidyverse)
library(mixtools)
library(doParallel)
library(foreach)



# registerDoParallel(detectCores())

set.seed(695432)

## Step 2: Minor logistics steps. Can be modified or avoided depending on use case
### Re-direct print statements to a text file for a permanent record of processing

# sink("results/output.txt")

setwd("C:/Users/jgalimanyskupham/Documents/2025_MI-Chris")


### Check that a results folder exists in the working directory
if(! ("results" %in% dir()) ) {
  stop("There is no 'results' folder in the working directory")
}

### The data directory is /results (user-defined):
data_dir <- "results"

### The "ID" that uniquely identifies this analysis (user-defined):
analysis_name <- 'MI'

## Step 3: Import DATA and VARIABLE INFORMATION

var_info <-  load_var_info('results/var_info_Stress_test.csv')
data_file <- 'results/data.csv'
dat <- reformat_mcp_data(data_file_path=data_file,
                         var_info=var_info)

## Step 4: Create Problem Files
### Note, user dictates type of problem file (ordinal vs. continuous or both)

main_problem <- generate_main_problem(data_file=dat, 
                                      var_info=var_info)

save_problem(data_dir=data_dir,
             analysis_name=analysis_name,
             problem=main_problem,
             is_folds=F)

ord_prob_list <- build_univariate_ord_problems(data_dir,
                                                                  analysis_name,
                                                                  add_folds=F)

# cont_prob_list <- build_univariate_cont_problems(data_dir,
#                                                  analysis_name,
#                                                  add_folds=F)

## Step 5: Optimize Univariate Models
### Continuous and Ordinal variables are optimized separately

base_seed <- 264528
set.seed(base_seed)
seed_vect <- sample.int(1000000, length(ord_prob_list), replace=F)

#### Ordinal Models

ord_success <-
  foreach::foreach(i=1:length(ord_prob_list), .combine=cbind) %do% {
    yada::solve_ord_problem(data_dir,
                            analysis_name,
                            ord_prob_list[[i]],
                            anneal_seed=seed_vect[i])
  }

#### Continuous Models

# cont_success <-
#   foreach::foreach(i=1:length(cont_prob_list), .combine=cbind) %dopar% {
#     yada::solve_cont_problem(data_dir, analysis_name, cont_prob_list[[i]])
#   }

## Step 6: Solve for X
### This step is necessary to solve for mututal information over X (age)

#### offset for the weibull fit

weib_offset <- 0.002

#### fit mixed weibull to age data
weib_fit <- weibullRMM_SEM(main_problem$x + weib_offset,
                                     k=3,
                                     maxit=2000)
#### save the results of the fit
theta_x <- list(fit_type='offset_weib_mix',
                fit=weib_fit,
                weib_offset=weib_offset)

#### save rds to used later
saveRDS(theta_x,build_file_path(data_dir, analysis_name, "solutionx"))

## Step 7: Model Selection + MI Calculation
### Here we use AIC. It is possible to use K-Fold CV or other technique.
### The model selection and MI calculation is completed in one step. 
### The continuous variables utilize a secondary KL Divergence function.

####  Ord MI Function
calc_univ_ord_mi <- function(th_v, mod_spec, th_x, x0, xcalc) {
  # Calculate the mutual information for a univariate ordinal model given the
  # baseline age, x0, and a vector at which to calculate the prior and
  # posterior densitites, xcalc.
  fprior <- calc_x_density(xcalc,th_x)
  
  M <- mod_spec$M
  pv     <- rep(NA, M+1)
  kl_div <- rep(NA, M+1)
  for (m in 0:M) {
    x_post_obj  <- calc_x_posterior(m,
                                    th_x,
                                    model,
                                    xcalc)
    kl_div[m+1] <- calc_kl_div(x_post_obj, th_x)
    pv[m+1] <- calc_q(x0, th_v, m, mod_spec)
  }
  return(sum(pv*kl_div))
}

#### Cont Kl function
calc_cont_kl_div_vect <- function(th_w, mod_spec, xcalc, wcalc, th_x) {
  # Calculat a vector of KL divergences for each entry in wcalc
  N1 <- length(wcalc)
  N2 <- length(xcalc)
  kl_div <- rep(NA, N1)
  
  # Calculate the prior
  fprior0 <- calc_x_density(xcalc,th_x)
  
  # Use a single for loop to calcualte the likelihood matrix (rather than
  # unwrapping matrices to make a single call to calc_neg_log_lik_cont).
  dx <- xcalc[2] - xcalc[1]
  for (n1 in 1:N1) {
    eta_vect <- calc_neg_log_lik_vect_cont(th_w, xcalc, rep(wcalc[n1], N2), 
                                           mod_spec)
    lik_vect <- exp(-eta_vect)
    fprior <- fprior0
    fpost <- fprior * lik_vect
    fpost <- fpost / dx / sum(fpost)
    ind_bad <- is.na(fprior) | is.na(fpost)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    ind_bad <- !is.finite(fprior) | !is.finite(fpost)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    ind_bad <- (fprior == 0) | (fpost == 0)
    fprior <- fprior[!ind_bad]
    fpost <- fpost[!ind_bad]
    kl_div[n1] <- sum(fpost * log2(fpost/fprior)) * dx
  }
  
  return(kl_div)
}

#### Cont MI function
calc_univ_cont_mi <- function(th_w, mod_spec, x0, wcalc, kl_div) {
  # Calculate the mutual information for a univariate cont. model given the
  # baseline age, x0, and a vector of KL divergences with length(wcalc) that
  # was output by calc_cont_kl_div_vect.
  
  N <- length(wcalc)
  dw <- wcalc[2] - wcalc[1]
  pw <- calc_neg_log_lik_vect_cont(th_w, rep(x0, N), wcalc, mod_spec)
  pw <- exp(-pw)
  pw <- pw / sum(pw) / dw
  return(dw*sum(pw*kl_div))
}

#### Ordinal MI Calculation

##### Pull out ordinal first
j_ord <- main_problem$mod_spec$J

##### vector to calculate prior and posterior densities
xcalc = seq(0,105,by=0.01)

##### baseline age
x0 = seq(0,105,by=.1)

##### create and empty list to store ordinal results (j = 44)
MI_ord <- setNames(vector("list", j_ord), main_problem$var_names[1:j_ord])

# Loop over each ord var, select the best model,print best model,do MI based on
#best model and store results in MI_ord


library(foreach)
library(doParallel)

ncores <- parallel::detectCores() -1
cl <- makeCluster(ncores)

registerDoParallel(cl)

clusterEvalQ(cl, {
  library(foreach)
  library(doParallel)
  library(yada)
  })


for (i in 1:j_ord) {
  
  var_name <- main_problem$var_names[i]
  
  print(paste0("AIC model selection for ", var_name))
  
  aic_output <- build_aic_output(data_dir, analysis_name, var_name,
                                       format_df=T, save_file=TRUE)
  
  c <- filter(aic_output, rank == 1)$model
  
  md <- paste0(data_dir,"/solutiony_",analysis_name,"_","ord","_","j_",i,
               "_",var_name,"_",c,".rds")
  
  file <- readRDS(md)
  
  th_v = file[[1]]
  
  mod_spec = file[[2]]
  
  print(c)
  
  MI_ord[[i]] <- foreach(n=1:length(x0), .multicombine = T,
                         .packages = c("yada")) %dopar%{
                           calc_univ_ord_mi(th_v, mod_spec, theta_x, x0[n],xcalc, file)
                         }
  
}

#### Continuous MI Calculation

##### Pull out # cont variables for loop below
k_cont <- main_problem$mod_spec$K

##### vector to calculate density over w
wcalc <- seq(0, 500, by=1)

##### create an empty list to store continuous results (k = 18)
MI_cont <- setNames(vector("list", k_cont), main_problem$var_names[45:62])

##### Loop over each cont var, select the best model with aic, print best model,  
#do KL div,do MI based on best model and kl results and store results in MI_cont
for (i in 1:k_cont) {
  
  var_name <- main_problem$var_names[45:62][i]
  
  print(paste0("AIC model selection for ", var_name))
  
  aic_output <- build_aic_output(data_dir, analysis_name, var_name,
                                       format_df=T, save_file=TRUE)
  
  c <- filter(aic_output, rank == 1)$model
  
  md<-paste0(data_dir,"/solutiony_",analysis_name,"_","cont","_","k_",i,
             "_",var_name,"_",c,".rds")
  
  file <- readRDS(md)
  
  th_w = file[[1]]
  
  mod_spec = file[[2]]
  
  print(c)
  
  kl <- calc_cont_kl_div_vect(th_w = th_w, mod_spec = mod_spec, xcalc = xcalc, 
                              wcalc = wcalc, th_x = theta_x)
  
  MI_cont[[i]] <- foreach(n=1:length(x0), .multicombine = T,
                          .packages = c("yada")) %dopar%{
                            calc_univ_cont_mi(th_w = th_w, mod_spec = mod_spec, x0 = x0[n], 
                                              wcalc = wcalc, kl_div = kl)
                          }
  
}

## Step 8: Save Results
### To be used in Script 3: Visualizations

write_rds(MI_ord, "results/MI_ord_allVars.rds")
write_rds(MI_cont,"results/MI_cont_allvars.rds")

## Step 9: End Parallel and Close Results

### Close all clusters used for paralle processing
stopImplicitCluster()

### End the re-directing of print statements to file
sink()

####################################END#########################################
