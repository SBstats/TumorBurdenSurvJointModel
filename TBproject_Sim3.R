# Load R packages
library(survival)
library(mvtnorm)
library(gridExtra)
#remove.packages(c("StanHeaders", "rstan"))
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
library(data.table)
library(MCMCpack)
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load R code
# ------ Define of constants (adjust to fit the data generating scenario) ------
true_gamma_A = -0.75
true_beta_A = 0

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J option

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)


txt.title = paste0("simu_results/Results_TBproj_Sim3_result.txt")
if (run_ID == 1) {
  df = data.frame(matrix(ncol = 25, nrow = 0))
  df_col_names = c("run_ID",
                   "gamma_X1_post_mean", "gamma_A_post_mean",
                   "beta_A_post_mean", "sig_err_post_mean", "beta_X1_post_mean", 
                   "b_mu1_post_mean", "b_mu2_post_mean", "b_sd1_post_mean", "b_sd2_post_mean",
                   "TB_point_est", "TB_point_est_BB", "TB_se", "TB_test_stat", "TB_pval",
                   "Surv_point_est", "Surv_point_est_BB", "Surv_se", "Surv_test_stat", "Surv_pval",
                   "Tot_point_est", "Tot_point_est_BB", "Tot_se","Tot_test_stat", "Tot_pval")
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------LOAD CODE FOR ONE REPLICATION HERE----------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

n <- 500  # Number of individuals

# Generate random effects

generate_random_effects = function(k,           #number of random effects
                                   n,           # number of units
                                   mean_vec,      # means for each random effect, length = k
                                   corr_mat,     # correlation matrix, size k*k
                                   sd_vec         # standard deviation vector for b_{ik}, length = k
){
  cholesky_decom = chol(corr_mat) #spits out  L^T, not L
  cat("cholesky decomposition matrix of the CORRELATION (not covariance) matrix, L: \n")
  print(t(cholesky_decom))
  LL_transpose = t(cholesky_decom) %*% cholesky_decom  
  cat("LL^T: \n")
  print(LL_transpose)
  #sd_diag_mat = diag(sd_vec)
  Lambda_mat =     diag(sd_vec) %*% t(cholesky_decom)    # like a square root of covariance matrix
  z= matrix(nrow = k, ncol = n)
  for(i in 1:k){
    z[i,] = rnorm(ncol(z),0,1)
  }
  b_mat = Lambda_mat %*% z + mean_vec    # b_i = L^{T} * z_i + \mu
  print(cbind(c("Mean of bi1", "Mean of bi2"), as.numeric(format(round(rowMeans(b_mat), 4), nsmall = 4))))
  print(cbind(c("SD of bi1", "SD of bi2" ), c(sd(b_mat[1,]),sd(b_mat[2,]))))
  return(b_mat)
}




generate_random_effects_indep_norm = function(k,           #number of random effects
                                              n,           # number of units
                                              mean_vec,      # means for each random effect, length = k
                                              sd_vec         # standard deviation vector for b_{ik}, length = k
){
  
  b_mat= matrix(nrow = k, ncol = n)
  for(i in 1:k){
    b_mat[i,] = rnorm(ncol(b_mat),mean_vec[i],sd_vec[i])
  }
  print(cbind(c("Mean of bi1", "Mean of bi2"), as.numeric(format(round(rowMeans(b_mat), 4), nsmall = 4))))
  print(cbind(c("SD of bi1", "SD of bi2" ), c(sd(b_mat[1,]),sd(b_mat[2,]))))
  return(b_mat)
}




### Generate enrollment times, visit times and end of  study times

# Simulate Enrollment Time
#
# Simulate enrollment time by total time

stb_tl_simu_enroll_by_dur <- function(n_pt, n_pt_tot, pt_dur_mth) {
  data.frame(mth_enroll = runif(n_pt,
                                0,
                                pt_dur_mth))
}

# Simulate Enrollment Time
#
# Simulate enrollment time by rate

stb_tl_simu_enroll_by_rate <- function(n_pt, n_pt_tot, pt_per_mth) {
  pt_dur_mth <- n_pt_tot / pt_per_mth
  data.frame(mth_enroll = runif(n_pt,
                                0,
                                pt_dur_mth))
}

# Simulate Enrollment Time
#
# Simulate enrollment time by center

stb_tl_simu_enroll_by_center <- function(n_pt,
                                         n_pt_tot,
                                         n_center,
                                         pt_per_center_per_mth,
                                         center_per_mth) {
  
  center_dur_mth <- n_center / center_per_mth
  pt_dur_mth     <- n_pt_tot / pt_per_center_per_mth
  
  en_center      <- runif(n_center, 0, center_dur_mth)
  en_center      <- sort(en_center)
  
  rst <- NULL
  for (i in 1:n_center) {
    cur_center <- en_center[i]
    cur_en_pt  <- cur_center + runif(n_pt, 0, pt_dur_mth)
    rst        <- rbind(rst, cbind(center            = i,
                                   mth_enroll_center = cur_center,
                                   mth_enroll        = cur_en_pt))
  }
  
  data.frame(rst) %>%
    arrange(mth_enroll) %>%
    slice_head(n = n_pt)
}

# Simulate Enrollment Time
#
# mth_fix_fu: fixed follow-up months for all patients
# mth_min_fu: minimum follow-up months for all patients

stb_tl_simu_enroll <- function(n_pt_arm    = 100,
                               n_pt_tot    = n_pt_arm * 2,
                               param_enroll  = list(type  = "by_rate",
                                                    pt_per_mth = 50),
                               mth_min_fu  = NULL,
                               mth_fix_fu  = NULL,
                               date_bos    = "2022-01-01",
                               mth_to_days = 30.4,
                               ...) {
  
  type                <- param_enroll$type
  param_enroll$type     <- NULL
  param_enroll$n_pt     <- n_pt_arm
  param_enroll$n_pt_tot <- n_pt_tot
  f_enroll            <- switch(type,
                                by_duration = stb_tl_simu_enroll_by_dur,
                                by_rate     = stb_tl_simu_enroll_by_rate,
                                by_center   = stb_tl_simu_enroll_by_center)
  
  rst <- do.call(f_enroll, param_enroll) %>%
    mutate(day_enroll = mth_to_days * mth_enroll) %>%
    arrange(day_enroll) %>%
    mutate(sid = row_number())
  
  ## set up dates in addition to days
  if (!is.null(date_bos)) {
    rst$date_bos    <- as.Date(date_bos)
    rst$date_enroll <- as.Date(date_bos) + rst$day_enroll
  }
  
  ## set up end of study time by minimum or fixed follow up
  day_chk <- NULL
  if (!is.null(mth_min_fu)) {
    day_chk <- max(rst$day_enroll)
    min_fu  <- mth_min_fu
  } else if (!is.null(mth_fix_fu)) {
    day_chk <- rst$day_enroll
    min_fu  <- mth_fix_fu
  }
  
  if (!is.null(day_chk)) {
    day_eos      <- day_chk + min_fu * mth_to_days
    day_eos      <- floor(day_eos)
    #rst$day_eos  <- day_eos - rst$day_enroll
    #rst$date_eos <- rst$date_bos + day_eos
  }
  
  ## return
  rst
}

# Simulate Enrollment Time
#
# Simulate Enrollment Time for all Arms

stb_tl_simu_enroll_arms <- function(n_by_arm,                  # this is a vector with size = number of arms
                                    ...,
                                    seed = NULL) {
  
  if (!is.null(seed))
    old_seed <- set.seed(seed)
  
  n_arm <- length(n_by_arm)
  rst   <- NULL
  for (i in seq_len(n_arm)) {
    cur_arm     <- stb_tl_simu_enroll(n_by_arm[i],
                                      n_pt_tot = sum(n_by_arm),
                                      ...)
    cur_arm$arm <- i - 1
    rst         <- rbind(rst, cur_arm)
  }
  
  ## reset
  if (!is.null(seed))
    set.seed(old_seed)
  
  ## return
  rst
}



generate_survdata_noRE = function(n, 
                                  baseline_data,
                                  Enrollment_times_df,
                                  gamma_X1,
                                  gamma_A,                     # treatment effect
                                  rate_param_cen
                                  
){
  ID <- baseline_data$ID
  Trt <- baseline_data$A
  X1 <- baseline_data$X1
  rate_param_vec_death <- c()
  time_Death_joint <- c()
  time_Cen_joint <- c()
  #rate_param_vec_cen <- c()
  for(id in ID){
    rate_param_vec_death[id] = exp( (gamma_X1* X1[id]) + (gamma_A *Trt[id]) )
    #rate_param_vec_death[id] = exp( gamma_0 + (gamma_A *Trt[id]) )
    
    time_Death_joint[id] <- rexp(1, rate = rate_param_vec_death[id])       # Death times (exponentially distributed)
    time_Cen_joint[id] <- rexp(1, rate = rate_param_cen)                  # Censored times (exponentially distributed)
  }
  
  time_Obs_joint <- pmin(time_Death_joint, time_Cen_joint) 
  status_joint <- ifelse(time_Death_joint <= time_Cen_joint,1,0)                # Event status (0 = censored, 1 = event)
  
  
  #time_Obs_joint <- time_Death_joint                      # Censored times (exponentially distributed)
  #status_joint <- rep(1,n)                # Event status (0 = censored, 1 = event)
  
  
  df_surv_joint <- data.frame(
    ID = ID,
    Trt = Trt,
    X1 = X1,
    Date_enroll = Enrollment_times_df$date_enroll,
    time_Death_joint = time_Death_joint,
    Delta = status_joint,
    ObsTime = time_Obs_joint,
    Date_obsTime = Enrollment_times_df$date_enroll + time_Obs_joint
  )
  
  #print(rate_param_vec_death)
  
  return(df_surv_joint)
  
}




get_endOfStudy_basedOnEvents = function(surv_data,
                                        target_event_num){
  num_death = sum(surv_data$Delta)
  if(target_event_num > num_death){
    warning("The total number of events is less than the target number.")
    return(NA)
  }
  surv_data_death_only = surv_data[surv_data$Delta==1,]
  surv_data_death_only = surv_data_death_only %>%
    arrange(Date_obsTime)
  date_temp = surv_data_death_only$Date_obsTime[target_event_num]
  return_df = surv_data %>%
    mutate(Date_eos = date_temp) %>%
    filter(Date_enroll <= Date_eos)
  return_df = return_df %>%
    mutate(Delta = if_else(Date_obsTime <= Date_eos,
                           Delta,
                           0),
           Date_obsTime = if_else(Date_obsTime <= Date_eos,
                                  Date_obsTime,
                                  Date_eos),
           ObsTime = as.numeric(Date_obsTime - Date_enroll),
           Day_eos = as.numeric(Date_eos - Date_enroll)
    )
  return_df = return_df %>%
    arrange(ID)
  
  return(return_df)  
}






generate_baseline_TB_data = function(num_subjects,baseline_data){
  ID = 1:num_subjects
  A = baseline_data$A
  Y0 = rnorm(num_subjects,15,0.5) 
  return(data.frame(ID = ID, A = A, Y0 = Y0))
}



generate_longitudinal_data = function(n,                            # number of individuals
                                      max_visits,                   #maximum number of visits per individual
                                      baseline_data,                #dataframe of baseline data
                                      visits_days_df,
                                      surv_time,                    #observed death time 
                                      true_death_time,              
                                      b_mat,                        #matrix if random effects: dim = k*n, k = # random effects
                                      beta_A,                       #treatment effect
                                      sd_error,
                                      beta_X1
){
  
  
  longitudinal <- vector("list", n)  # List to store individual longitudinal data
  all_visit_days_df = data.frame(visits_days_df[, 5:ncol(visits_days_df)])
  for (i in 1:n) {
    
    #visit <- 1:max_visits
    visit_days = all_visit_days_df[i,][all_visit_days_df[i,] < surv_time[i]]   # Truncate visit days at survival time (days) for each individual
    
    #visit = 1:length(visit_days)
    
    # Repeat baseline variables for each visit
    baseline_rep <- data.frame(
      ID = rep(baseline_data$ID[i], length(visit_days)),
      visit_days = visit_days,
      A = rep(baseline_data$A[i], length(visit_days)),
      X1 = rep(baseline_data$X1[i], length(visit_days)),
      #X2 = rep(baseline_data$X2[i], length(visit_days)),
      True_death_time = rep(true_death_time[i],length(visit_days)),
      Surv_time = rep(surv_time[i], length(visit_days))
    )
    
    
    #psi_ij = as.numeric(all_visit_days_df[i,]/surv_time[i])
    #psi_ij = visit_days/surv_time[i]
    psi_ij = visit_days/true_death_time[i]
    #print(psi_ij)
    
    #....................Generate TB data.......................\\
    epsilon <- rnorm(length(visit_days),0,sd_error)
    #fixed_effects_long <-  baseline_data$X1[i] * beta_X1 + baseline_data$X2[i] * beta_X2 + baseline_data$A[i] * beta_A
    #TB <- fixed_effects_long +  b_mat[1,i] * visit + b_mat[2,i] * (visit)^2 + epsilon
    #fixed_effects_long <-  beta_0 + baseline_data$X1[i] * beta_X1 +  baseline_data$A[i] * beta_A
    fixed_effects_long <-   baseline_data$X1[i] * beta_X1 +  baseline_data$A[i] * beta_A
    TB <- fixed_effects_long + psi_ij * b_mat[1,i]   + (psi_ij)^2 * b_mat[2,i]   + epsilon
    baseline_rep$TB <- as.vector(TB)
    longitudinal[[i]]<-baseline_rep
  }
  
  # Combine baseline and longitudinal data into a single dataframe
  df_long <- data.frame(do.call(rbind, longitudinal))
  df_long <- df_long %>%
    group_by(ID) %>%
    mutate(visit = row_number())
  
  return(data.frame(df_long))
  
}



library(dplyr)

combine_longitudinal_data = function(tb_baseline,df_long){
  new_rows <- data.frame(
    ID = tb_baseline$ID, 
    A = tb_baseline$A,
    visit = rep(0, times = nrow(tb_baseline)), 
    TB = tb_baseline$Y0 
  )
  
  # Bind the new rows to the original dataframe
  return_df <- bind_rows(df_long, new_rows)
  
  # Sort the dataframe by ID and visit
  return_df <- return_df %>%
    arrange(ID, visit)
  
  
  return(return_df)
}



generate_longitudinal_data_change_from_baseline = function(data_long_no_baseTB, base_TB){
  temp_base_TB_df = data.frame(ID = base_TB$ID, Y0 = base_TB$Y0)
  return_df = merge(data_long_no_baseTB, temp_base_TB_df, by = "ID", all.x = TRUE)
  return_df$TBfrombase = (return_df$TB - return_df$Y0)/return_df$Y0
  return(return_df)
}



plot_longitudinal_data = function(data_long){
  treatment <- ifelse(data_long$A == 1, "Drug","Control")
  data_long$treatment = treatment
  
  
  spider_plot <- ggplot(data_long, aes(x=visit, y=TB, group=ID)) +
    theme_bw(base_size=14) +
    theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(face="bold")) +
    theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    theme(plot.title = element_text(size=18, hjust=0.5)) +
    labs(list(x = "Visit", y = "Pseudo TB Data"))#+geom_smooth()
  ## Now plot
  spider_plot <- spider_plot + geom_line(aes(color=treatment)) +
    geom_point(aes(color=treatment), show.legend=FALSE) +
    scale_colour_discrete(name="Treatment", labels=c("Control", "Drug")) +
    scale_shape_manual(name = "cstatus", values = c("0"=3, "1"=16)) +
    coord_cartesian(xlim=c(1, max(data_long$visit)))+
    facet_grid(~treatment)
  
  spider_plot
}



plot_longitudinal_data_change_from_base = function(data_long){
  treatment <- ifelse(data_long$A == 1, "Drug","Control")
  data_long$treatment = treatment
  
  
  spider_plot <- ggplot(data_long, aes(x=visit, y=TBfrombase, group=ID)) +
    theme_bw(base_size=14) +
    theme(axis.title.x = element_text(face="bold"), axis.text.x = element_text(face="bold")) +
    theme(axis.title.y = element_text(face="bold"), axis.text.y = element_text(face="bold")) +
    theme(plot.title = element_text(size=18, hjust=0.5)) +
    labs(list(x = "Visit", y = "Pseudo TB Data From Baseline"))#+geom_smooth()
  ## Now plot
  spider_plot <- spider_plot + geom_line(aes(color=treatment)) +
    geom_point(aes(color=treatment), show.legend=FALSE) +
    scale_colour_discrete(name="Treatment", labels=c("Control", "Drug")) +
    scale_shape_manual(name = "cstatus", values = c("0"=3, "1"=16)) +
    coord_cartesian(xlim=c(1, max(data_long$visit)))+
    facet_grid(~treatment)
  
  spider_plot
}








simulate_visits_in_days = function(surv_data_with_eos,
                                   max_visits                   #maximum number of visits (9 weeks apart) per individual
){
  return_df = data.frame(ID = surv_data_with_eos$ID, Trt = surv_data_with_eos$Trt, Date_enroll = surv_data_with_eos$Date_enroll,
                         Date_eos = surv_data_with_eos$Date_eos)
  for(i in 1:max_visits){
    visit_col_name = paste("day_visit", i, sep = "")
    #return_df[visit_col_name] = as.numeric(surv_data_with_eos$Date_eos - surv_data_with_eos$Date_enroll) - (max_visits-i+1)*9.6*7
    return_df[visit_col_name] = i*9.6*7 
  }
  return(as.data.frame(return_df))
}


# Stan Model


model_code3 <- "
data {  
  
  int<lower=0> n;                              // Number of individuals 

  
  // Data for survival submodel
  int<lower=0> Nobs;                           // Number of individuals who experienced death before censoring
  int<lower=0> Ncen;                           // Number of censored individuals
  real<lower=0> ObsDeathTime[Nobs];            // Observed death times
  real<lower=0> ObsCenTime[Ncen];              // Observed Censoring times
  int ID_surv_unordered[n];                    // vactor with individual ID in the unordered survival data
  int ID_order_obs_cen[n];                     // vector with individual ID (1:number_of_indv) where the individuals are ordered (1:Nobs, Nobs+1:Nobs+Ncen = n)
  vector[n] X1; 
  vector[n] A_i;                               // Baseline treatment variable
  
  // Data for longitudinal submodel
  int<lower=0> N;                              // Number of data points
  int<lower=0> K;                              // Number of fixed effects
  //vector[n] number_of_visits;                  // vector with number of visits for each individual
  int ID[N];                                   // vector (length N) with visit ID (1:number_of_visits) for each individual
  //matrix[N, K] x_pred;                          // matrix of fixed baseline predictors for longitudinal sub-model (currently only two baseline predictors)
  vector[N] t;                                 // fixed visits (different for each individual)
  vector[N] y;                                 // longitudinal outcome
}

parameters {
  // Paramters for the survival submodel
  real gamma_X1;
  real gamma_A;
  vector<lower=0>[Ncen] T_miss_raw;            // vector for missing death times
  
  // Paramters for the longitudinal submodel
  //vector[K] beta_longitudinal;                 // Fixed effects parameters for the longitudinal submodel (We have K many fixed effects)
  real beta_X1;
  real beta_A;                                 // Baseline Treatment parametersreal 
  real<lower = 0> sig_err;                                    // error s.d.
  
  
  //group level parameters
  vector[2] b_mu;                   //group level mean
  vector<lower = 0>[2] b_sd;                   //group level sd
  matrix[2, n] b_i;
  
  
}

transformed parameters{
  vector[n] Psi_i;                                // min(del_dag_T_dag, Ti)
  vector[N] psi_ij;                                // tij/Psi_i
  vector<lower = 0>[n] alpha;                   // parameter for ObsTime
  
  vector<lower=0>[Ncen] T_miss;
  
  
  
  
  
  // Paramters for the survival submodel
  for(i in 1:Nobs){
    Psi_i[ID_order_obs_cen[i]] = ObsDeathTime[i];
  }
  for(j in 1:Ncen){
    T_miss[j] = T_miss_raw[j] + ObsCenTime[j];  // make sure that T_miss[j] is atleast as greater as censoring times
    Psi_i[ID_order_obs_cen[Nobs + j]] = T_miss[j];
  }
  
  
  for(k in 1:n){
   alpha[k] =   exp(gamma_X1*X1[k] + gamma_A* A_i[k]);
   //alpha[k] =   exp(gamma_0 + gamma_A* A_i[k]);
  }
  
  for(l in 1:N){
    psi_ij[l] = t[l] / Psi_i[ID[l]];
  } 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}


model {
  vector[N] xi_fixed;     // x_pred * beta_longitudinal
  vector[N] mu_long;       // Longitudinal means
  
  
  // Priors for the survival submodel
  gamma_A ~ normal(0,1);
  //gamma_A ~ double_exponential(0, 1);
  gamma_X1 ~ normal(0,10);
  
  
  // Priors for the longitudinal submodel
  beta_A ~ normal(0,10);
  
  //beta_longitudinal ~ normal(0,10);
  beta_X1 ~ normal(0,10);
  
  
  //group level terms
  to_vector(b_sd) ~ gamma(1,1);
  to_vector(b_mu) ~ normal(0,100);
  
  for(k in 1:n){
    b_i[1,k] ~ normal(b_mu[1],b_sd[1]);
    b_i[2,k] ~ normal(b_mu[2],b_sd[2]);
    
  }
  
  
  
 
  
  //Likelihood for the survival submodel
  for(i in 1:Nobs){
     target+= exponential_lpdf(ObsDeathTime[i]|alpha[ID_order_obs_cen[i]]);    // observed death time
    
  }
  for(j in 1:Ncen){
    target+= exponential_lpdf(T_miss[j]|alpha[ID_order_obs_cen[Nobs + j]]);     // contribution of censored death time to the likelihood once //we know what their values are
    
  }
  
  
   //Likelihood for the longitudinal submodel
  //xi_fixed = x_pred * beta_longitudinal;
  for (i in 1:N) { 
     mu_long[i] =  X1[ID[i]] * beta_X1 + A_i[ID[i]] * beta_A +   (b_i[1,ID[i]])* psi_ij[i] +  (b_i[2,ID[i]])* pow(psi_ij[i],2);
     target += normal_lpdf(y[i]|mu_long[i], sig_err);
  }
  
  
  
}
"



# Function to calculate TB AUC for each subject across Q iterations 
calculate_TB_area_change_from_base = function(longitudinal_data,
                                              survival_data,
                                              fitted_Stan_model,
                                              num_iterations){
  setDT(longitudinal_data)
  data_longitudinal_wide_format = dcast(longitudinal_data, ID + A+  X1+  Y0   ~ visit, 
                                        value.var = c("Surv_time", "TB"))
  parameter_samples_joint = rstan::extract(fitted_Stan_model)
  
  n = nrow(survival_data)                            # number of patients
  A = matrix(NA,n,num_iterations) 
  B= matrix(NA,n,num_iterations) 
  C = matrix(NA,n,num_iterations)
  t = matrix(NA,n,num_iterations)
  TB_area = matrix(NA,n,num_iterations)
  #TB_area_grp_diff = matrix(NA,1,num_iterations)
  
  #Y0 = base_TB_df$Y0
  #obs_surv_data = survival_data[survival_data$Delta == 1,]
  #cen_surv_data = survival_data[survival_data$Delta == 0,]
  
  for(k in 1:num_iterations){
    t[,k] = ifelse(survival_data$Delta == 1, survival_data$ObsTime,parameter_samples_joint$T_miss[k,])
    t[,k] = pmin(t[,k],survival_data$Day_eos) 
    
    # constant term
    C[,k] =  (parameter_samples_joint$beta_X1[k] * data_longitudinal_wide_format$X1)+ 
      (parameter_samples_joint$beta_A[k] * data_longitudinal_wide_format$A) -
      data_longitudinal_wide_format$Y0
    
    # linear term
    B[,k] =  parameter_samples_joint$b_i[k,1,]/(2*parameter_samples_joint$Psi_i[k,])
    
    # quadratic term
    A[,k] =  parameter_samples_joint$b_i[k,2,]/(3* (parameter_samples_joint$Psi_i[k,])^2)
    
    
    
    
    TB_area[,k] =  (C[,k] * t[,k] + B[,k] * (t[,k]^2) + A[,k] * (t[,k]^3))/data_longitudinal_wide_format$Y0
    #TB_area[,k] =  B[,k] * (t^2) + A[,k] * (t^3)
    
    
    #TB_area_grp_diff[1,k] = sum(TB_area[,k][survival_data$Trt==1]) - sum(TB_area[,k][survival_data$Trt==0])
  }
  
  
  
  colnames(TB_area) = paste0("TB Area at Iter", 1:num_iterations)
  return_df = cbind(ID = survival_data$ID, Trt_indicator = data_longitudinal_wide_format$A,  TB_area)
  return(as.data.frame(return_df))
  #return(TB_area_grp_diff)
  
}




# Function to calculate survival AUC across Q iterations


calculate_survival_area = function(survival_data, fitted_Stan_model, survival_util, num_iterations){
  n = nrow(survival_data)
  t_star = matrix(NA,n,num_iterations)
  Surv_area = matrix(NA,n,num_iterations)
  #obs_surv_data = survival_data[survival_data$Delta == 1,]
  #cen_surv_data = survival_data[survival_data$Delta == 0,]
  parameter_samples_joint = rstan::extract(fitted_Stan_model)
  #Surv_area_grp_diff = matrix(NA,1,num_iterations)
  
  
  for(k in 1:num_iterations){
    t_star[,k] = ifelse(survival_data$Delta == 1, survival_data$ObsTime,parameter_samples_joint$T_miss[k,])
    no_surv_area_condition = (pmin(t_star[,k],survival_data$Day_eos) == survival_data$Day_eos)
    surv_area_vec = (survival_data$Day_eos-t_star[,k])* survival_util
    #print((surv_data_temp$Day_eos-t_star[,k]))
    Surv_area[,k] = ifelse(no_surv_area_condition, 0 ,  
                           surv_area_vec)
    
    #Surv_area_grp_diff[1,k] = sum(Surv_area[,k][survival_data$Trt==1]) - sum(Surv_area[,k][survival_data$Trt==0])
  }
  
  
  colnames(Surv_area) = paste0("Surv Area at Iter", 1:num_iterations)
  return_df = cbind(ID = survival_data$ID,Trt_indicator = survival_data$Trt,  Surv_area)
  return(as.data.frame(return_df))
  #return(Surv_area_grp_diff)
}




# Total Area under the curve across Q iterations

calculate_total_area = function(TB_Area_df, Surv_Area_df){
  TB_Area_df_slice = TB_Area_df[,3:ncol(TB_Area_df)]
  Surv_Area_df_slice = Surv_Area_df[,3:ncol(Surv_Area_df)]
  
  Total_Area_df_slice = TB_Area_df_slice + Surv_Area_df_slice
  num_iterations = ncol(Total_Area_df_slice)
  colnames(Total_Area_df_slice) = paste0("Total Area at Iter", 1:num_iterations)
  return_df = cbind(ID = Surv_Area_df$ID,Trt_indicator = Surv_Area_df$Trt,  Total_Area_df_slice)
  
  return(return_df)
}







# Wald Test

library(data.table)


get_BB_weights_df = function(num_subj, num_iter){
  return_df = t(as.data.frame(rdirichlet(num_iter, rep(1, num_subj)))) 
  return(return_df)
}


calculate_hyp_test_data_across_iters = function(AUC_df){  
  AUC_slice = AUC_df[,3:ncol(AUC_df)] 
  #AUC_slice_BB = AUC_slice * BB_weights_df #multiply corresponding columns of two dataframes of the same dimensions element-wise 
  n1 = length(unique(AUC_df$ID[AUC_df$Trt_indicator==1]))
  n0 = length(unique(AUC_df$ID[AUC_df$Trt_indicator==0]))
  
  
  calculate_mean_trt = function(auc_vector_iter){ 
    grp_mean_AUC = mean(auc_vector_iter[AUC_df$Trt_indicator ==1])
    return(grp_mean_AUC)
  }
  calculate_mean_ctrl = function(auc_vector_iter){ 
    grp_mean_AUC = mean(auc_vector_iter[AUC_df$Trt_indicator ==0])
    return(grp_mean_AUC)
  }
  
  calculate_BBweighted_mean_trt = function(auc_vector_iter){ 
    dir_weights = rdirichlet(1, rep(1, n1))[1, ]
    grp_mean_AUC = sum(auc_vector_iter[AUC_df$Trt_indicator ==1]*dir_weights)
    return(grp_mean_AUC)
  }
  calculate_BBweighted_mean_ctrl = function(auc_vector_iter){ 
    dir_weights = rdirichlet(1, rep(1, n0))[1, ]
    grp_mean_AUC = sum(auc_vector_iter[AUC_df$Trt_indicator ==0]*dir_weights)
    return(grp_mean_AUC)
  }
  
  # sample mean
  mean_AUC_trt <- apply(AUC_slice, 2, calculate_mean_trt)
  mean_AUC_ctrl <- apply(AUC_slice, 2, calculate_mean_ctrl)
  
  diff_mean_AUC <- as.numeric(mean_AUC_trt - mean_AUC_ctrl)
  
  # BB weighted mean
  mean_AUC_trt_BB <- apply(AUC_slice, 2, calculate_BBweighted_mean_trt)
  mean_AUC_ctrl_BB <- apply(AUC_slice, 2, calculate_BBweighted_mean_ctrl)
  
  diff_mean_AUC_BB <- as.numeric(mean_AUC_trt_BB - mean_AUC_ctrl_BB)
  
  return_df = data.frame(avg_trt_effect_iter = diff_mean_AUC, avg_trt_effect_BB_iter =diff_mean_AUC_BB)
  return(return_df)
}







perform_Wald_test = function(hyp_test_vec_across_iters_df ){   #df with ncol=2 (avg trt effect with and without BB), nrow = num_iters
  point_estimate = mean(hyp_test_vec_across_iters_df[,1])  # mean of avg trt effect without BB
  point_estimate_BB = mean(hyp_test_vec_across_iters_df[,2])  # mean of avg trt effect with BB
  standard_err= sd(hyp_test_vec_across_iters_df[,2])       # sd of avg trt effect with BB
  test_stat = point_estimate/standard_err
  
  p_val = pnorm(test_stat, mean = 0, sd = 1, lower.tail = TRUE)
  #p_val = pchisq(test_stat, df = 1)
  return_df = data.frame(point_estimate = point_estimate, point_estimate_BB = point_estimate_BB, rep_se = standard_err,test_stat= test_stat, rep_pval = p_val)
  return(return_df)
}




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Function Calls

# Generate baseline variables
baseline <- data.frame(
  ID = 1:n,
  A =  rbinom(n,1,0.5),  # Baseline treatment variable (binary)
  X1 = rnorm(n,6,1)  # Baseline covariate1 (continuous) (1,0.4) median needs to be >1 for positive treatment effect
  #X2 = rbinom(n,1,0.8)  # Baseline covariate2 (binary)
)





# Generate random effects
mean_vec = c(-10, 11)
#sd_vec = c(0.75,0.75)
sd_vec = c(1,1)

b_mat = generate_random_effects_indep_norm(2, 500, mean_vec,sd_vec)





enrollment_times_df = stb_tl_simu_enroll_arms(c(250,250),
                                              param_enroll  = list(type  = "by_rate",
                                                                   pt_per_mth = 50),
                                              #mth_min_fu  = 5,   #minimum follow-up months for all patients
                                              mth_fix_fu  = 18)     # fixed follow-up months for all patients



# Generate Survival Data
true_gamma_X1 = -1.2
data_survnoRE = generate_survdata_noRE(n = 500,baseline_data = baseline, gamma_X1 = true_gamma_X1,
                                       gamma_A = true_gamma_A, Enrollment_times_df = enrollment_times_df, rate_param_cen = 0.00128)

data_survival = get_endOfStudy_basedOnEvents(data_survnoRE, 120)
simulate_visits_days = simulate_visits_in_days(data_survival,8)  # 8 visits (9 weeks apart), resulting in 72/4 = 18 months





# Generate Longitudinal Data
true_sd_error = 0.25
true_beta_X1 = 1.5
true_beta_0 = 2.25


base_TB = generate_baseline_TB_data(num_subjects = 500, baseline_data = baseline)




data_long = generate_longitudinal_data(n = 500, max_visits = 10, baseline_data = baseline,
                                       visits_days_df = simulate_visits_days, surv_time = data_survival$ObsTime, true_death_time = data_survival$time_Death_joint,
                                       b_mat = b_mat, beta_A = true_beta_A, sd_error = true_sd_error, 
                                       beta_X1 = true_beta_X1)




# data_long_baseTB has base TB data
data_long_baseTB = combine_longitudinal_data(base_TB, data_long) 

data_long_change_from_base = generate_longitudinal_data_change_from_baseline(data_long, base_TB)

#Add new ID variable
data_long <- data_long %>%
  group_by(ID) %>%
  mutate(Input_ID = cur_group_id())




data_survival = data_survival[data_survival$ID %in% unique(data_long$ID),]
data_survival$Input_ID = 1:nrow(data_survival)

ID_order_obs_cen = c(as.numeric(data_survival[data_survival$Delta == 1, ]$Input_ID),
                     as.numeric(data_survival[data_survival$Delta == 0, ]$Input_ID))






# Prepare the data for Stan 

stan_data_joint1 <- list(
  n = nrow(data_survival),      # Number of individuals
  
  
  
  #survival data
  Nobs = nrow(data_survival[data_survival$Delta == 1,]),
  Ncen = nrow(data_survival[data_survival$Delta == 0,]),
  ObsDeathTime = data_survival$ObsTime[data_survival$Delta == 1],
  ObsCenTime = data_survival$ObsTime[data_survival$Delta == 0],
  ID_surv_unordered = data_survival$Input_ID,
  ID_order_obs_cen = ID_order_obs_cen,
  X1 = data_survival$X1,
  A_i = data_survival$Trt,
  
  #longitudinal data
  N = nrow(data_long),       # Number of observations
  K =2,                      # one intercept, one covariate
  #number_of_visits = as.vector(table(data_long$ID)),
  ID = data_long$Input_ID,
  #x_pred = data.frame(intercept = rep(1,nrow(data_long)), X1 = data_long$X1),
  t = data_long$visit_days,
  y = data_long$TB
)


# Run Stan Model
model3 <- stan_model(model_code = model_code3)
fit3 <- sampling(model3, data = stan_data_joint1, cores = 4,  chains = 4, warmup = 1000, iter = 7250)
parameter_samples_joint <- rstan::extract(fit3)





TB_Area = calculate_TB_area_change_from_base(longitudinal_data = data_long_change_from_base,
                                             survival_data = data_survival, 
                                             fitted_Stan_model = fit3,
                                             num_iterations = 25000)
TB_AUC_hyp_test_data = calculate_hyp_test_data_across_iters(TB_Area)
TB_rep_Wald_test = perform_Wald_test(TB_AUC_hyp_test_data)


Surv_Area = calculate_survival_area(data_survival,fit3, 0.5,25000)
Surv_AUC_hyp_test_data = calculate_hyp_test_data_across_iters(Surv_Area)
Surv_rep_Wald_test = perform_Wald_test(Surv_AUC_hyp_test_data)



Total_Area = calculate_total_area(TB_Area, Surv_Area)
Total_AUC_hyp_test_data = calculate_hyp_test_data_across_iters(Total_Area)
Total_AUC_rep_Wald_test = perform_Wald_test(Total_AUC_hyp_test_data)




# Compute Average Posterior Means
# gamma_X1_post_mean = mean(parameter_samples_joint$gamma_X1)
# 
# beta_A_post_mean = mean(parameter_samples_joint$beta_A)
# 
# gamma_A_post_mean = mean(parameter_samples_joint$gamma_A)
# 
# sig_err_post_mean = mean(parameter_samples_joint$sig_err)
# 
# beta_X1_post_mean = mean(parameter_samples_joint$beta_longitudinal[,1])
# 
# beta_X2_post_mean = mean(parameter_samples_joint$beta_longitudinal[,2])
# 
# b_mu1_post_mean = mean(parameter_samples_joint$b_mu[,1])
# 
# b_mu2_post_mean = mean(parameter_samples_joint$b_mu[,2])
# 
# 
# b_sd1_post_mean = mean(parameter_samples_joint$b_sd[,1])
# 
# b_sd2_post_mean = mean(parameter_samples_joint$b_sd[,2])
# 
















# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

allinfo = data.frame( runID = run_ID,
                      gamma_X1_post_mean= mean(parameter_samples_joint$gamma_X1),
                      gamma_A_post_mean= mean(parameter_samples_joint$gamma_A),
                      beta_A_post_mean = mean(parameter_samples_joint$beta_A), 
                      sig_err_post_mean = mean(parameter_samples_joint$sig_err),
                      beta_X1_post_mean = mean(parameter_samples_joint$beta_X1), 
                      b_mu1_post_mean = mean(parameter_samples_joint$b_mu[,1]),
                      b_mu2_post_mean = mean(parameter_samples_joint$b_mu[,2]),
                      b_sd1_post_mean = mean(parameter_samples_joint$b_sd[,1]),
                      b_sd2_post_mean = mean(parameter_samples_joint$b_sd[,2]),
                      TB_point_est = TB_rep_Wald_test$point_estimate,
                      TB_point_est_BB = TB_rep_Wald_test$point_estimate_BB,
                      TB_se = TB_rep_Wald_test$rep_se,
                      TB_test_stat = TB_rep_Wald_test$test_stat,
                      TB_pval = TB_rep_Wald_test$rep_pval,
                      Surv_point_est = Surv_rep_Wald_test$point_estimate,
                      Surv_point_est_BB = Surv_rep_Wald_test$point_estimate_BB,
                      Surv_se = Surv_rep_Wald_test$rep_se,
                      Surv_test_stat = Surv_rep_Wald_test$test_stat,
                      Surv_pval = Surv_rep_Wald_test$rep_pval,
                      Tot_point_est = Total_AUC_rep_Wald_test$point_estimate,
                      Tot_point_est_BB = Total_AUC_rep_Wald_test$point_estimate_BB,
                      Tot_se = Total_AUC_rep_Wald_test$rep_se,
                      Tot_test_stat = Total_AUC_rep_Wald_test$test_stat,
                      Tot_pval = Total_AUC_rep_Wald_test$rep_pval)
write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
#save(allinfo, file = txt.title)
rm(list="parameter_samples_joint")
rm(list="TB_rep_Wald_test")
rm(list="Surv_rep_Wald_test")
rm(list="Total_AUC_rep_Wald_test")
gc()