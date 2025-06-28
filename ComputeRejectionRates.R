# Results after running the simulation for 1000 replicated datasets

## Load Sumilation Results (Bayesian Bootstrap) 


library(ggplot2)
library(gridExtra)

# df_col_names = c("run_ID",
#                  "gamma_X1_post_mean", "gamma_A_post_mean",
#                  "beta_A_post_mean", "sig_err_post_mean", "beta_X1_post_mean", 
#                  "b_mu1_post_mean", "b_mu2_post_mean", "b_sd1_post_mean", "b_sd2_post_mean",
#                  "TB_point_est", "TB_se", "TB_test_stat", "TB_pval",
#                  "Surv_point_est", "Surv_se", "Surv_test_stat", "Surv_pval",
#                  "Tot_point_est", "Tot_se","Tot_test_stat", "Tot_pval")


sim1_result_df = data.frame(data.table::fread("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Results_TBproj_Sim1_result.txt", sep = "\t"))

head(sim1_result_df)


sim2_result_df = data.frame(data.table::fread("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Results_TBproj_Sim2_result.txt", sep = "\t"))

head(sim2_result_df)


sim3_result_df = data.frame(data.table::fread("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Results_TBproj_Sim3_result.txt", sep = "\t"))


head(sim3_result_df)



sim4_result_df = data.frame(data.table::fread("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Results_TBproj_Sim4_result.txt", sep = "\t"))

head(sim4_result_df)






















## Parameter Estimation




parameter_estimation = function(sim_estimation_results_df,
                                true_sd_error = 0.25,
                                true_beta_A,
                                true_beta_X1 = 1.5, 
                                true_gamma_X1 = -1.2,
                                true_gamma_A,
                                true_b_mu1 = -10,
                                true_b_mu2 = 11,
                                true_b_sd1 = 1,
                                true_b_sd2 = 1){
  
  
  true_param = c(true_sd_error,true_beta_A, true_beta_X1, 
                 true_gamma_X1, true_gamma_A, true_b_mu1, true_b_mu2, true_b_sd1, true_b_sd2)
  est = c(mean(sim_estimation_results_df$sig_err_post_mean), mean(sim_estimation_results_df$beta_A_post_mean),
          mean(sim_estimation_results_df$beta_X1_post_mean), 
          mean(sim_estimation_results_df$gamma_X1_post_mean), mean(sim_estimation_results_df$gamma_A_post_mean),
          mean(sim_estimation_results_df$b_mu1_post_mean), mean(sim_estimation_results_df$b_mu2_post_mean),
          mean(sim_estimation_results_df$b_sd1_post_mean), mean(sim_estimation_results_df$b_sd2_post_mean))
  bias = est - true_param
  mean_sq_err = c(mean((true_sd_error - sim_estimation_results_df$sig_err_post_mean)^2) ,
                  mean((true_beta_A - sim_estimation_results_df$beta_A_post_mean)^2) ,
                  mean((true_beta_X1 - sim_estimation_results_df$beta_X1_post_mean)^2) ,
                  mean((true_gamma_X1 - sim_estimation_results_df$gamma_X1_post_mean)^2) ,
                  mean((true_gamma_A - sim_estimation_results_df$gamma_A_post_mean)^2),
                  mean(true_b_mu1 - sim_estimation_results_df$b_mu1_post_mean)^2,
                  mean(true_b_mu2 - sim_estimation_results_df$b_mu2_post_mean)^2,
                  mean(true_b_sd1 - sim_estimation_results_df$b_sd1_post_mean)^2,
                  mean(true_b_sd2 - sim_estimation_results_df$b_sd2_post_mean)^2)
  post_param_df <- data.frame( TrueValues= true_param,
                               Estimates= est, 
                               Bias= bias, 
                               MSE= mean_sq_err) 
  
  
  return_df = as.data.frame(lapply(post_param_df, function(x) {
    if(is.numeric(x)) {
      format(round(x, 3), nsmall = 3)
    } else {
      x
    }
  }))
  rownames(return_df) = c("sd_error","beta_A", "beta_X1","gamma_X1", "gamma_A" ,"b_mu1","b_mu2", "b_sd1", "b_sd2")
  return(return_df)
  
}

  





### Sim 1 (gamma_A = -0.75, beta_A = -2.25) parameter estimation

  
sim1_paramter_estimation = parameter_estimation(sim1_result_df, true_beta_A = -2.25, true_gamma_A = -0.75)
sim1_paramter_estimation

  

### Sim 2  (gamma_A = 0, beta_A = -2.25) parameter estimation


 
sim2_paramter_estimation = parameter_estimation(sim2_result_df, true_beta_A = -2.25, true_gamma_A = 0)
sim2_paramter_estimation

  

### Sim 3  (gamma_A = -0.75, beta_A = 0) parameter estimation




sim3_paramter_estimation = parameter_estimation(sim3_result_df, true_beta_A = 0, true_gamma_A = -0.75)
sim3_paramter_estimation

  

### Sim 4 (gamma_A = 0, beta_A = 0) parameter estimation




sim4_paramter_estimation = parameter_estimation(sim4_result_df, true_beta_A = 0, true_gamma_A = 0)
sim4_paramter_estimation












## Rejection rates



calculate_rejectionRate = function(sim_estimation_results_df){
  
  TB_rejectionRate =  mean(sim_estimation_results_df$TB_pval < 0.025)
  Surv_rejectionRate =  mean(sim_estimation_results_df$Surv_pval < 0.025)
  Tot_rejectionRate =  mean(sim_estimation_results_df$Tot_pval < 0.025)
  
  
  rejectionRate = c(TB_rejectionRate, Surv_rejectionRate,  Tot_rejectionRate)
  error_rate_df <- data.frame( rejection_rate= rejectionRate)     
  
  rownames(error_rate_df) = c("TB_AUC","Surv_AUC", "Total_AUC")  
  return(error_rate_df)
  
}







### Sim 1 (gamma_A = -0.75, beta_A = -2.5) error rates

  
sim1_err_rate = calculate_rejectionRate(sim1_result_df)
sim1_err_rate





Sim1_test_stat_plot1 = ggplot(sim1_result_df, aes(x  = TB_test_stat)) +
  geom_density() +
  labs(title = "TB AUC test stat ") +
  theme_minimal()


Sim1_test_stat_plot2 = ggplot(sim1_result_df, aes(x = Surv_test_stat)) +
  geom_density() +
  labs(title = "Surv AUC test stat ") +
  theme_minimal()



Sim1_test_stat_plot3 = ggplot(sim1_result_df, aes(x =  Tot_test_stat)) +
  geom_density() +
  labs(title = "Total AUC test stat ") +
  theme_minimal()







grid.arrange(Sim1_test_stat_plot1, Sim1_test_stat_plot2, Sim1_test_stat_plot3,  nrow = 1)








# Generate theoretical quantiles from Normal(0,1) distribution
theoretical_quantiles1_ts = qnorm(ppoints(nrow(sim1_result_df)))

y11_ts = sort(sim1_result_df$TB_test_stat)
y12_ts = sort(sim1_result_df$Surv_test_stat)
y13_ts = sort(sim1_result_df$Tot_test_stat)

Sim1_test_stat_QQplot1 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1_ts, y = y11_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()


Sim1_test_stat_QQplot2 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1_ts, y = y12_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()

Sim1_test_stat_QQplot3 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1_ts, y = y13_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()





grid.arrange(Sim1_test_stat_QQplot1, Sim1_test_stat_QQplot2, Sim1_test_stat_QQplot3,  nrow = 1)










Sim1_pval_plot1 = ggplot(sim1_result_df, aes(x = run_ID, y = TB_pval)) +
  geom_point() +
  labs(title = "TB p-value",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


Sim1_pval_plot2 = ggplot(sim1_result_df, aes(x = run_ID, y = Surv_pval)) +
  geom_point() +
  labs(title = "Survival p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()



Sim1_pval_plot3 = ggplot(sim1_result_df, aes(x =  run_ID, y = Tot_pval)) +
  geom_point() +
  labs(title = "Total p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


grid.arrange(Sim1_pval_plot1, Sim1_pval_plot2, Sim1_pval_plot3,  nrow = 1)



# Generate theoretical quantiles from Uniform(0,1) distribution
theoretical_quantiles1 = qunif(ppoints(nrow(sim1_result_df)))




y11 = sort(sim1_result_df$TB_pval)
y12 = sort(sim1_result_df$Surv_pval)
y13 = sort(sim1_result_df$Tot_pval)

Sim1_pval_QQplot1 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1, y = y11)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()


Sim1_pval_QQplot2 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1, y = y12)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()

Sim1_pval_QQplot3 = ggplot(sim1_result_df, aes(x = theoretical_quantiles1, y = y13)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()





grid.arrange(Sim1_pval_QQplot1, Sim1_pval_QQplot2, Sim1_pval_QQplot3,  nrow = 1)

  

### Sim 2 (gamma_A = 0, beta_A = -2.5) error rates


 
sim2_err_rate = calculate_rejectionRate(sim2_result_df)
sim2_err_rate




Sim2_test_stat_plot1 = ggplot(sim2_result_df, aes(x  = TB_test_stat)) +
  geom_density() +
  labs(title = "TB AUC test stat ") +
  theme_minimal()


Sim2_test_stat_plot2 = ggplot(sim2_result_df, aes(x = Surv_test_stat)) +
  geom_density() +
  labs(title = "Surv AUC test stat ") +
  theme_minimal()



Sim2_test_stat_plot3 = ggplot(sim2_result_df, aes(x =  Tot_test_stat)) +
  geom_density() +
  labs(title = "Total AUC test stat ") +
  theme_minimal()







grid.arrange(Sim2_test_stat_plot1, Sim2_test_stat_plot2, Sim2_test_stat_plot3,  nrow = 1)








# Generate theoretical quantiles from Normal(0,1) distribution
theoretical_quantiles2_ts = qnorm(ppoints(nrow(sim2_result_df)))

y21_ts = sort(sim2_result_df$TB_test_stat)
y22_ts = sort(sim2_result_df$Surv_test_stat)
y23_ts = sort(sim2_result_df$Tot_test_stat)

Sim2_test_stat_QQplot1 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2_ts, y = y21_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()


Sim2_test_stat_QQplot2 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2_ts, y = y22_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()

Sim2_test_stat_QQplot3 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2_ts, y = y23_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()





grid.arrange(Sim2_test_stat_QQplot1, Sim2_test_stat_QQplot2, Sim2_test_stat_QQplot3,  nrow = 1)










Sim2_pval_plot1 = ggplot(sim2_result_df, aes(x = run_ID, y = TB_pval)) +
  geom_point() +
  labs(title = "TB p-value",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


Sim2_pval_plot2 = ggplot(sim2_result_df, aes(x = run_ID, y = Surv_pval)) +
  geom_point() +
  labs(title = "Survival p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()



Sim2_pval_plot3 = ggplot(sim2_result_df, aes(x =  run_ID, y = Tot_pval)) +
  geom_point() +
  labs(title = "Total p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


grid.arrange(Sim2_pval_plot1, Sim2_pval_plot2, Sim2_pval_plot3,  nrow = 1)





# Generate theoretical quantiles from Uniform(0,1) distribution
theoretical_quantiles2 = qunif(ppoints(nrow(sim2_result_df)))




y21 = sort(sim2_result_df$TB_pval)
y22 = sort(sim2_result_df$Surv_pval)
y23 = sort(sim2_result_df$Tot_pval)

Sim2_pval_QQplot1 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2, y = y21)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()


Sim2_pval_QQplot2 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2, y = y22)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()

Sim2_pval_QQplot3 = ggplot(sim2_result_df, aes(x = theoretical_quantiles2, y = y23)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()





grid.arrange(Sim2_pval_QQplot1, Sim2_pval_QQplot2, Sim2_pval_QQplot3,  nrow = 1)


  

### Sim 3 (gamma_A = -0.75, beta_A = 0) error rates




sim3_err_rate = calculate_rejectionRate(sim3_result_df)
sim3_err_rate






Sim3_test_stat_plot1 = ggplot(sim3_result_df, aes(x  = TB_test_stat)) +
  geom_density() +
  labs(title = "TB AUC test stat ") +
  theme_minimal()


Sim3_test_stat_plot2 = ggplot(sim3_result_df, aes(x = Surv_test_stat)) +
  geom_density() +
  labs(title = "Surv AUC test stat ") +
  theme_minimal()



Sim3_test_stat_plot3 = ggplot(sim3_result_df, aes(x =  Tot_test_stat)) +
  geom_density() +
  labs(title = "Total AUC test stat ") +
  theme_minimal()







grid.arrange(Sim3_test_stat_plot1, Sim3_test_stat_plot2, Sim3_test_stat_plot3,  nrow = 1)






# Generate theoretical quantiles from Normal(0,1) distribution
theoretical_quantiles3_ts = qnorm(ppoints(nrow(sim3_result_df)))

y31_ts = sort(sim3_result_df$TB_test_stat)
y32_ts = sort(sim3_result_df$Surv_test_stat)
y33_ts = sort(sim3_result_df$Tot_test_stat)

Sim3_test_stat_QQplot1 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3_ts, y = y31_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()


Sim3_test_stat_QQplot2 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3_ts, y = y32_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()

Sim3_test_stat_QQplot3 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3_ts, y = y33_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()





grid.arrange(Sim3_test_stat_QQplot1, Sim3_test_stat_QQplot2, Sim3_test_stat_QQplot3,  nrow = 1)






Sim3_pval_plot1 = ggplot(sim3_result_df, aes(x = run_ID, y = TB_pval)) +
  geom_point() +
  labs(title = "TB p-value",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


Sim3_pval_plot2 = ggplot(sim3_result_df, aes(x = run_ID, y = Surv_pval)) +
  geom_point() +
  labs(title = "Survival p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()



Sim3_pval_plot3 = ggplot(sim3_result_df, aes(x =  run_ID, y = Tot_pval)) +
  geom_point() +
  labs(title = "Total p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


grid.arrange(Sim3_pval_plot1, Sim3_pval_plot2, Sim3_pval_plot3,  nrow = 1)




# Generate theoretical quantiles from Uniform(0,1) distribution
theoretical_quantiles3 = qunif(ppoints(nrow(sim3_result_df)))




y31 = sort(sim3_result_df$TB_pval)
y32 = sort(sim3_result_df$Surv_pval)
y33 = sort(sim3_result_df$Tot_pval)

Sim3_pval_QQplot1 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3, y = y31)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()


Sim3_pval_QQplot2 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3, y = y32)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()

Sim3_pval_QQplot3 = ggplot(sim3_result_df, aes(x = theoretical_quantiles3, y = y33)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()





grid.arrange(Sim3_pval_QQplot1, Sim3_pval_QQplot2, Sim3_pval_QQplot3,  nrow = 1)

  

### Sim 4 (gamma_A = 0, beta_A = 0) error rates



sim4_err_rate = calculate_rejectionRate(sim4_result_df)
sim4_err_rate







Sim4_test_stat_plot1 = ggplot(sim4_result_df, aes(x  = TB_test_stat)) +
  geom_density() +
  labs(title = "TB AUC test stat ") +
  theme_minimal()


Sim4_test_stat_plot2 = ggplot(sim4_result_df, aes(x = Surv_test_stat)) +
  geom_density() +
  labs(title = "Surv AUC test stat ") +
  theme_minimal()



Sim4_test_stat_plot3 = ggplot(sim4_result_df, aes(x =  Tot_test_stat)) +
  geom_density() +
  labs(title = "Total AUC test stat ") +
  theme_minimal()







grid.arrange(Sim4_test_stat_plot1, Sim4_test_stat_plot2, Sim4_test_stat_plot3,  nrow = 1)










# Generate theoretical quantiles from Normal(0,1) distribution
theoretical_quantiles4_ts = qnorm(ppoints(nrow(sim4_result_df)))

y41_ts = sort(sim4_result_df$TB_test_stat)
y42_ts = sort(sim4_result_df$Surv_test_stat)
y43_ts = sort(sim4_result_df$Tot_test_stat)

Sim4_test_stat_QQplot1 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4_ts, y = y41_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()


Sim4_test_stat_QQplot2 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4_ts, y = y42_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()

Sim4_test_stat_QQplot3 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4_ts, y = y43_ts)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total test_stat",
       x = "Theoretical Quantiles",
       y = "test_stat") +
  theme_minimal()





grid.arrange(Sim4_test_stat_QQplot1, Sim4_test_stat_QQplot2, Sim4_test_stat_QQplot3,  nrow = 1)










Sim4_pval_plot1 = ggplot(sim4_result_df, aes(x = run_ID, y = TB_pval)) +
  geom_point() +
  labs(title = "TB p-value",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


Sim4_pval_plot2 = ggplot(sim4_result_df, aes(x = run_ID, y = Surv_pval)) +
  geom_point() +
  labs(title = "Survival p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()



Sim4_pval_plot3 = ggplot(sim4_result_df, aes(x =  run_ID, y = Tot_pval)) +
  geom_point() +
  labs(title = "Total p-value ",
       x = "Replication",
       y = "p-value") +
  theme_minimal()


grid.arrange(Sim4_pval_plot1, Sim4_pval_plot2, Sim4_pval_plot3,  nrow = 1)




# Generate theoretical quantiles from Uniform(0,1) distribution
theoretical_quantiles4 = qunif(ppoints(nrow(sim4_result_df)))




y41 = sort(sim4_result_df$TB_pval)
y42 = sort(sim4_result_df$Surv_pval)
y43 = sort(sim4_result_df$Tot_pval)

Sim4_pval_QQplot1 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4, y = y41)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: TB p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()


Sim4_pval_QQplot2 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4, y = y42)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Survival p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()

Sim4_pval_QQplot3 = ggplot(sim4_result_df, aes(x = theoretical_quantiles4, y = y43)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot: Total p-value",
       x = "Theoretical Quantiles",
       y = "p-value") +
  theme_minimal()





grid.arrange(Sim4_pval_QQplot1, Sim4_pval_QQplot2, Sim4_pval_QQplot3,  nrow = 1)



























