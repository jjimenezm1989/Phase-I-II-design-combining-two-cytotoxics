

trial_simulation = function(true_rho00, true_rho01, true_rho10, true_alpha3, true_beta0, true_beta1, true_beta2,
                            true_beta3, true_beta4, true_beta5, path_JAGS){

    
  jags_txt_stage1 = paste0(path_JAGS,"/stage_one.bug.txt")
  jags_txt_stage2 = paste0(path_JAGS,"/marginal_models.bug.txt")
    
  
  # Declaration of the output statistics
  interim_prob_dlt = numeric()
  interim_prob_eff = numeric()
  
  # Beginning of study simulation
  
  #First 2 patients
  
  X = c(0,0)
  Y = c(0,0)
  n = 2
  sample_size_stopping_safety = 2
  stopping_safety_ind = 0
  
  
  outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                          beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[1], y = Y[1])
  
  outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                          beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[2], y = Y[2])
  
  Z = c(outcomes_1[1],outcomes_2[1])
  E = c(outcomes_1[2],outcomes_2[2])
    
  jags <- jags.model(jags_txt_stage1,data = list('theta_Z'= theta_Z,'N' = n,'X' = X,'Y' = Y,'Z' = Z), n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(jags,c('mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3"),mm,progress.bar = "none")
  ss=as.data.frame(s[[1]])
  
  #Calculating posterior probability of toxicity at x = 0, y = 0.
  #We create, for each MCMC sample, the correspondent P(DLT). 
  #We then calculate the mean posterior P(DLT) at x = 0, y = 0
  ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
                      temp2 = ifelse(temp1 > theta_Z + 0.1, 1, 0))
  
  interim_prob_dlt[1] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
  
  #We assess whether we will enroll a new cohort or not based on safety stopping rule
  stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[1] <= delta_Z1), 0, 1)
  sample_size_stopping_safety = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[1] <= delta_Z1), sample_size_stopping_safety + 2, sample_size_stopping_safety)
  
  
  trcmtdx1 = ss$mtdx1[ss$mtdx1 > lbx]
  xx1 = max(0,quantile(trcmtdx1,alpha[1]))
  xx1 = min(xx1,1)
  if ((xx1 - X[1]) > delta1)
    xx1 = X[1]+delta1
  
  trcmtdy2 = ss$mtdy2[ss$mtdy2 > lby]
  yy2 = max(0,quantile(trcmtdy2,alpha[1]))
  yy2 = min(yy2,1)
  if ((yy2 - Y[2]) > delta1)
    yy2 = Y[2]+delta1
  
  # Update data
  
  X = c(X,0,xx1)
  Y = c(Y,yy2,0)
  n = length(X)
  
  #Patients 3 and 4
  
  outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                          beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[3], y = Y[3])
  
  outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                          beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                          beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[4], y = Y[4])
  
  Z = c(Z,outcomes_1[1],outcomes_2[1])
  E = c(E,outcomes_1[2],outcomes_2[2])
  
  
  for (c1 in 3:(NN_1/2)){
    
    jags <- jags.model(jags_txt_stage1, data = list('theta_Z'= theta_Z,'N' = n,'X' = X,'Y' = Y,'Z' = Z), n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(jags,c('mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3"),mm,progress.bar = "none")
    ss=as.data.frame(s[[1]])
    
    #Calculating posterior probability of toxicity at x = 0, y = 0.
    #We create, for each MCMC sample, the correspondent P(DLT). 
    #We then calculate the mean posterior P(DLT) at x = 0, y = 0
    ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
             temp2 = ifelse(temp1 > theta_Z + 0.1, 1, 0))
    
    interim_prob_dlt[c1-1] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
    
    #We assess whether we will enrol a new cohort or not based on safety stopping rule
    stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[c1-1] <= delta_Z1), 0, 1)
    sample_size_stopping_safety = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[c1-1] <= delta_Z1), sample_size_stopping_safety + 2, sample_size_stopping_safety)
    
    
    trcmtdx1 = ss$mtdx1[ss$mtdx1 > lbx]
    xx1 = max(0,quantile(trcmtdx1,alpha[c1-1]))
    xx1 = min(xx1,1)
    if ((xx1 - X[2*c1-3]) > delta1)
      xx1 = X[2*c1-3]+delta1
    
    trcmtdx2 = ss$mtdx2[ss$mtdx2 > lbx]
    xx2 = max(0,quantile(trcmtdx2,alpha[c1-1]))
    xx2 = min(xx2,1)
    if ((xx2 - X[2*c1-2]) > delta1)
      xx2 = X[2*c1-2]+delta1
    
    trcmtdy1 = ss$mtdy1[ss$mtdy1 > lby]
    yy1 = max(0,quantile(trcmtdy1,alpha[c1-1]))
    yy1 = min(yy1,1)
    if ((yy1 - Y[2*c1-3]) > delta1)
      yy1 = Y[2*c1-3]+delta1
    
    trcmtdy2 = ss$mtdy2[ss$mtdy2 > lby]
    yy2 = max(0,quantile(trcmtdy2,alpha[c1-1]))
    yy2 = min(yy2,1)
    if ((yy2 - Y[2*c1-2]) > delta1)
      yy2 = Y[2*c1-2]+delta1
    
    
    if (X[2*c1-3] == X[2*c1-5]){
      X = c(X,xx1,X[2*c1-2])
      Y = c(Y,Y[2*c1-3],yy2)
    }else{
      X = c(X,X[2*c1-3],xx2)
      Y = c(Y,yy1,Y[2*c1-2])
    }
    
    
    n = length(X)
    
    outcomes_1 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                            beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                            beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[2*c1-1], y = Y[2*c1-1])
    
    outcomes_2 = generate_bivariate_outcome(rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                                            beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, 
                                            beta4 = true_beta4, beta5 = true_beta5, corr = 0, x = X[2*c1-2], y = Y[2*c1-2])
    
    Z = c(Z,outcomes_1[1],outcomes_2[1])
    E = c(E,outcomes_1[2],outcomes_2[2])
    
  }
  
  #STAGE II
  
  stopping_futility_ind = 0
  sample_size_stopping_futility = NN_1 + m_cohort
  
  n_dlts = sum(Z)
  n_doses = length(Z)
    
  jags <- jags.model(jags_txt_stage2,data = list('theta_Z'= theta_Z,'N' = n,'X' = X,'Y' = Y,'E' = E,'Z' = Z, 'n_dlts' = n_dlts, 'n_doses' = n_doses), n.chains=chains,n.adapt=burn,quiet = TRUE)
  s=coda.samples(jags,c('safety_stopping','dlt_rate','mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3", "beta0", "beta1", "beta2", "beta3","beta4","beta5"),mm,progress.bar = "none")
  ss=as.data.frame(s[[1]])
  
  #Calculating posterior probability of toxicity at x = 0, y = 0.
  #We create, for each MCMC sample, the correspondent P(DLT). 
  #We then calculate the mean posterior P(DLT) at x = 0, y = 0
  ss2 = ss %>% mutate(temp1 = pdlt(rho00,rho01,rho10,alpha3, x = 0, y = 0),
                      temp2 = ifelse(temp1 > theta_Z + 0.1, 1, 0))
  
  interim_prob_dlt[NN_1/2] = dplyr::pull(ss2 %>% summarise(m = mean(temp2)))
  
  #We assess whether we will enroll a new cohort or not based on safety stopping rule
  stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[NN_1/2] <= delta_Z1), 0, 1)
  sample_size_stopping_safety = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[NN_1/2] <= delta_Z1), sample_size_stopping_safety + m_cohort, sample_size_stopping_safety)
  
  
  #Safety and efficacy model parameters after Stage 1.
  post_beta0 = median(ss$beta0)
  post_beta1 = median(ss$beta1)
  post_beta2 = median(ss$beta2)
  post_beta3 = median(ss$beta3)
  post_beta4 = median(ss$beta4)
  post_beta5 = median(ss$beta5)
  
  post_rho00 = median(ss$rho00)
  post_rho10 = median(ss$rho10)
  post_rho01 = median(ss$rho01)
  post_alpha3 = median(ss$alpha3)
  
  #In this algorithm, we generate samples from the standardized density F(beta0+beta1*x+beta2*y+beta3*x*y)
  
  for(c2 in 1:(NN_2/m_cohort)){
    
    target_function = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta_Z,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      p = peff(post_beta0,post_beta1,post_beta2,post_beta3,post_beta4,post_beta5,x,y_mtd)
      return(p)
    }
    
    c_constant = dplyr::pull(data.frame(x_grid = x_grid) %>% 
                               mutate(y = target_function(x_grid)) %>%
                               summarise(tempmax = max(y)))
    
    #Sample from the standardized density based on the rejection-sampling principle.
    accepted_sample_x = numeric()
    accepted_sample_y = numeric()
    accepted_sample = 0
    
    #We need m_cohort samples that are admissible and are accepted by the rejection sample algorithm.
    
    while(length(accepted_sample_x) < m_cohort){
      
      sample_x = runif(n = 1, min = 0, max = 1)
      U = runif(n = 1, min = 0, max = 1)
      
      acceptance_ind = ifelse(U <= (target_function(sample_x)/(c_constant*dunif(sample_x,0,1))), 1, 0)
      
      if(acceptance_ind == 1){
        accepted_sample = accepted_sample + 1
        accepted_sample_x[accepted_sample] = sample_x  
        accepted_sample_y[accepted_sample] = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta_Z,sample_x)
      }
    }
    
    
    accepted_sample_y = ifelse(accepted_sample_y < 0 , 0, ifelse(accepted_sample_y > 1, 1, accepted_sample_y))
    
    xx1 = accepted_sample_x
    yy1 = accepted_sample_y

    outcomes = mapply(generate_bivariate_outcome, rho00 = true_rho00, rho01 = true_rho01, rho10 = true_rho10, alpha3 = true_alpha3,
                      beta0 = true_beta0, beta1 = true_beta1, beta2 = true_beta2, beta3 = true_beta3, beta4 = true_beta4, beta5 = true_beta5,
                      corr = 0, x = xx1, y = yy1)
    
    Z = c(Z,outcomes[1,])
    E = c(E,outcomes[2,])
    X = c(X,xx1)
    Y = c(Y,yy1)
    
    n = length(X)
    
    n_dlts = sum(Z)
    n_doses = length(Z)

      
    jags <- jags.model(jags_txt_stage2,data = list('theta_Z'= theta_Z,'N' = n,'X' = X,'Y' = Y,'E' = E,'Z' = Z, 'n_dlts' = n_dlts, 'n_doses' = n_doses), n.chains=chains,n.adapt=burn,quiet = TRUE)
    s=coda.samples(jags,c('safety_stopping','dlt_rate','mtdx1','mtdx2','mtdy1','mtdy2',"rho00","rho10","rho01","alpha3", "beta0", "beta1", "beta2", "beta3","beta4","beta5"),mm,progress.bar = "none")
    ss=as.data.frame(s[[1]])
    
    #Calculating posterior probability of toxicity in Stage II.
    interim_prob_dlt[c2 + (NN_1/2)] = dplyr::pull(ss %>% summarise(stopping_rule_stage2 = mean(safety_stopping)))
    
    #We assess whether we will enroll a new cohort or not based on safety stopping rule
    stopping_safety_ind = ifelse(stopping_safety_ind == 0 & (interim_prob_dlt[c2 + (NN_1/2)] <= delta_Z2), 0, 1)
    sample_size_stopping_safety = ifelse(stopping_safety_ind == 0 & interim_prob_dlt[c2 + (NN_1/2)] <= delta_Z2 & (sample_size_stopping_safety < (NN_1 + NN_2)), sample_size_stopping_safety + m_cohort, sample_size_stopping_safety)
    
    post_beta0  = median(ss$beta0)
    post_beta1  = median(ss$beta1)
    post_beta2  = median(ss$beta2)
    post_beta3  = median(ss$beta3)
    post_beta4  = median(ss$beta4)
    post_beta5  = median(ss$beta5)
    
    post_rho00  = median(ss$rho00)
    post_rho10  = median(ss$rho10)
    post_rho01  = median(ss$rho01)
    post_alpha3 = median(ss$alpha3)
    
    post_prob_eff = function(x){
      y_mtd = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta_Z,x)
      y_mtd = ifelse(y_mtd < 0, 0, ifelse(y_mtd > 1, 1, y_mtd))
      
      df_posterior = data.frame(beta0 = ss$beta0, beta1 = ss$beta1, beta2 = ss$beta2, 
                                beta3 = ss$beta3, beta4 = ss$beta4, beta5 = ss$beta5, x = x, y = y_mtd) %>%
        mutate(ind = 1:n()) %>%
        group_by(ind) %>%
        mutate(peff = peff(beta0,beta1,beta2,beta3,beta4,beta5,x,y),
               condition = ifelse(peff > theta_E,1,0)) %>%
        ungroup() %>%
        summarise(out = mean(condition))
      return(df_posterior$out)
    }
    
    post_prob_efficacy_estimated_MTD = as.numeric(lapply(x_grid, post_prob_eff))
    
    interim_prob_eff[c2] = max(post_prob_efficacy_estimated_MTD)
    
    #We assess whether we will enroll a new cohort or not based on futility stopping rule
    stopping_futility_ind = ifelse(stopping_futility_ind == 0 & (interim_prob_eff[c2] >= delta_E0), 0, 1)
    sample_size_stopping_futility = ifelse(stopping_futility_ind == 0 & interim_prob_eff[c2] >= delta_E0 & (sample_size_stopping_futility < (NN_1 + NN_2)), sample_size_stopping_futility + m_cohort, sample_size_stopping_futility)
    
  }
  
  
  # ------------ Produce operating characteristics ------------- #
  
  #Pointwise MTD bias & Pointwise percentage of correct MTD recommendation:
  
  fdist<-function(x,MTD_point){
    
    # MTD_point is a point (x,y) on true MTD curve
    alpha0<-qnorm(post_rho00)
    alpha1<-qnorm(post_rho10) - qnorm(post_rho00)
    alpha2<-qnorm(post_rho01) - qnorm(post_rho00)
    
    d <- ((MTD_point[1]-x)^2+(MTD_point[2]-(qnorm(theta_Z)-alpha0-alpha1*x)/(alpha2+post_alpha3*x))^2)^0.5
    
    return(d)
  }
  
  true_mtdx_y0 = (qnorm(theta_Z) - qnorm(true_rho00)) / (qnorm(true_rho10) - qnorm(true_rho00))
  
  auxiliary_df = data.frame(x_true_MTD = seq(0,by = 0.01, min(1,true_mtdx_y0))) %>%
    mutate(y_true_MTD = twodimmtd2(true_rho00,true_rho01,true_rho10,true_alpha3,theta_Z, x_true_MTD)) %>%
    filter(y_true_MTD <= 1) %>%
    group_by(x_true_MTD) %>%
    mutate(minimum = optimize(f=fdist,interval = c(0,1), tol=0.0001,
                              MTD_point = c(x_true_MTD,y_true_MTD))$minimum,
           objective = optimize(f=fdist,interval = c(0,1), tol=0.0001,
                                MTD_point = c(x_true_MTD,y_true_MTD))$objective) %>%
    ungroup() %>%
    mutate(objective_bias = ifelse(y_true_MTD > twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta_Z,minimum), objective*(-1), objective)) 
  
  MTD_bias_df = auxiliary_df %>% dplyr::select(x_true_MTD, objective_bias) 
  
  percent_MTD_recomm_df = auxiliary_df %>% 
    dplyr::select(-c(objective_bias,minimum)) %>%
    mutate(tolerability_threshold_1 = Pointwise_percent_recomm_threshold_1*((x_true_MTD^2+y_true_MTD^2)^0.5),
           tolerability_threshold_2 = Pointwise_percent_recomm_threshold_2*((x_true_MTD^2+y_true_MTD^2)^0.5),
           correct_recomm_threshold_1 = ifelse(objective <= tolerability_threshold_1, 1, 0),
           correct_recomm_threshold_2 = ifelse(objective <= tolerability_threshold_2, 1, 0)) %>%
    dplyr::select(x_true_MTD, correct_recomm_threshold_1, correct_recomm_threshold_2)
  
  #Proportion of patients allocated in doses with true prob of efficacy > theta_E
  
  proportion_allocation_efficacious_doses_df = data.frame(x = X[(NN_1+1):(NN_1+NN_2)],
                                                                  y = Y[(NN_1+1):(NN_1+NN_2)]) %>%
    mutate(peff = peff(true_beta0, true_beta1, true_beta2, true_beta3, true_beta4, true_beta5, x, y),
           indicator = ifelse(peff > theta_E, 1, 0)) %>%
    summarise(mean_allocation = mean(indicator))
  
  
  post_prob_efficacious_estimated_MTD = data.frame(x = x_grid) %>%
    mutate(y = twodimmtd2(post_rho00, post_rho01, post_rho10, post_alpha3, theta_Z, x),
           post_prob_efficacious = post_prob_efficacy_estimated_MTD) %>%
    filter(y >= 0 & y <= 1) %>%
    mutate(obs_ind = 1:n()) %>%
    group_by(obs_ind) %>%
    mutate(true_peff = peff(true_beta0, true_beta1, true_beta2, true_beta3, true_beta4, true_beta5, x, y)) %>%
    ungroup()
  
  optimal_dose_combination = post_prob_efficacious_estimated_MTD %>%
    filter(post_prob_efficacious == max(post_prob_efficacious)) %>%
    dplyr::select(x,y,true_peff)
  
  #-----------------------------------------------------------------------------------------------#
  
  #Global output
  
  dt = list(Z = Z,
            E = E,
            X = X,
            Y = Y,
            proportion_allocation_efficacious_doses = proportion_allocation_efficacious_doses_df$mean_allocation,
            post_rho00 = post_rho00,
            post_rho10 = post_rho10,
            post_rho01 = post_rho01,
            post_alpha3 = post_alpha3,
            post_beta0 = post_beta0,
            post_beta1 = post_beta1,
            post_beta2 = post_beta2,
            post_beta3 = post_beta3,
            post_beta4 = post_beta4,
            post_beta5 = post_beta5,
            MTD_curve = twodimmtd2(true_rho00,true_rho01,true_rho10,true_alpha3,theta_Z, seq(0,1,length.out = 101)),
            interim_prob_dlt = interim_prob_dlt,
            interim_prob_eff = interim_prob_eff,
            max_post_prob_efficacy_higher_theta_E = max(post_prob_efficacy_estimated_MTD),
            absolute_difference_rho00 = post_rho00 - true_rho00,
            absolute_difference_rho10 = post_rho10 - true_rho10,
            absolute_difference_rho01 = post_rho01 - true_rho01,
            absolute_difference_alpha3 = post_alpha3 - true_alpha3,
            absolute_difference_beta0 = post_beta0 - true_beta0,
            absolute_difference_beta1 = post_beta1 - true_beta1,
            absolute_difference_beta2 = post_beta2 - true_beta2,
            absolute_difference_beta3 = post_beta3 - true_beta3,
            absolute_difference_beta4 = post_beta4 - true_beta4,
            absolute_difference_beta5 = post_beta5 - true_beta5,
            pointwise_MTD_bias = MTD_bias_df,
            percent_correct_MTD_recommentation = percent_MTD_recomm_df,
            sample_size_stopping_futility = sample_size_stopping_futility,
            sample_size_stopping_safety = sample_size_stopping_safety,
            optimal_dose_combination = optimal_dose_combination)
  
  dt
  
}








