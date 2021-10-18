#'Design parameters


# Target Probability of DLT 
theta_Z = 1/3 

#Probability of standard treatment of care
theta_E = 0.15

# Number of patients in a trial 
NN_1 = 30
NN_2 = 30
m_cohort = 5

# Controls size of jump
delta1 = 0.2
Xlow = 0
Ylow = 0

#Early stopping thresholds
delta_Z1 = 0.5
delta_Z2 = 0.7
delta_E0 = 0.1
delta_E1 = 0.95

# Minimum and maximum doses available in the trial
# X dose is for Cabazitaxel
Xmin = 10  
Xmax = 25
# Y dose is for Cisplatin
Ymin = 50
Ymax = 100

lbx  =  (0 - Xmin)/(Xmax - Xmin)
lby  =  (0 - Ymin)/(Ymax - Ymin)

qq = 101
x_grid = seq(0,1,length.out = qq)

# Feasibility bound (EWOC)
alpha = numeric()
alphaa = seq(0.4,0.5,by=0.05)
for (ind in 1:3){
  alpha[ind] = alphaa[ind]
}
for (ind in 4:(NN_1-1)){
  alpha[ind] = alphaa[3]
}


# mcmc parameters
chains = 1
burn = 5000
mm = 2500


#Pointwise percent correct recommentation tolerability thresholds (operating characteristics)

Pointwise_percent_recomm_threshold_1 = 0.1
Pointwise_percent_recomm_threshold_2 = 0.2

#Dose-toxicity and dose-efficacy model parameters

#Dose-toxicity scenario 1

true_rho00_DLT_sc1 = 1e-7
true_rho01_DLT_sc1 = 0.3
true_rho10_DLT_sc1 = 0.3
true_alpha3_DLT_sc1  = 2
true_corr_DLT_sc1 = 0

#Dose-toxicity scenario 2

true_rho00_DLT_sc2 = 1e-5
true_rho01_DLT_sc2 = 0.01
true_rho10_DLT_sc2 = 0.005
true_alpha3_DLT_sc2 = 9
true_corr_DLT_sc2 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 1, H1

true_beta0_DLT_sc1_EFF_sc1_H1 = -5.51
true_beta1_DLT_sc1_EFF_sc1_H1 = 2
true_beta2_DLT_sc1_EFF_sc1_H1 = 4.3
true_beta3_DLT_sc1_EFF_sc1_H1 = 10
true_beta4_DLT_sc1_EFF_sc1_H1 = 0
true_beta5_DLT_sc1_EFF_sc1_H1 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 2, H1

true_beta0_DLT_sc1_EFF_sc2_H1 = -5.51
true_beta1_DLT_sc1_EFF_sc2_H1 = 4.3
true_beta2_DLT_sc1_EFF_sc2_H1 = 2
true_beta3_DLT_sc1_EFF_sc2_H1 = 10
true_beta4_DLT_sc1_EFF_sc2_H1 = 0
true_beta5_DLT_sc1_EFF_sc2_H1 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 3, H1

true_beta0_DLT_sc1_EFF_sc3_H1 = -6.5
true_beta1_DLT_sc1_EFF_sc3_H1 = 6.17
true_beta2_DLT_sc1_EFF_sc3_H1 = 5.5
true_beta3_DLT_sc1_EFF_sc3_H1 = 0
true_beta4_DLT_sc1_EFF_sc3_H1 = 0
true_beta5_DLT_sc1_EFF_sc3_H1 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 4, H1

true_beta0_DLT_sc1_EFF_sc4_H1 = -4
true_beta1_DLT_sc1_EFF_sc4_H1 = 1.25
true_beta2_DLT_sc1_EFF_sc4_H1 = 1.25
true_beta3_DLT_sc1_EFF_sc4_H1 = 12
true_beta4_DLT_sc1_EFF_sc4_H1 = 0
true_beta5_DLT_sc1_EFF_sc4_H1 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 1, H1

true_beta0_DLT_sc2_EFF_sc1_H1 = -2
true_beta1_DLT_sc2_EFF_sc1_H1 = 0.05
true_beta2_DLT_sc2_EFF_sc1_H1 = 1.57
true_beta3_DLT_sc2_EFF_sc1_H1 = 1
true_beta4_DLT_sc2_EFF_sc1_H1 = 0
true_beta5_DLT_sc2_EFF_sc1_H1 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 2, H1

true_beta0_DLT_sc2_EFF_sc2_H1 = -2
true_beta1_DLT_sc2_EFF_sc2_H1 = 1.55
true_beta2_DLT_sc2_EFF_sc2_H1 = 0.05
true_beta3_DLT_sc2_EFF_sc2_H1 = 1
true_beta4_DLT_sc2_EFF_sc2_H1 = 0
true_beta5_DLT_sc2_EFF_sc2_H1 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 3, H1

true_beta0_DLT_sc2_EFF_sc3_H1 = -5.8
true_beta1_DLT_sc2_EFF_sc3_H1 = 4.63
true_beta2_DLT_sc2_EFF_sc3_H1 = 4.73
true_beta3_DLT_sc2_EFF_sc3_H1 = 0
true_beta4_DLT_sc2_EFF_sc3_H1 = 0
true_beta5_DLT_sc2_EFF_sc3_H1 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 4, H1

true_beta0_DLT_sc2_EFF_sc4_H1 = -6.49
true_beta1_DLT_sc2_EFF_sc4_H1 = 0.2
true_beta2_DLT_sc2_EFF_sc4_H1 = 0.2
true_beta3_DLT_sc2_EFF_sc4_H1 = 26
true_beta4_DLT_sc2_EFF_sc4_H1 = 0
true_beta5_DLT_sc2_EFF_sc4_H1 = 0


#Dose-toxicity scenario 1, dose-efficacy scenario 1, H0

true_beta0_DLT_sc1_EFF_sc1_H0 = -6.3
true_beta1_DLT_sc1_EFF_sc1_H0 = 2
true_beta2_DLT_sc1_EFF_sc1_H0 = 4.3
true_beta3_DLT_sc1_EFF_sc1_H0 = 10
true_beta4_DLT_sc1_EFF_sc1_H0 = 0
true_beta5_DLT_sc1_EFF_sc1_H0 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 2, H0

true_beta0_DLT_sc1_EFF_sc2_H0 = -6.3
true_beta1_DLT_sc1_EFF_sc2_H0 = 4.3
true_beta2_DLT_sc1_EFF_sc2_H0 = 2
true_beta3_DLT_sc1_EFF_sc2_H0 = 10
true_beta4_DLT_sc1_EFF_sc2_H0 = 0
true_beta5_DLT_sc1_EFF_sc2_H0 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 3, H0

true_beta0_DLT_sc1_EFF_sc3_H0 = -7.3
true_beta1_DLT_sc1_EFF_sc3_H0 = 6.17
true_beta2_DLT_sc1_EFF_sc3_H0 = 5.5
true_beta3_DLT_sc1_EFF_sc3_H0 = 0
true_beta4_DLT_sc1_EFF_sc3_H0 = 0
true_beta5_DLT_sc1_EFF_sc3_H0 = 0

#Dose-toxicity scenario 1, dose-efficacy scenario 4, H0

true_beta0_DLT_sc1_EFF_sc4_H0 = -4.8
true_beta1_DLT_sc1_EFF_sc4_H0 = 1.25
true_beta2_DLT_sc1_EFF_sc4_H0 = 1.25
true_beta3_DLT_sc1_EFF_sc4_H0 = 12
true_beta4_DLT_sc1_EFF_sc4_H0 = 0
true_beta5_DLT_sc1_EFF_sc4_H0 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 1, H0

true_beta0_DLT_sc2_EFF_sc1_H0 = -2.8
true_beta1_DLT_sc2_EFF_sc1_H0 = 0.05
true_beta2_DLT_sc2_EFF_sc1_H0 = 1.57
true_beta3_DLT_sc2_EFF_sc1_H0 = 1
true_beta4_DLT_sc2_EFF_sc1_H0 = 0
true_beta5_DLT_sc2_EFF_sc1_H0 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 2, H0

true_beta0_DLT_sc2_EFF_sc2_H0 = -2.8
true_beta1_DLT_sc2_EFF_sc2_H0 = 1.55
true_beta2_DLT_sc2_EFF_sc2_H0 = 0.05
true_beta3_DLT_sc2_EFF_sc2_H0 = 1
true_beta4_DLT_sc2_EFF_sc2_H0 = 0
true_beta5_DLT_sc2_EFF_sc2_H0 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 3, H0

true_beta0_DLT_sc2_EFF_sc3_H0 = -6.6
true_beta1_DLT_sc2_EFF_sc3_H0 = 4.63
true_beta2_DLT_sc2_EFF_sc3_H0 = 4.73
true_beta3_DLT_sc2_EFF_sc3_H0 = 0
true_beta4_DLT_sc2_EFF_sc3_H0 = 0
true_beta5_DLT_sc2_EFF_sc3_H0 = 0

#Dose-toxicity scenario 2, dose-efficacy scenario 4, H0

true_beta0_DLT_sc2_EFF_sc4_H0 = -7.28
true_beta1_DLT_sc2_EFF_sc4_H0 = 0.2
true_beta2_DLT_sc2_EFF_sc4_H0 = 0.2
true_beta3_DLT_sc2_EFF_sc4_H0 = 26
true_beta4_DLT_sc2_EFF_sc4_H0 = 0
true_beta5_DLT_sc2_EFF_sc4_H0 = 0
  