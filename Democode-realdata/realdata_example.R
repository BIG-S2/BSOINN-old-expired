set.seed(as.numeric(Sys.time()))
library(BSOINN)  ##load the library

####file names
Y_file = "Y.txt"
Cov_file = "COV.txt"  
R_file = "R.txt"
Xi_file = "XI.txt"

####define constants
N_REPS = 1
N_subj=802     ##sample size
n_iter=20000    ##total iterations
n_burnin=10000  ##burn-in phase
K_x = 5        ##number of eigenimages retained   
q_cov = 7      ##dimension of covariates in SoI
q_cov_NN = 6   ##dimension of covariates in missing mechanism after removing the instrument 
ind_cov_NN = c(T,T,F,T,T,T,T)  ##we define educational level (the third) as an instrument


####define objects to record the estimation results
Ealpha_IN = array(0, dim = c(N_REPS,1))
Egamma_IN = array(0,dim = c(N_REPS,q_cov))
Ebeta_IN = array(0, dim = c(N_REPS,K_x))
Edelta_IN = array(0, dim = c(N_REPS,1))

SDalpha_IN = array(0, dim = c(N_REPS,1))
SDgamma_IN = array(0,dim = c(N_REPS,q_cov))
SDbeta_IN = array(0, dim = c(N_REPS,K_x))
SDdelta_IN = array(0, dim = c(N_REPS,1))

Ealpha_NN = array(0, dim = c(N_REPS,1))
Egamma_NN = array(0,dim = c(N_REPS,q_cov))
Ebeta_NN = array(0, dim = c(N_REPS,K_x))
Edelta_NN = array(0, dim = c(N_REPS,1))
Ealpha_r_NN = array(0, dim = c(N_REPS,1))
Egamma_r_NN = array(0,dim = c(N_REPS,q_cov_NN))
Ebeta_r_NN =array(0, dim = c(N_REPS,K_x))
Ephi_NN = array(0, dim = c(N_REPS,1))

SDalpha_NN = array(0, dim = c(N_REPS,1))
SDgamma_NN = array(0,dim = c(N_REPS,q_cov))
SDbeta_NN = array(0, dim = c(N_REPS,K_x))
SDdelta_NN = array(0, dim = c(N_REPS,1))
SDalpha_r_NN = array(0, dim = c(N_REPS,1))
SDgamma_r_NN = array(0,dim = c(N_REPS,q_cov_NN))
SDbeta_r_NN =array(0, dim = c(N_REPS,K_x))
SDphi_NN = array(0, dim = c(N_REPS,1))

varphi_acceptrate_NN = array(0, dim = c(N_REPS,1))
ym_acceptrate_NN = array(0,dim=c(N_REPS, 1))

####analyze the real data
for(CIR in 1:N_REPS){
	####read the values
	Y = matrix(scan(Y_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1, byrow = T)   ##outcome subject to missingness
	X_cov = matrix(scan(Cov_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=q_cov, byrow = T)  ##covariates in SoI regression
	R = matrix(scan(R_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1,byrow = T)  ##missing indicator
	xi = matrix(scan(Xi_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=K_x, byrow=T)  ##eigenscores from FPCA, and the FPCA could be easily conducted through matlab or the fast.svd function in R
	cat("\n")
	print(paste("Replication:", CIR))
	
	X_cov_NN = X_cov[,ind_cov_NN]
	N_mis = sum(R==1)
	
	####SoI_IN procedure
	print("BSOI_IN:")
	##hyperparameters
	beta_k0 = rep(0,K_x)
	sigma2_betak0 = rep(100,K_x)
	gamma_q0 = rep(0.01,q_cov)
	sigma2_gammaq0 = rep(100,q_cov)
	a_delta0 = 9
	b_delta0 = 3
    
	##initial values for the patameters
	alpha = 0.1
	beta = rep(0.1,K_x)
	gamma_cov = rep(0.1,q_cov)
	sigma2_delta = 1
	Y[R==1] = rnorm(N_mis,mean = 0, sd = 1)
	
	##conduct the Bayesian analysis using the function SoIIN()
	##Only the data inputs are compulsory, and we may implement the analysis as
    ##model_IN <- BSOIIN(Y, R, xi, X_cov, n_iter)
	model_IN <- BSOIIN(Y, R, xi, X_cov, n_iter,  beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha, beta,  gamma_cov, sigma2_delta)
	
	##record the results
	Ealpha_IN[CIR,] = mean(model_IN$alpha[(n_burnin+1):n_iter])
	Ebeta_IN[CIR,] = apply(model_IN$beta[(n_burnin+1):n_iter,],2,mean)
	Egamma_IN[CIR,] = apply(model_IN$gamma_cov[(n_burnin+1):n_iter,],2,mean)
	Edelta_IN[CIR,] = mean(model_IN$sigma2_delta[(n_burnin+1):n_iter])
	
	SDalpha_IN[CIR,] = sd(model_IN$alpha[(n_burnin+1):n_iter])
	SDbeta_IN[CIR,] = apply(model_IN$beta[(n_burnin+1):n_iter,],2,sd)
	SDgamma_IN[CIR,] = apply(model_IN$gamma_cov[(n_burnin+1):n_iter,],2,sd)
	SDdelta_IN[CIR,] = sd(model_IN$sigma2_delta[(n_burnin+1):n_iter])
	
	####SoI_NN procedure
	cat("\n")
	print("BSOI_NN:")
	##variance parameters that control acceptance rate in MH algorithm
	var_ym = 20
	var_varphi = 1.5
	
	##hyperparameters
	beta_k0 = rep(0,K_x)
	sigma2_betak0 = rep(100,K_x)
	gamma_q0 = rep(0.01,q_cov)
	sigma2_gammaq0 = rep(100,q_cov)
	a_delta0 = 9
	b_delta0 = 3
	varphi_0 = c(0, rep(0,q_cov_NN), rep(0,K_x), 0)
	sigma_varphi0 = diag(2+q_cov_NN+K_x)
	
	##initial values for the parameters
	alpha = 0.1
	beta = rep(0.1,K_x)
	gamma_cov = rep(0.1,q_cov)
	sigma2_delta = 1
	alpha_r = 0.1
	beta_r = rep(0.1,K_x)
	gamma_cov_r = rep(0.1,q_cov_NN)
	phi = 0.1
	Y[R==1] = rnorm(N_mis,mean = 0, sd = 1)
	
	##conduct the Bayesian analysis using the function SoINN()
	##Only the data inputs are compulsory, and we may implement the analysis as
    ##model_NN <- BSOINN(Y, R, xi, X_cov,  X_cov_NN, n_iter)
	model_NN <- BSOINN(Y, R, xi, X_cov, X_cov_NN, n_iter, beta_k0, sigma2_betak0, gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, varphi_0, sigma_varphi0,  var_ym,var_varphi, alpha, beta, gamma_cov, sigma2_delta,  alpha_r, beta_r, gamma_cov_r, phi)
	
	##record the results
	Ealpha_NN[CIR,] = mean(model_NN$alpha[(n_burnin+1):n_iter])
	Ebeta_NN[CIR,] = apply(model_NN$beta[(n_burnin+1):n_iter,],2,mean)
	Egamma_NN[CIR,] = apply(model_NN$gamma_cov[(n_burnin+1):n_iter,],2,mean)
	Edelta_NN[CIR,] = mean(model_NN$sigma2_delta[(n_burnin+1):n_iter])
	
	Ealpha_r_NN[CIR,] = mean(model_NN$alpha_r[(n_burnin+1):n_iter])
	Ebeta_r_NN[CIR,] = apply(model_NN$beta_r[(n_burnin+1):n_iter,],2,mean)
	Egamma_r_NN[CIR,] = apply(model_NN$gamma_cov_r[(n_burnin+1):n_iter,],2,mean)
	Ephi_NN[CIR,] = mean(model_NN$phi[(n_burnin+1):n_iter])
	
	SDalpha_NN[CIR,] = sd(model_NN$alpha[(n_burnin+1):n_iter])
	SDbeta_NN[CIR,] = apply(model_NN$beta[(n_burnin+1):n_iter,],2,sd)
	SDgamma_NN[CIR,] = apply(model_NN$gamma_cov[(n_burnin+1):n_iter,],2,sd)
	SDdelta_NN[CIR,] = sd(model_NN$sigma2_delta[(n_burnin+1):n_iter])
	
	SDalpha_r_NN[CIR,] = sd(model_NN$alpha_r[(n_burnin+1):n_iter])
	SDbeta_r_NN[CIR,] = apply(model_NN$beta_r[(n_burnin+1):n_iter,],2,sd)
	SDgamma_r_NN[CIR,] = apply(model_NN$gamma_cov_r[(n_burnin+1):n_iter,],2,sd)
	SDphi_NN[CIR,] = sd(model_NN$phi[(n_burnin+1):n_iter])
	
	varphi_acceptrate_NN[CIR,] = (model_NN$acce_varphi_n)/n_iter
	ym_acceptrate_NN[CIR,] =  sum(model_NN$acce_ym)/n_iter/N_mis
} 