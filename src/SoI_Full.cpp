// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List BSOIFull(vec Y, mat xi, mat Z_cov, int n_iter, vec beta_k0, vec sigma2_betak0, vec gamma_q0, vec sigma2_gammaq0,double a_delta0, double b_delta0,  double alpha, vec beta, vec gamma_cov, double sigma2_delta) {
	
	
	 // Initializing
	 int N_subj = Y.n_elem;
	 int K_x = xi.n_cols;
	 int q_cov = Z_cov.n_cols;
	 
	 
	 vec data_alpha(n_iter);
	 mat data_beta(n_iter,K_x);
	 mat data_gamma_cov(n_iter,q_cov);
	 vec data_sigma2_delta(n_iter);
	 
	 
	 
	 double sigma_beta_temp = 0;
	 double mean_beta_temp = 0;
	 double sigma_gamma_temp = 0;
	 double mean_gamma_temp = 0;
	 vec xi_sum_square = zeros<vec>(K_x);
	 vec XCov_sum_square = zeros<vec>(q_cov);
	 int current_posi;
	 uvec xi_col_indices(K_x - 1);
	 uvec gamma_col_indices(q_cov - 1);
	 
	 
	 for (int kk =0; kk < K_x; kk++){
		xi_sum_square[kk] = sum(pow(norm(xi.col(kk)),2));
	    }
		
	 for (int kk =0; kk < q_cov; kk++){
		XCov_sum_square[kk] = sum(pow(norm(Z_cov.col(kk)),2));
	    }
		
	 
	 for (int it = 0; it < n_iter; it++){
		 
		 alpha = sum(Y - xi*beta -Z_cov*gamma_cov)/N_subj + rnorm(1,0,1)[0]*sqrt(sigma2_delta/N_subj);
		 
		 for (int kk =0; kk < K_x; kk++){
			 current_posi = 0;
			 sigma_beta_temp = 1/(xi_sum_square[kk]/sigma2_delta + 1/sigma2_betak0[kk]);
			 if(K_x ==1){
				mean_beta_temp = sigma_beta_temp*(beta_k0[kk]/sigma2_betak0[kk] + sum((Y - alpha - Z_cov*gamma_cov)%xi.col(kk))/sigma2_delta);
			 }else{
				 for(int kk2 = 0; kk2 < K_x; kk2++){
					 if(kk2!=kk){
						xi_col_indices[current_posi] = kk2;
						current_posi++;
					 }
				 }
				 mean_beta_temp = sigma_beta_temp*(beta_k0[kk]/sigma2_betak0[kk] + sum((Y - alpha - Z_cov*gamma_cov - xi.cols(xi_col_indices)*beta(xi_col_indices))%xi.col(kk))/sigma2_delta);
			 }
			 beta[kk] = mean_beta_temp + sqrt(sigma_beta_temp)*rnorm(1,0,1)[0];
		 }
		 
		 for (int kk = 0; kk < q_cov; kk++){
			 current_posi = 0;
			 sigma_gamma_temp = 1/(XCov_sum_square[kk]/sigma2_delta + 1/sigma2_gammaq0[kk]);
			  if(q_cov ==1){
				mean_gamma_temp = sigma_gamma_temp*(gamma_q0[kk]/sigma2_gammaq0[kk] + sum((Y - alpha - xi*beta)%Z_cov.col(kk))/sigma2_delta);
			 }else{
				 for(int kk2 = 0; kk2 < q_cov; kk2++){
					 if(kk2!=kk){
						gamma_col_indices[current_posi] = kk2;
						current_posi++;
					 }
				 }
				 mean_gamma_temp = sigma_gamma_temp*(gamma_q0[kk]/sigma2_gammaq0[kk] + sum((Y - alpha - xi*beta - Z_cov.cols(gamma_col_indices)*gamma_cov(gamma_col_indices))%Z_cov.col(kk))/sigma2_delta);
			 }
			 gamma_cov[kk] = mean_gamma_temp + sqrt(sigma_gamma_temp)*rnorm(1,0,1)[0];
			 
			
		 }
		 sigma2_delta = 1 / rgamma(1, a_delta0 + 0.5 * N_subj, 1 / (b_delta0 + 0.5 *
               sum(pow(norm(Y - alpha - xi*beta -Z_cov*gamma_cov),2))))[0];
		 
		data_alpha[it] = alpha;
		data_beta.row(it) = beta.t();
		data_gamma_cov.row(it) = gamma_cov.t();
		data_sigma2_delta[it] = sigma2_delta;
 
		if (it > 0 && it % (n_iter/100) == 0)
			Rprintf("\rRunning MCMC, %d%% completed...", it / (n_iter/100));
	 }		  
			  
	 return List::create(Named("alpha") = data_alpha,
                         Named("beta") = data_beta,
                         Named("sigma2_delta") = data_sigma2_delta,
                         Named("gamma_cov") = data_gamma_cov);	  

}

